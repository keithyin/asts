use lazy_static::lazy_static;
use std::{
    collections::HashMap,
    sync::{Arc, Mutex},
    thread,
    time::Instant,
};

use crate::{
    align_sbr_to_smc, params, read_info::ReadInfo, reporter::Reporter,
    sbr_and_ref_to_cs::build_msa_result_from_records,
};
use crate::{sbr_and_ref_to_cs::MsaResult, SubreadsAndSmc};
use crossbeam::channel::Sender;
use mm2::gskits::{phreq::phreq_list_2_quality, utils::ScopedTimer};

use tracing;

pub fn align_sbr_and_fake_cs_to_cs_worker(
    recv: crossbeam::channel::Receiver<SubreadsAndSmc>,
    sender: Sender<MsaResult>,
    target_idx: &HashMap<String, (usize, usize)>,
    map_params: &params::MapParams,
    align_params: &params::AlignParams,
    oup_params: &params::OupParams,
    reporter: Arc<Mutex<Reporter>>,
    low_q: Option<u8>,
) {
    let mut scoped_timer = ScopedTimer::new();

    let mut max_time = 0;
    let mut max_time_qname = String::new();

    let mut inp_sbrs = 0;
    let mut out_sbrs = 0;
    let mut out_smc = 0;
    let mut fallback_num = 0;
    let mut fallback_rescued_num = 0;

    for mut subreads_and_smc in recv {
        // tracing::info!("sbr_cnt:{}-{}", subreads_and_smc.smc.name, subreads_and_smc.subreads.len());
        let timer = scoped_timer.perform_timing();
        let start = Instant::now();

        let pred_q = phreq_list_2_quality(subreads_and_smc.smc.qual.as_ref().unwrap());

        // let fake_cs = generate_fake_cs_from_cs(&subreads_and_smc.smc);

        let ref_read_info =
            ReadInfo::new_fa_record("ref/-1".to_string(), subreads_and_smc.smc.seq.clone());
        subreads_and_smc.subreads.push(ref_read_info);
        subreads_and_smc.subreads.iter_mut().for_each(|read_info| {
            read_info.dw = None;
            read_info.ar = None;
            read_info.cr = None;
            read_info.nn = None;
        });

        inp_sbrs += subreads_and_smc.subreads.len();
        let (mut align_infos, no_hit_indices) = align_sbr_to_smc(
            &subreads_and_smc,
            (0..subreads_and_smc.subreads.len()).collect(),
            target_idx,
            map_params,
            align_params,
            oup_params,
            false,
        );

        fallback_num += no_hit_indices.len();
        if align_infos.len() < 20 && !no_hit_indices.is_empty() {
            tracing::warn!(
                "fallback: smc_name:{}, num_sbrs:{}, num_fallback_sbrs:{}",
                subreads_and_smc.smc.name,
                subreads_and_smc.subreads.len(),
                no_hit_indices.len(),
            );
            let (fallback_align_infos, _) = align_sbr_to_smc(
                &subreads_and_smc,
                no_hit_indices,
                target_idx,
                map_params,
                align_params,
                oup_params,
                true,
            );

            fallback_rescued_num += fallback_align_infos.len();

            align_infos.extend(fallback_align_infos);
        }

        let elapsed_secs = start.elapsed().as_secs();
        if elapsed_secs > max_time {
            max_time = elapsed_secs;
            max_time_qname = subreads_and_smc.smc.name.clone();
        }

        out_sbrs += align_infos.len();

        let align_res = build_msa_result_from_records(
            align_infos,
            &subreads_and_smc.smc.seq,
            &subreads_and_smc.smc.name,
            subreads_and_smc.smc.qual.as_deref(),
            low_q,
        );
        if align_res.is_none() {
            continue;
        }
        let mut align_res = align_res.unwrap();
        align_res.identity = pred_q;
        align_res = align_res.extract_lowq_region();

        timer.done_with_cnt(1);
        out_smc += 1;
        sender.send(align_res).unwrap();
    }

    let filter_by_alignment = inp_sbrs - out_sbrs;
    {
        let mut reporter_ = reporter.lock().unwrap();
        reporter_.sbr_reporter.filter_by_alignment += filter_by_alignment;
        reporter_.channel_reporter.out_num += out_smc;
        reporter_.sbr_reporter.fallback_num += fallback_num;
        reporter_.sbr_reporter.fallback_resuced_num += fallback_rescued_num;
    }

    tracing::info!(
        "align_sbr_to_smc_worker-{:?}. speed:{:.4}/s. single_smc_max_time:{}secs, max_time_qname:{}",
        thread::current().id(),
        scoped_timer.speed(Some(1000_000_000)),
        max_time,
        max_time_qname
    );
}

lazy_static! {
    static ref LOWER20_BASE_MAPPING: HashMap<u8, u8> = {
        let mut m = HashMap::new();
        m.insert(b'A', b'T');
        m.insert(b'T', b'A');
        m.insert(b'C', b'G');
        m.insert(b'G', b'C');
        m
    };
    static ref LOWER10_BASE_MAPPING: HashMap<u8, u8> = {
        let mut m = HashMap::new();
        m.insert(b'A', b'C');
        m.insert(b'T', b'C');
        m.insert(b'C', b'A');
        m.insert(b'G', b'A');
        m
    };
    static ref LOWER5_BASE_MAPPING: HashMap<u8, u8> = {
        let mut m = HashMap::new();
        m.insert(b'A', b'G');
        m.insert(b'T', b'G');
        m.insert(b'C', b'T');
        m.insert(b'G', b'T');
        m
    };
}
