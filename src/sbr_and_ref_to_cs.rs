use std::{
    collections::{HashMap, HashSet},
    sync::{Arc, Mutex},
    thread,
    time::Instant,
};

use crate::SubreadsAndSmc;
use crate::{align_sbr_to_smc, reporter::Reporter};
use crossbeam::channel::Sender;
use minimap2::{Aligner, Built, PresetSet};
use mm2::gskits::{
    self,
    ds::ReadInfo,
    gsbam::{
        bam_record_ext::BamRecordExt, plp_counts_from_records::compute_max_ins_of_each_ref_position,
    },
    utils::ScopedTimer,
};
use mm2::{
    mapping_ext::MappingExt,
    params::{InputFilterParams, TOverrideAlignerParam},
    NoMemLeakAligner,
};
use rust_htslib::bam::{ext::BamRecordExtensions, Read, Record};
use tracing;

use serde::{Serialize, Deserialize};

fn align_cs_to_ref(cs: &str, ref_aligner: &Aligner<Built>, ref_seq: &str) -> Option<String> {
    let hits = ref_aligner
        .map(
            cs.as_bytes(),
            false,
            false,
            None,
            Some(&[68719476736]), // 67108864 eqx, 68719476736 secondary
            None,
        )
        .unwrap();

    return if hits.len() == 1 {
        let hit = &hits[0];
        let target_start = if hit.target_start > 10 {
            hit.target_start - 10
        } else {
            0
        } as usize;

        let target_end = (hit.target_end + 10).min(ref_seq.len() as i32) as usize;

        Some(ref_seq[target_start..target_end].to_string())
    } else {
        None
    };
}

#[derive(Debug, Serialize)]
pub struct MsaResult {
    pub identity: f32,
    pub mm: u32,
    pub ins: u32,
    pub del: u32,
    pub msa_seqs: Vec<String>,
    pub names: Vec<String>,
}

pub fn align_sbr_and_ref_to_cs_worker(
    recv: crossbeam::channel::Receiver<SubreadsAndSmc>,
    sender: Sender<MsaResult>,
    ref_aligner: &Aligner<Built>,
    ref_seq: &str,
    target_idx: &HashMap<String, (usize, usize)>,
    map_params: &mm2::params::MapParams,
    align_params: &mm2::params::AlignParams,
    oup_params: &mm2::params::OupParams,
    reporter: Arc<Mutex<Reporter>>,
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

        let ref_seq = align_cs_to_ref(&subreads_and_smc.smc.seq, ref_aligner, ref_seq);
        if ref_seq.is_none() {
            continue;
        }
        let ref_read_info = ReadInfo::new_fa_record("ref/-1".to_string(), ref_seq.unwrap());
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
        );
        if align_res.is_none() {
            continue;
        }
        let align_res = align_res.unwrap();

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

fn build_msa_result_from_records(
    mut records: Vec<Record>,
    ref_seq: &str,
    ref_name: &str,
) -> Option<MsaResult> {
    if records.len() == 0 {
        return None;
    }
    let ref_seq = ref_seq.as_bytes();
    records.sort_by_key(|record| {
        BamRecordExt::new(record)
            .get_qname()
            .rsplit_once("/")
            .unwrap()
            .1
            .parse::<i32>()
            .unwrap()
    });

    let first_record = records.first().unwrap();
    let first_record = BamRecordExt::new(first_record);
    if !first_record.get_qname().starts_with("ref") {
        return None;
    }

    let major_pos_ins = compute_max_ins_of_each_ref_position(&records, None, None, None);
    let mut major_pos_ins_vec = major_pos_ins
        .iter()
        .map(|(&k, &v)| (k as usize, v as usize))
        .collect::<Vec<_>>();
    major_pos_ins_vec.sort_by_key(|v| v.0);
    let tot_len = major_pos_ins_vec.iter().map(|v| v.1 + 1).sum::<usize>();

    let mut cursor = 0;
    let major_start_point = major_pos_ins_vec
        .iter()
        .map(|&(_, max_ins)| {
            let cur_point = cursor;
            cursor += max_ins + 1;
            cur_point
        })
        .collect::<Vec<_>>();

    let major_pos2major_starting_point = major_pos_ins_vec
        .iter()
        .map(|&(ma, _)| ma)
        .zip(major_start_point.into_iter())
        .collect::<HashMap<_, _>>();

    let mut anchor_aligned = vec!['-' as u8; tot_len];
    major_pos_ins_vec.iter().for_each(|(major, _)| {
        anchor_aligned[*major_pos2major_starting_point.get(major).unwrap()] = ref_seq[*major];
    });

    let anchor_aligned = String::from_utf8(anchor_aligned).unwrap();

    let msa = records
        .iter()
        .map(|record| {
            let mut aligned = vec!['-' as u8; tot_len];
            build_one_record_of_msa(record, &major_pos2major_starting_point, &mut aligned);
            String::from_utf8(aligned).unwrap()
        })
        .collect::<Vec<_>>();

    let mut msa_seqs = vec![anchor_aligned];
    msa_seqs.extend(msa.into_iter());
    let mut names = vec![ref_name.to_string()];
    names.extend(
        records
            .iter()
            .map(|record| BamRecordExt::new(record).get_qname()),
    );

    Some(MsaResult {
        identity: 0.0,
        mm: 0,
        ins: 0,
        del: 0,
        msa_seqs: msa_seqs,
        names: names,
    })
}

fn build_one_record_of_msa(
    record: &Record,
    major_pos2major_starting_point: &HashMap<usize, usize>,
    result: &mut Vec<u8>,
) {
    let record_ext = BamRecordExt::new(record);
    let ref_start = record_ext.reference_start() as i64;
    let ref_end = record_ext.reference_end() as i64;
    let q_start = record_ext.query_alignment_start() as i64;
    let q_end = record_ext.query_alignment_end() as i64;
    let query = record_ext.get_seq();
    let query = query.as_bytes();
    let mut q_pos_cursor = None;
    let mut r_pos_cursor = None;
    let mut delta = 0;
    for [qpos, rpos] in record.aligned_pairs_full() {
        if qpos.is_some() {
            q_pos_cursor = qpos;
        }
        if rpos.is_some() {
            r_pos_cursor = rpos;
        }

        if q_pos_cursor.is_none() || r_pos_cursor.is_none() {
            continue;
        }

        if q_pos_cursor.unwrap() < q_start || r_pos_cursor.unwrap() < ref_start {
            continue;
        }

        if q_pos_cursor.unwrap() >= q_end || r_pos_cursor.unwrap() >= ref_end {
            break;
        }

        if rpos.is_some() {
            delta = 0;
        } else {
            delta += 1;
        }
        let r_cursor = r_pos_cursor.map(|v| v as usize).unwrap();
        let r_cursor = &r_cursor;
        if !major_pos2major_starting_point.contains_key(r_cursor) {
            continue;
        }

        let base_pos = *major_pos2major_starting_point.get(r_cursor).unwrap();

        if let Some(qpos) = qpos.map(|v| v as usize) {
            result[base_pos + delta] = query[qpos];
        }
    }
}
