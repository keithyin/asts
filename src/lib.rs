use std::{
    collections::{HashMap, HashSet},
    env,
    hash::Hash,
    sync::{Arc, Mutex},
    thread,
    time::Instant,
};

use ::gskits::gsbam::{bam_header_ext::BamHeaderExt, bam_record_ext::BamRecordExt};
use crossbeam::channel::Sender;
use gskits::{self, utils::ScopedTimer};
use params::{AlignParams, InputFilterParams, MapParams, OupParams, TOverrideAlignSetup};
use read_info::ReadInfo;
use reporter::Reporter;
use rust_htslib::bam::{ext::BamRecordExtensions, Read};
use tracing;
use tracing_subscriber::EnvFilter;

use crate::align::{
    build_aligner, cigar, cigar_adjust::cigar_adjust_poly_gap_left_align, AlignSetup, RefAligner,
};

pub mod align;
pub mod bam_writer;
pub mod mapping2record;
pub mod params;
pub mod read_info;
pub mod reporter;
pub mod sbr_and_cs_to_cs;
pub mod sbr_and_ref_to_cs;
pub mod utils;

// cCsSiIf int8, uint8, int16, uint16, int32, uint32, float
type BamRecord = rust_htslib::bam::record::Record;

/// 一批比对结果（每个 SMC channel 一次）。
pub struct AlignResult {
    pub records: Vec<BamRecord>,
}

pub struct SubreadsAndSmc {
    pub smc: ReadInfo,
    smc_len: usize,
    pub subreads: Vec<ReadInfo>,
}

impl SubreadsAndSmc {
    pub fn new(smc_record: &rust_htslib::bam::Record) -> Self {
        let smc_len = smc_record.seq_len();
        Self {
            smc: ReadInfo::from_bam_record(smc_record, None, &HashSet::new()),
            smc_len: smc_len,
            subreads: vec![],
        }
    }

    pub fn from_seq_qual(seq: &str, qual: Vec<u8>) -> Self {
        let smc_len = seq.len();
        Self {
            smc: ReadInfo::new_fq_record(format!("smc"), seq.to_string(), qual),
            smc_len: smc_len,
            subreads: vec![],
        }
    }

    pub fn add_subread(
        &mut self,
        record: &rust_htslib::bam::Record,
        tags: &HashSet<String>,
    ) -> bool {
        let sbr_len = record.seq_len() as f64;
        let smc_len = self.smc_len as f64;
        let max_len = smc_len * 3.0;
        let min_len = smc_len * 0.0;

        if sbr_len > min_len && sbr_len < max_len {
            self.subreads
                .push(ReadInfo::from_bam_record(record, None, tags));
            true
        } else {
            false
        }
    }

    pub fn add_subread_str(&mut self, name: String, seq: String) -> bool {
        let sbr_len = seq.len() as f64;
        let smc_len = self.smc_len as f64;
        let max_len = smc_len * 3.0;
        let min_len = smc_len * 0.0;

        if sbr_len > min_len && sbr_len < max_len {
            self.subreads.push(ReadInfo::new_fa_record(name, seq));
            true
        } else {
            false
        }
    }
}

pub struct SingleChannelAlignRes {
    pub records: Vec<BamRecord>,
}

pub fn subreads_and_smc_generator(
    sorted_sbr_bam: &str,
    sorted_smc_bam: &str,
    input_filter_params: &InputFilterParams,
    oup_params: &OupParams,
    sender: crossbeam::channel::Sender<SubreadsAndSmc>,
    reporter: Arc<Mutex<Reporter>>,
    whitelist: Option<&HashSet<usize>>,
) {
    tracing::info!("{input_filter_params:?}");

    let n_threads = 2;
    let mut smc_bam_reader = rust_htslib::bam::Reader::from_path(sorted_smc_bam).unwrap();
    smc_bam_reader.set_threads(n_threads).unwrap();
    let mut subreads_bam_reader = rust_htslib::bam::Reader::from_path(sorted_sbr_bam).unwrap();
    subreads_bam_reader.set_threads(n_threads).unwrap();

    let mut subreads_records = subreads_bam_reader.records();

    let mut sbr_record_opt = subreads_records.next();

    let mut scoped_timer = ScopedTimer::new();

    let mut timer = scoped_timer.perform_timing();
    let mut cnt = 0;
    let mut channel_filter_by_input_filter = 0;
    let mut sbr_inp_cnt = 0;
    let mut sbr_filter_by_length = 0;

    let is_debug = env::var("RUST_LOG")
        .unwrap_or("INFO".to_string())
        .to_uppercase()
        == "DEBUG";

    for smc_record in smc_bam_reader.records() {
        let smc_record = smc_record.unwrap();

        if sbr_record_opt.is_none() {
            tracing::info!("sbr_record_opt is none. break");
            break;
        }

        if !input_filter_params.valid(&smc_record) {
            continue;
        }

        let smc_record_ext = BamRecordExt::new(&smc_record);
        if let Some(whitelist) = whitelist {
            let ch = smc_record_ext.get_ch().unwrap();
            if !whitelist.contains(&(ch as usize)) {
                channel_filter_by_input_filter += 1;
                continue;
            }
        }

        let mut subreads_and_smc = SubreadsAndSmc::new(&smc_record);
        let smc = gskits::gsbam::bam_record_ext::BamRecordExt::new(&smc_record);

        loop {
            if sbr_record_opt.is_none() {
                tracing::info!("sbr_record_opt is none. break");
                break;
            }

            let sbr = sbr_record_opt.as_ref().unwrap().as_ref().unwrap();
            let sbr_ext = gskits::gsbam::bam_record_ext::BamRecordExt::new(sbr);

            if sbr_ext.get_ch().unwrap() > smc.get_ch().unwrap() {
                break;
            }

            if sbr_ext.get_ch().unwrap() == smc.get_ch().unwrap() {
                if !subreads_and_smc.add_subread(sbr, &oup_params.pass_through_tags) {
                    sbr_filter_by_length += 1;
                } else {
                    sbr_inp_cnt += 1;
                }
            }

            sbr_record_opt = subreads_records.next();
        }

        if subreads_and_smc.subreads.is_empty() {
            tracing::warn!("no subreads for {}", subreads_and_smc.smc.name);
            continue;
        }

        timer.done_with_cnt(1);
        cnt += 1;

        // tracing::de
        if is_debug {
            let all_names = subreads_and_smc
                .subreads
                .iter()
                .map(|sbr| sbr.name.clone())
                .collect::<Vec<String>>();
            let all_names = all_names.join("    \n");
            tracing::debug!("ch={:?}\nsbrnames:\n{}", subreads_and_smc.smc.ch, all_names);
        }

        sender.send(subreads_and_smc).unwrap();
        timer = scoped_timer.perform_timing();
    }

    {
        let mut reporter_ = reporter.lock().unwrap();
        reporter_.channel_reporter.inp_num += cnt;
        reporter_.channel_reporter.filter_by_input_filter += channel_filter_by_input_filter;
        reporter_.sbr_reporter.inp_num += sbr_inp_cnt;
        reporter_.sbr_reporter.filter_by_length += sbr_filter_by_length;
    }

    tracing::info!(
        "subreads_and_smc_generator: speed:{:.4}/s",
        scoped_timer.speed(Some(1000_000_000))
    );
}

pub fn align_sbr_to_smc_worker<T>(
    recv: crossbeam::channel::Receiver<SubreadsAndSmc>,
    sender: Sender<AlignResult>,
    target_idx: &HashMap<T, (usize, usize)>,
    map_params: &params::MapParams,
    align_params: &params::AlignParams,
    oup_params: &params::OupParams,
    reporter: Arc<Mutex<Reporter>>,
    max_subreads: Option<usize>,
) where
    T: std::borrow::Borrow<str> + Eq + Hash,
{
    let max_subreads = max_subreads.unwrap_or(20);
    let mut scoped_timer = ScopedTimer::new();

    let mut max_time = 0;
    let mut max_time_qname = String::new();

    let mut inp_sbrs = 0;
    let mut out_sbrs = 0;
    let mut out_smc = 0;
    let mut fallback_num = 0;
    let mut fallback_rescued_num = 0;

    for subreads_and_smc in recv {
        // tracing::info!("sbr_cnt:{}-{}", subreads_and_smc.smc.name, subreads_and_smc.subreads.len());
        let timer = scoped_timer.perform_timing();
        let start = Instant::now();
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
        if align_infos.len() < max_subreads && !no_hit_indices.is_empty() {
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

        let align_res = AlignResult {
            records: align_infos,
        };
        out_sbrs += align_res.records.len();
        timer.done_with_cnt(1);
        if align_res.records.len() > 0 {
            out_smc += 1;
            sender.send(align_res).unwrap();
        }
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

pub fn align_sbr_to_smc<T>(
    subreads_and_smc: &SubreadsAndSmc,
    sbr_indices: Vec<usize>,
    target_idx: &HashMap<T, (usize, usize)>,
    map_params: &MapParams,
    align_params: &AlignParams,
    oup_params: &OupParams,
    fallback: bool,
) -> (Vec<BamRecord>, Vec<usize>)
where
    T: std::borrow::Borrow<str> + Eq + Hash,
{
    let setup = build_sbr_to_smc_setup(
        subreads_and_smc.smc.seq.len() < 200,
        fallback,
        map_params,
        align_params,
    );

    let ref_aligner: Box<dyn RefAligner> = build_aligner(
        &setup,
        subreads_and_smc.smc.seq.as_bytes(),
        &subreads_and_smc.smc.name,
    );

    let smc_seq = subreads_and_smc.smc.seq.as_bytes();

    let mut bam_records = vec![];
    let mut no_hit_subreads_idx = vec![];

    sbr_indices.into_iter().for_each(|i| {
        let subread = &subreads_and_smc.subreads[i];
        let query_seq = subread.seq.as_bytes();

        let mut hits = ref_aligner.map(&subread.name, query_seq);

        let has_hit = !hits.is_empty();
        for hit in hits.iter_mut() {
            // no supp is needed !!
            if align_params.poly_n_gap_left_align {
                cigar_adjust_poly_gap_left_align(hit, smc_seq, query_seq);
            }

            let identity = cigar::identity_gap_compressed(hit);
            let coverage = cigar::query_coverage(hit, query_seq.len()).max(cigar::target_coverage(hit));
            if coverage < oup_params.oup_coverage_threshold
                || identity < oup_params.oup_identity_threshold
            {
                continue;
            }

            if !subread.name.starts_with("ref") && align_params.fwd_only && hit.is_reverse {
                continue;
            }

            if !subread.name.starts_with("ref") && align_params.rev_only && !hit.is_reverse {
                continue;
            }

            if hit.is_primary && !hit.is_supplementary {
                let bam_record =
                    mapping2record::build_bam_record_from_mapping(hit, subread, target_idx);
                bam_records.push(bam_record);
                break;
            }
        }

        if !has_hit {
            no_hit_subreads_idx.push(i);
        }
    });

    (bam_records, no_hit_subreads_idx)
}

/// https://github.com/PacificBiosciences/actc/blob/main/src/PancakeAligner.cpp#L128
/// https://github.com/lh3/minimap2/blob/master/minimap.h
///
/// 把 subread→SMC 比对所需的引擎参数翻译成与引擎无关的 [`AlignSetup`]。
fn build_sbr_to_smc_setup(
    short_insert: bool,
    fallback: bool,
    map_params: &MapParams,
    align_params: &AlignParams,
) -> AlignSetup {
    let (k, w) = if short_insert {
        if !fallback {
            (4, 1)
        } else {
            (3, 1)
        }
    } else if !fallback {
        (12, 5)
    } else {
        (7, 5)
    };

    let mut setup = AlignSetup {
        preset: crate::align::Preset::MapOnt,
        k: Some(k),
        w: Some(w),
        // 短插入不加 hpc
        hpc: !short_insert,
        best_n: Some(1),
        short_insert,
        // subread→smc 需要 cs 与 md tag
        emit_cs: true,
        emit_md: true,
        sam_hit_only: true,
        secondary_seq: false,
        ..Default::default()
    };

    // this is all map-ont default
    // zdrop = 400; zdrop_inv = 200; bw = 500; end_bonus = 0;

    map_params.modify_setup(&mut setup);
    align_params.modify_setup(&mut setup);

    setup
}

pub fn draw_aligned_seq(
    record: &BamRecord,
    ref_seq: &[u8],
    r_start: Option<usize>,
    r_end: Option<usize>,
) -> (String, String) {
    let mut ref_aligned_seq = String::new();
    let mut query_aligned_seq = String::new();

    let query_seq = record.seq().as_bytes();
    let mut rpos_cursor = None;
    for [qpos, rpos] in record.aligned_pairs_full() {
        if rpos.is_some() {
            rpos_cursor = rpos;
        }

        if let Some(r_start) = r_start {
            if let Some(rpos_cursor) = rpos_cursor {
                if (rpos_cursor as usize) < r_start {
                    continue;
                }
            } else {
                continue;
            }
        }

        let q_char = if let Some(qpos_) = qpos {
            unsafe { (*query_seq.get_unchecked(qpos_ as usize)) as char }
        } else {
            '-'
        };

        let r_char = if let Some(rpos_) = rpos {
            unsafe { (*ref_seq.get_unchecked(rpos_ as usize)) as char }
        } else {
            '-'
        };

        ref_aligned_seq.push(r_char);
        query_aligned_seq.push(q_char);

        if let Some(r_end) = r_end {
            if let Some(rpos_cursor) = rpos_cursor {
                if (rpos_cursor as usize) >= (r_end - 1) {
                    break;
                }
            }
        }
    }

    (ref_aligned_seq, query_aligned_seq)
}

#[cfg(test)]
mod tests {
    use gskits::fastx_reader::fasta_reader::FastaFileReader;

    use crate::mapping2record::build_bam_record_from_mapping;

    use super::*;

    #[test]
    fn test_align() {
        let smc_seq = "CACTCTTTAAAAGGGGGGTTGAGG";
        let subreads_and_smc = SubreadsAndSmc {
            smc: ReadInfo::new_fa_record("hello".to_string(), smc_seq.to_string()),
            smc_len: smc_seq.len(),
            subreads: vec![
                ReadInfo::new_fa_record("".to_string(), "CACTCTTTAAAAGGGGGGGGTTGAGG".to_string()),
                ReadInfo::new_fa_record("".to_string(), "CACTATTTAAAAGGGGGGTTGAGG".to_string()),
                ReadInfo::new_fa_record("".to_string(), "CACTCTTTAAAAGGGGGGTTGAGA".to_string()),
                ReadInfo::new_fa_record("".to_string(), "CACTCTTTAAAAGGGGGGAAGAGG".to_string()),
            ],
        };

        let mut target2idx = HashMap::new();
        target2idx.insert("hello".to_string(), (0, subreads_and_smc.smc.seq.len()));
        let (records, _) = align_sbr_to_smc(
            &subreads_and_smc,
            (0..subreads_and_smc.subreads.len()).collect(),
            &target2idx,
            &params::MapParams::default(),
            &params::AlignParams::default(),
            &params::OupParams::default(),
            false,
        );
        for record in records {
            println!("{:?}", record.cigar());
        }
    }

    #[test]
    fn test_align2() {
        let query_fa = FastaFileReader::new("test_data/query.fa".to_string());
        let query_fa = gskits::fastx_reader::read_fastx(query_fa);
        let query_fa = &query_fa[0];

        let ref_fa = FastaFileReader::new("test_data/ref.fa".to_string());
        let ref_fa = gskits::fastx_reader::read_fastx(ref_fa);
        let ref_fa = &ref_fa[0];

        let setup = AlignSetup {
            k: Some(7),
            w: Some(5),
            emit_cs: true,
            emit_md: true,
            sam_hit_only: true,
            ..Default::default()
        };

        let aligner = build_aligner(&setup, ref_fa.seq.as_bytes(), &ref_fa.name);
        for hit in aligner.map(&query_fa.name, query_fa.seq.as_bytes()) {
            println!("query2ref {:?}", hit)
        }

        let aligner = build_aligner(&setup, query_fa.seq.as_bytes(), &query_fa.name);
        for hit in aligner.map(&ref_fa.name, ref_fa.seq.as_bytes()) {
            println!("ref2query {:?}", hit)
        }
    }

    #[test]
    fn test_build_aligner() {
        let setup = build_sbr_to_smc_setup(
            true,
            false,
            &params::MapParams::default(),
            &params::AlignParams::default(),
        );
        println!("{:#?}", setup);

        let setup = build_sbr_to_smc_setup(
            false,
            false,
            &params::MapParams::default(),
            &params::AlignParams::default(),
        );
        println!("{:#?}", setup);
        assert_eq!(setup.k, Some(12));
        assert_eq!(setup.w, Some(5));
        assert!(setup.hpc);
    }

    #[test]
    fn test_short_insert() {
        let setup = build_sbr_to_smc_setup(
            true,
            false,
            &params::MapParams::default(),
            &params::AlignParams::default(),
        );
        println!("{:#?}", setup);
        assert!(!setup.hpc);
        assert!(setup.short_insert);
        assert_eq!(setup.k, Some(4));
        assert_eq!(setup.w, Some(1));

        let ref_seq = "CCCCTTTTAAAAGGGGGAGAGA";
        let q_seq = "CCCCTTTTAAAAGGGGGAGAGA";
        let aligner = build_aligner(&setup, ref_seq.as_bytes(), "ref");

        for hit in aligner.map("query", q_seq.as_bytes()) {
            println!("{:?}", hit);
        }
    }

    #[test]
    fn test_align_two_seq() {
        let icing_reads = "CCGCGCCAATTATTGTCGATGGCGACGGCGCATTGCCCGTCAGGATTACCTTGTGGATCAATCATTTCCACCATCAGCAATTTGTCATGCGCCAGTCCGTGATGGCGTACGGTACAAACAATTTGTCCAGTGACGACTGCCAGTTTCATACCCGCCTCCGTGGCGTATTTCAGGTAAAAGCTCCCCCTACCCTCCGCAGAAGGTAAAATGAAAAAGGAGAGAGCGTGACGCCCGAATCGACGTCACACAGGGTGATTACAGGTTGCTGCTATCGCCTTTCAGGCCGATCGGGAAGACTTCTTCCAGGTCGCCGTGCGGACGCGGGATAACGTGTACAGATACCAGCTCGCCGATACGCTGTGCCGCTGCGGCACCTGCATCGGTCGCGGCTTTACAGGCTGCAACATCACCACGAACCATGGCAGTACACAGGCCGCCGCCAATCTGCTTAACACCAACCAGCTTGACGCGCGCAGCTTTCACCATCGCGTCAGAAGCCTCAATCAGTGCAACCAGGCCCCGGGTTTCGATCATTCCTAATGCTTCCATTGTGCTTTCCTCTTTATCAGGGTCCAGAACGGGACCGTTCATTCAACCAGTGTTTGTAACTGCTTTCGCGGTTCACTTCTGTCTGACGCGGCACAGCTGCCACCAGCGCCAGCTCGATAATTTCCTGCACGCTACAACCACGAGAGAGATCGTGCATCGGCGCGGCAAGTCCTTGTATCAGTGGCCCGACGGCACGATATCCGCCGAGTCGTTGTGCGATTTTGTAACCAATATTTCCGGCTTCCAGCGACGGAAAAACCATCACATTGGCCTTGCCCTGTAGCGGGCTGGCAGGCGCTTTTTGCGCCGCCACTTCCGGCACGAAGGCGGCGTCAAACTGTAACTCGCCATCCACCACCAGCTTTGGTGCGCGCTCACGGACGATTTCTGTCGCCTGCTGGACGTTAGCAACACAGGGGTGACGGGCGCTACCATTGCTGGAAAACGACAGCATCGCCACGCGCGGCTCTTCTCCGGTGATGGCGCGCCAGGTTTCGGCACTGGCAA";
        let sbr = "CCGCGCCATATTTTTTGATAATAACTTGGTATTGGACTTTTTGATTGTCGGATGGCGACGGCGCATTGCCGTCAGGATTACCTTGTGATCAATCATTCCACCATAACAAATTTGTCATGCGCCAGTCCGTGATGGCGTAACGTTGTATTGTCCCGTGACGACTGCCAAGTTTCAACCCGCTCCGTGGCGTATTTCAGGTAAAAGCATTCCCCCTA
CCCTCCGCAGAAGGTAATAAAATGAAAAAGGAGAGAGACGTGACAGCCCGGAATCGCCCTGTCACATCTGGGTGATTAACAGGTTGCTGCTATTCGCCTTTCAAAGGGCCGATCGAGGAACACCTTTCTCCAGTCCGTGGACGCGTGGATACGTGTACAGATACCAGCTCGCCGATACGCTGTGCCGTGACGCACCTGCATCGGTCGCGGCTTTACAGCAACATCACCACGAACCATGCGTACACAGCCGCCGCATCTGCTAAAACCCATAACCACCAGCTTGACGCGCGCAGCTTTCCATCGCGTCCAGA
AAGCCTCAATCAGTGCAACCAGGCCCCGGGTTTCGATCATTCCTAATGCTTCCATTGTGTTTCCTCTTATATCAGGTCCAGAACGGACCGTTCATTCACAACCGTGTTTGTAAACGGCTTTCGCGGTCACTTCCTGTCTGACGCGGCACGCTGGCCACCAGCGCCAGCTCGATTATTTCTTGCACGCTACAACCACGAGAGAGAGAGATCGTGCATCGGCCGGCAAGTCCCTTGTATCAGTGGCCGACGGCACGATATCCGCCGGTCGTTGTGCGATTTTTAACCATATTTCCGGGCTTCCAGCGACAGGA
AAAACCATCACATTGGCCTTGCCTGTAGCGGGCTGGCAGGCAGCTTTTTGCGCCCCACTTCCGCACGAGGCAGGCGTCACAACTGTAACTCGCCTCCACACCAGCTTTGGTGCAGCGCTCACGGACTGATTTCTGGCCGCTGCTGGACGTTAGCAACAACCGGGGTGACGGGCGCTACCATTGCTGGAAAACGACACGCCACGCGTCGGCTCTTCCCGGTGATGGCGCGCCAGGTTTCGGCACTCAGCAAA";
        let setup = build_sbr_to_smc_setup(
            false,
            false,
            &params::MapParams::default(),
            &params::AlignParams::default(),
        );
        let aligner = build_aligner(&setup, icing_reads.as_bytes(), "icing");
        let mut target_idx = HashMap::new();
        target_idx.insert("icing".to_string(), (0, icing_reads.len()));
        let query_len = sbr.len();
        for hit in aligner.map("sbr", sbr.as_bytes()) {
            if hit.is_primary {
                println!(
                    "identity_gap_compressed:{}",
                    cigar::identity_gap_compressed(&hit)
                );
                println!("identity:{}", cigar::identity(&hit));
                println!("query_coverage:{}", cigar::query_coverage(&hit, query_len));
                println!("target_coverage:{}", cigar::target_coverage(&hit));
                let query_record = ReadInfo::new_fa_record("name".to_string(), sbr.to_string());
                let record = build_bam_record_from_mapping(&hit, &query_record, &target_idx);

                let (ref_aligned, query_aligned) =
                    draw_aligned_seq(&record, icing_reads.as_bytes(), Some(590), Some(620));
                println!("\n{}\n{}", ref_aligned, query_aligned)
            }
        }
    }
}
// GGTGCTGGCCTATCTAGTTGTTGGGTTTGGGATGGGGCTGCTCTGCCTATGCAGTAGTCATTCGTTTATTTCTTTACCTCTCAAGAAATTTGCTTTGACTTAAAAAAAACAAAAACAAAAAAGAAAGAAAGGCATATTCTGAAACGGTAGTGACTGCATTTGAAGTGGGACCTTTTATGATATGATGGAGCCAGGCAGAGATGTTGGAGAAGTTCCTGAAAGAAGAAGGGCATGTGCCAAATTCTGAGGCTGAGGAGAAAAAGAAAAAAGAAAAAAAGAATAAAGAACTTTACATTTCACTGTATGTAAGGACTTACAAGGCTAGAGTAAAGCATGTTGAAGTAAAAATAGGAGAAATCAAGTTAGAGAGAAGGCGCAGGCTTATTATTTGGGTCTTATAGATGAGAGAGTAGAGTAGGTAATTTATACTGAAACATAGGGGACGAGGGGGAGCTTTGCAAACCTGAGATAAACAGGTCTGTAGTCCACTCAAACCCAGGCTTATTTTTTTCTCTTCTTCATCTGCCTTACTCATGTAAGGTGCTATAACAAAATACAAAGACTGGGTGGCTTCAACACAGAGTTTTTTCTCACATTTCTGAAGGCTGACTCGTGGTCAGGGTTCTCTTCCTGGATTGCAGATGGCGGACTTCTCACTGTGTCCTCACCGGGAGGAAAGAGAAAGGCTTGTCTCTCTGCTTCTGATAAGGAAAATAATTTTATGATGGGGGTCTGCTCTTATGAGCTCAGCTAAACTAATTACTTTTCAAAAGCCTCCCCCACAGATAAGTCATTGGGAGTGGGAATTCAGCATGTGAATTTTGAGGAGACACAAACATTAGGTCCGTAACACCATCATATTAAACTCTTACTTATAGAGATGCACTAACATCCAACTATACAAAAATTCCATTAAGTAAAACCTTCCAATTAATGTTGTACAAATTCTTCCCACTACATTTTGATCTCACAGTGCTGGTCTGTTTTCAGACATAACTGCTCATGCTGACTTGTGAGCCTCTGCTGATTCATTTCTTACACCTAAGTCCCCACCTCCCCCAGAGTCTTACAATGAGGTCTCTCCCTTGAAATGCTGTGACCAACTGCACACTTAGGCATGGGATTAATTTCACTTTCTTAGGCATCCACAGGGCGGTGAAAAGCTAAGTGCCATTAGAATAGTATTCCTCCTAATAAAATAGAATTGCTACACCTGACTAAACCTTGAATTCCTAAAGCTTGGAACACTTTCCCTTCATTAAGAACCATCCTTGCTACTCAGCTGCAATCAATCAGCCCCCAGGTCTTCACTGAACCTTTTCCATCTCTTCCAAAACATCTGTTTCTGAGAAGTCCTGTCCTATAGAGGTCTTTCTTCCCACCGGATTTCTCCTACACATTTACCCCCACTTGCAGAACTCCCGTGTACAAGTGTCTTTACTGCTTTTATTTGCTCAACAAAAGCACATCTCATATAAAAAAAATGAGGAGCATGCACACACCACAAACACAAACAGGCATGCAGAAATACACATACACACTTCCCTCAATATAAACCCTTTGTGGCTCATATATTTAAAAAGGGAAATGTAAAACAAAAAAGAGCTGAAGAAACATGTGTGATCTCCTCAGCAGAATAGATTTATTATTTGTATTGCTTGCAGATAAAGCCTATCCTTGAAAGCTCTGAATCATGGGCAAGAGGCTCAGTGGTATCTGGAGGACAGGGCACTGCCACTGCAGTCACCATCTTCTGCCAGGAAGCCTGCACCTCAGGGGTGAATTCTTTGCCAAATGGATTGCCAAAACGGTCACAGCACATTTCCCAGGAGCTGTTGAGATGAAAGGAGACAATAAAGATGAACCCATAGTGAGCTGAGAGCTCCAGCCTGGCCTCCAGATAACTACACACCCAGCTTCTACCCAAATCAAGCCTATGTTAACTCCCTCAAAGCCTGAGATTTTGCCTTCCCATTAAATGCAGGGTAGTTGTTCCCCTCAACACTAGTCACTGCCATAATTTAAATCTTGCTATCTTCTTGCCACCATGAACCCTGTATGTTGTAGGCTGAAGACGTAAAAGAACACACGCTGACACACACACACACACACGCGCGCGCGCACACACACACACACACAGAGCTGACTTTCAAAATCTACTCCAGCCCAAATGTTTCATTGTTCTCACCCTGGACATACTTTGCCCCCATCTGGAATTAAAGGATAAAGTTTGTAATGAAGCATTAGCAGCATTCTATATGGTCCAGCTGATATAGGAATAGCCTTAGCAATGTATGTTTGGCCACCAAGTTCCCCATTTGACTGAGCCAATATATGCCTTCTGCCTGCACTTTTTAACGACATACTGTCCTGCTCAGATAGATGTTTTAAAAAACAAAAATGAGGGAAAGATGAAAGTTCTTTCTACTGGAATCTAATAAAAGTCATTTTCCTCATTCCACCTCTCTTTTCTCAAAGTCAAAATGTCCATCTAGATTTTTAGAGGCGCTCCTTAGGCCCTAAAACATTACCACTGGTCTCAGCCAGCTAGTCCTCTGCAGTTTCTTCACCCCAACCCCAGTATCTTCAAACAGCTCACACCCTGGTGCTCGATCATACTCAGTTGTCTAAGTTGCCTCGAGACTAAAGGCAACAGGCTGAACATCTCCTGACTCACCTTGAAGTTCTCAGGATCCACATGCAGCTTGTCACAGTGCAGTTCACTCAGCTGGGCAAAGGTGCCCTTGAGATCATCCAGGTGCTTTATGGCATCTCCCAAGGAAGTCAGCACCTTCTTGCCATGTGCCTTGACTTTGGGGTTGCCCATGATGGCAGAGGAGAGGACAGTTGCCAAAGCTGTCAAAAACACCCTGGGTCATGGGTAGACAACCAGGACCTGTGAGATTGACAAGAACAGTTTGACAGTCAGAAGGTGCCACAAATCCTGAGAAGCAACCTGGACTTTTGCCAGGCACAGGGTCCTTCCTTCCCTCCCTTGTCCTGGTCACAAGCCTACCTCCAGGGTTTCTCCCCAGCATCTTCCACATTCCCTTGCCCCACAGGCTTGTGATAGTAGCCTTGTCTCCTCTGTGAAATGACCCATGGCGTCTGACTAGGAGCTTATTGATAACCCAGACGTTCCAGAAGCGAGTGTGTGGAACTGCTGAAGGTGCTTCCTTTTATTCTTCATCCCTAGCCAGCCGCCGGCCCCTGGCCTCACTGGATATACTCTAAGACTATTGGTCAAGTTTGCCTTGTCAAGGCTATTGGTCAGGGCAAGGCTGGCCAACCCATGGGTGGAGTTTAGCCAGGGACCGTTTCAAACAATATTTGCATTGAGATAGTGTGGGGAAGGGGCCCCCAAGAGGATACTGCTGCTTAATTTTTTTATAGCCTTTGCCTTGTTCCGATTCAGTCATTCCAGTTTTTCTCTAATTTATTCTTCCTTTAGCTAGTTTCCTTCTCCCATCATAGAGGAACCAGGACTTCTTTTGCAAGCCGTTTTTAACCTCTTGCCTCTAGCTCCAGTGAGGCCTGTAGTTTAAAGCTAAAGCATGTACCAATTTTTGAAAAGTCAGGGATTGTGAAATGTGTTTTAGGCATAGGTCCAGGATTTTTGACGGGACAATCTTAGTCTCTTTCAGTTAGCAGTGGTTTCAAAGGAAAAGTGCTATACTCTTTTTGAAATATACTCTTTGTGACTTTTGCCATATCTCTTAATCTCTCATAGTGCAGTGAAACAATTTCTATAAAGCCACAGTTTCAGCGCAGTAATAGATTAGTGTTACATAATATAGACCTAATGCTTACCTCAAATCTACTTATCCGTACCTATTTGAACATAAAATCATGACTGTTTCATCTTAGAAAAATAATTTGATTCCATATTCAGGTATGTAAGTATACACCAGATGATGTGTATTTACCACTGATAAGTGTGTGTGCTGGCTGATGACCCAGGGTTTTGGCGTAGCTCTTCTATGCTCAGTAAAATGAGGGAGAATGTTCTTTGGCAGGTACTGTGGATTAGAATTAATTATCTTGTATAAATGCTAGGTTCACTTCTCAGGGAATCTTACTCTAAGAATAAGATGTGCGTGTACATGGAAAACAACTCTAAAGAGGCAGGGTTGTTTTTATTGACTAATAGTCCACACACTATTATAACTCGAAATAAGTGTACTTTAGACAGCTTTATTTCTAACACAGCGCTGTTCTGACATATTGGACCATTAACAGGGTAAGAAGTATTTATGGTGGTTTTTTGGTTCTGTTTTGCTTTGGTTAGTTTGTTTTTGTTTTTCTCCGAAAGTGATCCAGATCTCTAACCTTGCTAGATTATAATGCCAGAAGCTCTGGAATTCTGGCTTATCGGAGGCAAGCTGTATCTTCAAATTAGTTTATCCCCTAAGCTATCAGGTTGATTGAAATTATATAATATTGGTGAAATTCTTTCATCCTTCATGATCCTGTGTAAAGCTTATACCTGCTCATTGAACAGATGGCTAAAAGGCCAAAAGACAACACCATTAACCCATGTAACAAACCTGCACATCCTGCACATGTACCTTCAAAATTAGAATGAAATGAAATAAAATTAAAAGGAAAAAATGCCATGGCAACATCAGAAGTATCTTTTAGTCTAAAAATGGGGCAGCATAAATAATCTACCCCTTGTTTAGCATACTATTGAAAAATAACAATAAAAATCGGCAACCAGTAGCCCTTGCGTCTACTCTGCCTATGAAGTAGCCATTCATTTATTCCTTCAATTTTTTAAAACTTGTTTACTAAAAAAAAAGGACAAAGAAAGAAAGTGCGAATTGTGAAATGGTAGTGAGTGATGGCATTTGAAAGTGGGTCCTTTATGATTTGATGGAGCCAGGCAAAGACGTTTGGAGAAAGAAGTTCCTGAAAGTAGGAAGGGCATGTGGAAAACTCTGAGGCTGGCCTATCTAGGTTGTTGGGTTAACCAAACCCAACAACCTAGATAGGCCAGCCTCAGAGTTTCCACATGCCCTTCCTTCTTTCAGAACTTCTTTCTCCAAACGTCTTCTGCTGCTCCATCAAATCATAAGGACCCACTTCAAATGCCATCACTCACTACCATTTCACAATTCGCACTTTCTTTCTTTGTCCTTTTTTTTTAGTAAAACAAGTTTATAAAAAATTGAAGGAATAATGAATGGCTACTTCATAGGCAGAGTAGACGCAAGGCTATTGTTGCCATTTTTATTGTTATTTTTTCAATAGTAGCTAAACAAGGGGTAGATTATTTATGCTGCCCATTTTAGACCATAAAAGTAACTTCCTGATGTTGCCATGGCATTTTTTCCTTTTAATTTTATTTCATTTTCATTCTATTTCGAAGGTACATGTGCAGGATGTGCAGGTTTGTTACATGGTTAAATGTGTGTCTTTCGGCCTTTTAGCCATCTGTATCAATGAGCAGATATAAGCTTTACACAGGATCATGAGGTGAAAGAATTTCACCAATATTATAATAATTTCAATCAACCTGATAGCTTAGGGATAAACTAATTTGAAGATACAGCTTGCCTCCGATAAGCCAGAATTCCAGAGCTTCTGGCATTATAATCAGCAAGGTTAGAGATCATGGATCACTTTCGGAGAAAACAAAAACAAACTAACCAAAAGCAAAACAGAACCAAAAAACCACCATAAATACTTCCTACCTGTTATAATGGTCCAATATGTCAGAAACAGCGCTGTGTTAGAAATAAAGCTGTCTAAGTACACTAATATTCGAGTTATAATAGTGTGTGGACTATTAGTCAAAAAAAACAACCCTTGCCTCTTAGAGTTGTTTTCCATGTACACGCACATCTTATGTCTTAGAGTAAATTCCCTGAGAAGTGAACCTAGCATTTATACAAGATAATTAATTCTAATCCACAGTACCTGCCAAGAACATTCTACCATCATCTTTACTGAGCATAGAAGAGCTACGCCAAACCCTGGTCATCAGCAGCACACACACTTATCCAGTGGTAAATACACATCATCTGGTGTATACATACATACCTGAATATGGAATCAAATATTTTTCTAAATGAAACAGTCATGATTTATTTCAAATAGGTACGGATAAGTAGATATTGAGGTAAGCATTAGGTCTTATATTATGTAACACTAATCTATTACTGCGCTGAAACTGTGGCTTTATAGAATTGTTTTCACTGCACTATTGAGAGATTAAGAGATAATGGCAAAAGTCACAAGAGTATATTCAAAAAGAAGTATAGCACTTTTTCCTTAGAAACCACTGCTAACTGAAAGAGACTAAGATTTGTCCCGTCAAATCCTGGACCTATGCCTAAAACACATTTCACAATCCTGACTTTCAAAAATTGGTACATGCTTTAGCTTTAAACTACAGGCCCCACTGGAGCTAGAGGCAAGAGGTAAAAAACGGCTGACAAAAGAAGTCCTGGTATCCTCTATGATGGGAGAAGGAAACTAGCTAAAGGGAAGAATAAATTAGAGAAAACGGAATGACTGAATCGGAACAAGGCAAGGCTATAAAAAATTAAGCAGCAGTATCCTCCTGGGGGCCCCTTCCCCACACTATCTCAATGCAAATATCTGTCTGAAACGGTCCCTGGCTAAACTCCACCCATGGGTTGGCCAGCCTTGCCCTGACCAATAGCCTTGACAAGGCAAACTTGACCAATAGTCTTAGAGTATCCAGTGAGGCCAGGGCCGCGCTGCTAGGATGAAGAATAAAAGGAAGCACCTTCAGCAGTTCCACACTCGCTTCTGAACGTCTGAGGTTATCATAAGCTCCTAGTCCAGACGCCATGGTCATTCACAGAGGAGGACAAGGCTACTATCAAAGCCTGTGGGCAGGTGAATGTGGAAGATGCTGGAGGAGAAACCCTGGAAGGTAGGCTCTGTGACAGGACAAGGGAGGGAAGGAAGGACCCTGTGCCTGGCAAAAGTCCAGGTTGCTTCTCAGGTTTGTGGCACCTTCTGACTGTCAAACTGTTCTTGTCAATCTCACAGCTCCTGGTTGTCTACCCATGACCCAGGGGTTCTTTGACAGCTTTGCAACTGTCCTCTGCCTCTGCCATCATGGGCAACCCCAAGTCAAGGCACATGGCAAGAAGGTGCTGACTTCCTTGGAGATGCCTAAAGCACCTGGATGATCTCAGGCACCTTTGCCCAGCTGAGTGAACTGCACTGTGACAAGCTGCATGTGGAATCCTGAGAACTTCAAGTGAGTCCAGGAGAGTTTCAGCCCTGTTGCCTTTAGTCTCGAGGCTACTTAGACAACTGAGTATTGATCTGAGCACAGCAGGGTGTGAGCTGTTTGAAGATACTGGGTTGGGGGTGAAGAAACTGCAGAGGACTAGCTGGGCTGAGACCCAGTGGTAAATGTTTTAGGGCCTAAGGAGCGCCTCTAAAATCTAGATGGACAATTTTGACTTTGAGAAAAGAGAGATGGAAATGAGAAAAATGACTTTTATTAGATTCCAGTAGAAAGACTTTCATCTTTCCCTCATTTTGTTGTTTAAAACATCTATCTGGAGGCAGGACAAGTATGGTCGTTAAAAAGATGCAGGCAGAAGGCATATATTGGCTCAGTCAAGTGGGAACTTTGGTGGCCAAACATACATTGCTAGGCTATTCCTATATCAGCTGGACACATATAGAATGCTGCTAATGCTTCATTACAAACTATATCCTTTAATTCCAGATGGGGCAAAGTATGTCCAGGGTGAGGAACAATTGAAACATTTGGCTGGAGTAGATTTTGAAAGTCAGCTCTGTGTGTGTGTGTGTGTGCGCGCGCGCGTGTGTGTGTGTGTGTGTCAGCGTGTGTTTCTTTTAACGTCTTCACCTACACATACAGGGTTCATGGTGGCAAGAAGATAGCAAGATTTAAATTATGGCCAGTGACTAGTGCTTGAAGGGAACAACTACCTGCATTTAATGGGAAGGCAAAATCTCAGGCTTTGAGGGAAGTTAACATAGGCTTGATTCTGGGTAGAAGCTGGGTGTGTAGTTATCTGAGGCAGGCTGGAGCTCTCAGCTCACTATGGGTTCATCTTTATTGTCCCTTTCATCTCAACAGCTCCTGGGAAATGTGCTGGTGCCGTTTTGGCAATCCATTTCGGCAAAGAATTCACCCCTGAGTGCAGGCTTCCTGGCAGAAGATGGTGACTGCAGTGGCCAGTGCCCTGTCCTCAGATACCACTGAGCCTCTTGCCATGATTCAGACTTCAAGGATAGGCTTTATCTGCAAGCAATACAAATAATAAATCTATCTGCTGAGAGATCACACATGATTTTCTTCAGCTCTTTTTTTTTAAATCTTTTTAAATATATGAGCCACAAAGGTTTATATTGAGGGAAGTGTGTAATGTGTATTTTCTGCATGCCTGTTTGTGTTTGTGGTGTGTGCATGCTCCTCATTTATTTTATATGAGATGTGCATTTTGTTGAGCAATAAAAGCAGTAAAGACACTTGTACACGGGAGTTCTGCAAGTGGGGGTAAATGGTGTAGGGAAATCCGTGGAAGAAAGACCTCTATAGGACAGGACTTCTCAGAAACAGAATGTTTTGGAAGAGATGGGGAAAAGGTTCAGTGAAGACGGGGGCTGGTTGATTGCAGCTGAGTAGCAAGGATGGTTCTTAATGAAGGGAAAGTGTTCCAAGCTTTAGGAATTCAAGGTTTAGTCAGGTATGCAATTCCTATTTTATTAGGAGGAATACTATTTCAATGGCACTTAGCTTTTCACCGCCTGTGATGCCTAAGAAAGTGAAATTAATCCCATGCCCTCAAGTGTGCAGATTGGTCACAGCATTTCAAGGGAGAACCTCATTGTAAGACTCTGGGGAGGTGGGGACTTAGGGTAAGAAATGAATCAGCAGAGGCTCACAAAAGTCAGCATGAGCATGTTATGTCTGAGAAACAGACCAGCACTGTGAGATCAAAATGTAGTGGGAAGAATTTGTACAACAATTAATTGGAAGGTTTACTTAATGGAATTTTTGTATAGTTGGATGTTAGTGCATCTCTATAAGTAAGAGTTTAATATGATGGTGTACGGACCTAATGTTTGTGTCTCCTCAAAATTCACATGCTGAATCCAACTCCCAACTGACCTATCTGTGGGGGAGGCTTTTGAAAAGTAATTAGGTTTAGCTGAGCTCATAAGAGCAGACCCCATCATAAAAATTATTTCCTTATCAAGAAAGAAGAGACAAGCTCCATTTCTCTTTCCTCCCGTGAGGACACAAGTGAGAAGTCCGCCATCTGCAATCCAGGAAGAGGAACCCTGACCACGAGTCAGCCTTCAGAATTGTGAGAAAAACTCTGTTGTTGAAGCCCACCCAGTCTTTTGTATTTTGTTATAGCACCTTACACTGAGTAAGGCAGATGAAGAAGGAGAAAAAAATTAAGCCTGGGTTTTTGAGTGGACTAACAGACCATGTTTTTATCTTCAAGGGTTTGCATAAAGCCTCCCTCCGTCCCCTATGTTTCCAGTATAAAATACCTACTCCTAAAAAACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCTTTTTTTCCCCCCCCCCCCCCCCCCCCCCCCCTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCTTTTTAAAACCCCCCCCCCTTCCCTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCATCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCTTTTTTTATAAGACCCAAAAATAAGCCTGCGCCCCCCCCCCCCCTTTTTTTTTTTTCTCTCTATACTTTTTGATTTCTCCTATTTTTACTTCAACTATGCTTTTACTCTAGCCTTGTAATGTCCTTACATACAGGAAATGTAAAGTTCTTTATTCTTTTTTTCTTTCTTTCTTTTTTTCTCCTCAGCCTCAGAAATTTGGCACATGCCCTTCCCTCTTTCAGGAACTTCTCCACATCTCTGCCTGGCTTCCATCATTCATAAAGGTCCCACTTCAAATGCAGTCACTACCGTTTCAGAATATGCACTTTCTTCTTTTTTGTTTGTTTTTTTTAATAGTCAAAGCAATATATTTTCTTGAGAGAGTAAAATATAAACTGATGACACTGCATACCAAGCATCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCAAAAAAAAAAATTTTTTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCCCCCTTTTAAAGGGGGGGGAAAAAAAAAAAAAAAAAAAAAAATTTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCCCCCCCCAAGGCCCCCCCAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCCCAAATTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCTCCCCCCCCCCCCCCCCCCCCTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCTTTTTTTTCCCCCCCCCCCCCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTATTTTTTTTTTTTTTTTTTTAAAAAATTTTATAATATTAAAAAAA
