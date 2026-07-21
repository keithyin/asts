use std::{
    collections::HashMap,
    sync::{Arc, Mutex},
    thread,
    time::Instant,
};

use crate::{
    align_sbr_to_smc, mapping2record, read_info::ReadInfo, reporter::Reporter, utils::Range,
};
use crate::{params, SubreadsAndSmc};
use crossbeam::channel::Sender;
use gskits::dna::reverse_complement;
use minimap2::{Aligner, Built, Mapping};
use mm2::gskits::{
    gsbam::{
        bam_record_ext::BamRecordExt, plp_counts_from_records::compute_max_ins_of_each_ref_position,
    },
    utils::ScopedTimer,
};
use mm2::mapping_ext::MappingExt;
use rust_htslib::bam::{ext::BamRecordExtensions, Record};
use tracing;

use serde::Serialize;

#[derive(Debug, Default)]
struct Cs2RefAlnRes {
    identity: f32,
    query_align_len: u32,
    mm: u32,
    ins: u32,
    homo_ins: u32,
    del: u32,
    homo_del: u32,
    pub is_reverse: bool,
    pub ref_sub_seq: Option<String>,
}

impl Cs2RefAlnRes {
    fn new(hit: &Mapping, ref_seq: &str) -> Self {
        let target_start = if hit.target_start > 10 {
            hit.target_start - 10
        } else {
            0
        } as usize;

        let target_end = (hit.target_end + 10).min(ref_seq.len() as i32) as usize;
        let ref_sub_seq = Some(ref_seq[target_start..target_end].to_string());

        let identity = MappingExt(hit).identity();

        Self {
            identity: identity,
            query_align_len: (hit.query_end - hit.query_start) as u32,
            mm: 0,
            ins: 0,
            homo_ins: 0,
            del: 0,
            homo_del: 0,
            is_reverse: hit.strand.eq(&minimap2::Strand::Reverse),
            ref_sub_seq: ref_sub_seq,
        }
    }
}

fn align_cs_to_ref(
    cs: &str,
    ref_aligner: &Aligner<Built>,
    ref_seq: &str,
    query_name: &str,
    ref_range: Option<&mut Range<usize>>,
) -> Option<Cs2RefAlnRes> {
    let hits = ref_aligner
        .map(
            cs.as_bytes(),
            false,
            false,
            None,
            Some(&[67108864, 68719476736]), // 67108864 eqx, 68719476736 secondary
            None,
        )
        .unwrap();

    let ref_bytes = ref_seq.as_bytes();

    if let Some(ref_range) = &ref_range {
        let mut has_error_in_ref_range = false;

        for hit in &hits {
            if !hit.is_primary {
                continue;
            }
            let read_info = ReadInfo::new_fa_record("cs".to_string(), cs.to_string());
            let mut target_idx = HashMap::new();
            target_idx.insert(hit.target_name.as_ref().unwrap().as_str(), (0, 0));
            let align_record =
                mapping2record::build_bam_record_from_mapping(hit, &read_info, &target_idx);

            let mut r_cursor = None;
            let query = BamRecordExt::new(&align_record).get_seq();
            let query_bytes = query.as_bytes();

            for [qpos, rpos] in align_record.aligned_pairs_full() {
                if rpos.is_some() {
                    r_cursor = rpos;
                }
                if r_cursor.is_none() {
                    continue;
                }

                let r_cursor_value = r_cursor.map(|v| v as usize).unwrap();
                if rpos.is_none() || qpos.is_none() {
                    if ref_range.within_range(r_cursor_value) {
                        has_error_in_ref_range = true;
                        break;
                    }
                }

                if let (Some(qpos_), Some(rpos_)) = (qpos, rpos) {
                    let qpos_ = qpos_ as usize;
                    let rpos_ = rpos_ as usize;
                    if ref_bytes[rpos_] != query_bytes[qpos_] {
                        if ref_range.within_range(r_cursor_value) {
                            has_error_in_ref_range = true;
                            break;
                        }
                    }
                }
            }
        }

        if !has_error_in_ref_range {
            tracing::info!(
                "no error in interested ref range: query_name:{} ",
                query_name,
            );

            return None;
        }
    }

    return if hits.len() == 1 {
        let hit = &hits[0];

        ref_range.map(|range| range.shift(-hit.query_start as i64));

        Some(Cs2RefAlnRes::new(hit, ref_seq))
    } else {
        let mut hits = hits;
        hits.sort_by_key(|hit| hit.query_start);

        let infos = hits
            .iter()
            .map(|hit| {
                format!(
                "qstart:{}, qend:{}, rstart:{}, rend:{}, primary:{}, supp:{}, rev:{:?}, score:{:?}, identity:{:.2}",
                hit.query_start,
                hit.query_end,
                hit.target_start,
                hit.target_end,
                hit.is_primary,
                hit.is_supplementary,
                hit.strand,
                hit.alignment.as_ref().unwrap().alignment_score,
                MappingExt(hit).identity()
            )
            })
            .collect::<Vec<String>>()
            .join("\n");

        tracing::info!(
            "align_cs_to_ref: query_name:{} \n{}\n------------------------------",
            query_name,
            infos
        );

        None
    };
}

#[derive(Debug, Serialize, Default)]
pub struct MsaResult {
    pub identity: f32,
    pub query_aln_len: u32,
    pub mm: u32,
    pub ins: u32,
    pub del: u32,
    pub msa_seqs: Vec<String>,
    pub names: Vec<String>,
    pub positions: Vec<i32>,
    pub low_q: Option<u8>,
    pub strands: Vec<String>,
    pub smc_ranges: Option<Range<usize>>,
}

impl MsaResult {
    pub fn extract_error_region(self) -> Self {
        assert!(self.names[2].starts_with("ref"), "{}", self.names[2]);

        let smc_aligned = self.msa_seqs[1].as_bytes();
        let ref_aligned = self.msa_seqs[2].as_bytes();

        let tot_len = smc_aligned.len();
        let mut regions = vec![];

        if let Some(smc_ranges) = &self.smc_ranges {
            (0..tot_len).for_each(|pos| {
                let major = self.positions[pos];

                if smc_ranges.within_range(major as usize) {
                    if regions.is_empty() {
                        regions.push((pos, pos));
                    } else {
                        let last = regions.last_mut().unwrap();
                        if (last.1 + 1) == pos {
                            last.1 += 1;
                        } else {
                            regions.push((pos, pos));
                        }
                    }
                }
            });
        } else {
            (0..tot_len).for_each(|pos| {
                if smc_aligned[pos] != ref_aligned[pos] {
                    let start = if pos > 20 { pos - 20 } else { 0 };
                    let end = (pos + 20).min(tot_len);
                    regions.push((start, end));
                }
            });
        }

        let regions = merge_intervals(regions);

        let mut region_positions = regions
            .iter()
            .flat_map(|&(s, e)| {
                self.positions[s..e]
                    .iter()
                    .copied()
                    .chain(vec![-1].into_iter())
            })
            .collect::<Vec<_>>();
        if let Some(last) = region_positions.last().copied() {
            if last == -1 {
                region_positions.pop();
            }
        }

        let new_msa_seqs = self
            .msa_seqs
            .iter()
            .map(|ori_str| {
                regions
                    .iter()
                    .map(|&(s, e)| ori_str[s..e].to_string())
                    .collect::<Vec<_>>()
                    .join("#")
            })
            .collect::<Vec<_>>();

        if !region_positions.is_empty() {
            assert!(region_positions.len() == new_msa_seqs[0].len());
        }

        Self {
            msa_seqs: new_msa_seqs,
            positions: region_positions,
            ..self
        }
    }

    pub fn extract_lowq_region(self) -> Self {
        assert!(self.names[2].starts_with("ref"), "{}", self.names[2]);

        let qual = self.msa_seqs[0]
            .as_bytes()
            .iter()
            .map(|&v| if v == '-' as u8 { 10 } else { v - '0' as u8 })
            .collect::<Vec<_>>();

        let tot_len = qual.len();
        let mut regions = vec![];
        (0..tot_len).for_each(|pos| {
            if qual[pos] < self.low_q.unwrap_or(10) {
                let start = if pos > 10 { pos - 10 } else { 0 };
                let end = (pos + 10).min(tot_len);
                regions.push((start, end));
            }
        });
        let regions = merge_intervals(regions);

        let mut region_positions = regions
            .iter()
            .flat_map(|&(s, e)| {
                self.positions[s..e]
                    .iter()
                    .copied()
                    .chain(vec![-1].into_iter())
            })
            .collect::<Vec<_>>();
        if let Some(last) = region_positions.last().copied() {
            if last == -1 {
                region_positions.pop();
            }
        }

        let new_msa_seqs = self
            .msa_seqs
            .iter()
            .map(|ori_str| {
                regions
                    .iter()
                    .map(|&(s, e)| ori_str[s..e].to_string())
                    .collect::<Vec<_>>()
                    .join("#")
            })
            .collect::<Vec<_>>();

        if !region_positions.is_empty() {
            assert!(region_positions.len() == new_msa_seqs[0].len());
        }

        Self {
            msa_seqs: new_msa_seqs,
            positions: region_positions,
            ..self
        }
    }
}

pub fn align_sbr_and_ref_to_cs_worker(
    recv: crossbeam::channel::Receiver<SubreadsAndSmc>,
    sender: Sender<MsaResult>,
    ref_aligner: &Aligner<Built>,
    ref_seq: &str,
    target_idx: &HashMap<String, (usize, usize)>,
    map_params: &params::MapParams,
    align_params: &params::AlignParams,
    oup_params: &params::OupParams,
    reporter: Arc<Mutex<Reporter>>,
    mut ref_range: Option<Range<usize>>,
) {
    let mut scoped_timer = ScopedTimer::new();

    let mut max_time = 0;
    let mut max_time_qname = String::new();

    let mut inp_sbrs = 0;
    let mut out_sbrs = 0;
    let mut out_smc = 0;
    let mut fallback_num = 0;
    let mut fallback_rescued_num = 0;
    let mut channel_filter_by_no_align = 0;
    let mut channel_filter_by_cs2ref_align_fail = 0;

    for mut subreads_and_smc in recv {
        // tracing::info!("sbr_cnt:{}-{}", subreads_and_smc.smc.name, subreads_and_smc.subreads.len());
        let timer = scoped_timer.perform_timing();
        let start = Instant::now();

        let cs2ref_aln_res = align_cs_to_ref(
            &subreads_and_smc.smc.seq,
            ref_aligner,
            ref_seq,
            &subreads_and_smc.smc.name,
            ref_range.as_mut(),
        );

        if cs2ref_aln_res.is_none() {
            channel_filter_by_cs2ref_align_fail += 1;
            continue;
        }
        let mut cs2ref_aln_res = cs2ref_aln_res.unwrap();

        // 为了输出的结果更直观，这里将 smc 的结果和 reference的方向进行对齐
        if cs2ref_aln_res.is_reverse {
            subreads_and_smc.smc.seq =
                String::from_utf8(reverse_complement(subreads_and_smc.smc.seq.as_bytes())).unwrap();
        }

        // if cs2ref_aln_res.identity >= 0.9999 {
        //     continue;
        // }

        let ref_read_info = crate::read_info::ReadInfo::new_fa_record(
            "ref/-1".to_string(),
            cs2ref_aln_res.ref_sub_seq.take().unwrap(),
        );
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
            None,
            ref_range.as_ref(),
        );
        if align_res.is_none() {
            channel_filter_by_no_align += 1;
            continue;
        }
        let mut align_res = align_res.unwrap();
        align_res.identity = cs2ref_aln_res.identity;
        align_res.query_aln_len = cs2ref_aln_res.query_align_len;
        align_res = align_res.extract_error_region();

        timer.done_with_cnt(1);
        out_smc += 1;
        sender.send(align_res).unwrap();
    }

    let filter_by_alignment = inp_sbrs - out_sbrs;
    {
        let mut reporter_ = reporter.lock().unwrap();
        reporter_.sbr_reporter.filter_by_alignment += filter_by_alignment;
        reporter_.channel_reporter.out_num += out_smc;
        reporter_.channel_reporter.filter_by_no_align += channel_filter_by_no_align;
        reporter_.channel_reporter.filter_by_cs2ref_align_fail +=
            channel_filter_by_cs2ref_align_fail;
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

#[tracing::instrument(skip(records, ref_seq, ref_name, qual, low_q))]
pub fn build_msa_result_from_records(
    mut records: Vec<Record>,
    ref_seq: &str,
    ref_name: &str,
    qual: Option<&[u8]>,
    low_q: Option<u8>,
    ref_range: Option<&Range<usize>>,
) -> Option<MsaResult> {
    if records.len() == 0 {
        tracing::warn!("no align records for msa. {}", ref_name);
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

    records = records.into_iter().take(20).collect();

    let first_record = records.first().unwrap();
    let first_record = BamRecordExt::new(first_record);
    if !first_record.get_qname().starts_with("ref") {
        tracing::warn!(
            "first records not start with ref. ref_name={}, firstRecordName:{}",
            ref_name,
            first_record.get_qname()
        );
        return None;
    }

    let strands = records
        .iter()
        .map(|record| {
            if record.is_reverse() {
                "rev".to_string()
            } else {
                "fwd".to_string()
            }
        })
        .collect::<Vec<String>>();

    // 基于 ref-range 构建出来 smc range

    let smc_ranges = ref_range.map(|v| get_ref_range_from_query_range(&records[0], v));

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

    let mut major_positions = vec![-1; tot_len];
    major_pos_ins_vec.iter().for_each(|&(major, _)| {
        major_positions[*major_pos2major_starting_point.get(&major).unwrap()] = major as i32
    });
    (0..tot_len).for_each(|idx| {
        if idx > 0 && major_positions[idx] == -1 {
            major_positions[idx] = major_positions[idx - 1];
        }
    });

    let mut anchor_aligned = vec!['-' as u8; tot_len];
    major_pos_ins_vec.iter().for_each(|(major, _)| {
        anchor_aligned[*major_pos2major_starting_point.get(major).unwrap()] = ref_seq[*major];
    });
    let anchor_aligned = String::from_utf8(anchor_aligned).unwrap();

    let mut anchor_qual = vec!['-'; tot_len];
    if let Some(qual) = qual {
        major_pos_ins_vec.iter().for_each(|(major, _)| {
            anchor_qual[*major_pos2major_starting_point.get(major).unwrap()] =
                ((qual[*major] / 5u8).min(9) + '0' as u8) as char;
        });
    }
    let anchor_qual = anchor_qual.iter().collect::<String>();

    let msa = records
        .iter()
        .map(|record| {
            let mut aligned = vec!['-' as u8; tot_len];
            build_one_record_of_msa(record, &major_pos2major_starting_point, &mut aligned);
            String::from_utf8(aligned).unwrap()
        })
        .collect::<Vec<_>>();

    let mut msa_seqs = vec![anchor_qual, anchor_aligned];
    msa_seqs.extend(msa.into_iter());
    let mut names = vec!["qual".to_string(), ref_name.to_string()];
    names.extend(
        records
            .iter()
            .map(|record| BamRecordExt::new(record).get_qname()),
    );

    Some(MsaResult {
        identity: 0.0,
        query_aln_len: 0,
        mm: 0,
        ins: 0,
        del: 0,
        msa_seqs: msa_seqs,
        names: names,
        positions: major_positions,
        low_q: low_q,
        strands: strands,
        smc_ranges: smc_ranges,
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

fn merge_intervals(mut intervals: Vec<(usize, usize)>) -> Vec<(usize, usize)> {
    if intervals.is_empty() {
        return vec![];
    }
    intervals.sort_by_key(|item| item.0);

    let mut merged = Vec::new();
    let mut current = intervals[0].clone(); // 初始化为第一个区间

    for interval in intervals.into_iter().skip(1) {
        if interval.0 <= current.1 {
            // 如果当前区间的 start 小于等于上一个区间的 end，进行合并
            current.1 = current.1.max(interval.1); // 更新 end 为两者中的较大值
        } else {
            // 没有重叠，保存当前区间，更新 current 为新的区间
            merged.push(current);
            current = interval;
        }
    }

    // 最后将最后一个区间添加到 merged 中
    merged.push(current);

    merged
}

pub fn get_ref_range_from_query_range(record: &Record, range: &Range<usize>) -> Range<usize> {
    let mut qpos_cursor = None;
    let mut rpos_cursor = None;
    let mut ref_ranges = vec![];
    for [qpos, rpos] in record.aligned_pairs_full() {
        if qpos.is_some() {
            qpos_cursor = qpos;
        }
        if rpos.is_some() {
            rpos_cursor = rpos;
        }

        if qpos_cursor.is_none() || rpos_cursor.is_none() {
            continue;
        }

        if range.within_range(qpos_cursor.map(|v| v as usize).unwrap()) {
            if ref_ranges.is_empty() {
                ref_ranges.push((
                    rpos_cursor.map(|v| v as usize).unwrap(),
                    rpos_cursor.map(|v| v as usize).unwrap(),
                ));
            } else {
                let last = ref_ranges.last_mut().unwrap();
                if (last.1 + 1) == rpos_cursor.map(|v| v as usize).unwrap() {
                    last.1 += 1;
                } else {
                    ref_ranges.push((
                        rpos_cursor.map(|v| v as usize).unwrap(),
                        rpos_cursor.map(|v| v as usize).unwrap(),
                    ));
                }
            }
        }
    }

    Range::new_range(ref_ranges)
}

#[cfg(test)]
mod test {
    use std::collections::HashMap;

    use crate::{
        align_sbr_to_smc,
        params::{AlignParams, MapParams, OupParams},
        sbr_and_ref_to_cs::build_msa_result_from_records,
        SubreadsAndSmc,
    };

    #[test]
    fn merge_intervals() {
        let intervals = vec![(0, 65), (64, 146), (145, 355)];
        println!("{:?}", super::merge_intervals(intervals));
    }

    #[test]
    fn test_msa() {
        // 20260311_240601Y0012_Run0003/181179
        // 1041 CCCCTT
        let smc_seq = "ACAGATCTTTCCCTCGGGCATTCTCAAGACGGATCCCCATTTCCAGACGATAAGGCTGCATTAAATCGAGCGGGCGGAGTACGCCATACAAGCCGGAAGCATTCGCAAATGCTGTTGGCAAAATCGAAATCGTCTTCGCTGAAGGTTTCGGCCTGCAAGCCGGTGTAGACATCACCTTTAAACGCCAGAATCGCCTGGCGGGCATTCGCCGGCGTGAAATCTGGCTGCCAGTCATGAAAGCGAGCGGCGTTGATACCCGCCAGTTTGTCGCTGATGCGCATCAGCGTGCTAATCTGCGGAGGCGTCAGTTTCCGCGCCTCATGGATCAACTGCTGGGAATTGTCTAACAGCTCCGGCAGCGTATAGCGCGTGGTGGTCAACGGGCTTTGGTAATCAAGCGTTTTCGCAGGTGAAATAAGAATCAGCATATCCAGTCCTTGCAGGAAATTTATGCCGACTTTAGCAAAAAATGAGAATGAGTTGATCGATAGTTGTGATTACTCCTGCGAAACATCATCCCACGCGTCCGGAGAAAGCTGGCGACCGATATCCGGATAACGCAATGGATCAAACACCGGGCGCACGCCGAGTTTACGCTGGCGTAGATAATCACTGGCAATGGTATGAACCACAGGCGAGAGCAGTAAAATGGCGGTCAAATTGGTAATAGCCATGCAGGCCATTATGATATCTGCCAGTTGCCACATCAGCGGAAGGCTTAGCAAGGTGCCGCCGATGACCGTTGCGAAGGTGCAGATCCGCAAACACCAGATCGCTTTAGGGTTGTTCAGGCGTAAAAAGAAGAGATTGTTTTCGGCATAAATGTAGTTGGCAACGATGGAGCTGAAGGCAAACAGAATAACCCACAAGGGTAACAAACTCAGCACCCCAGGAACCCATTAGCACCCGCATCGCCTTCTGGATAAGCTGAATACCTTCCAGCGGCATGTAGGTTGTGCCGTTACCCGCCAGTAATATCAGCATGGCGCTTGCCGTACAGATGACCAGGGTGTCGATAAAAATGCCAATCATCTGGACAATCCCCTTGCGCTGCCGGATGCGGAGGCCAGGACGCCGCTGCCGCTGCCGCGTTTGGCGTCGAACCCATTCCCGCCTCATTGGAAAACATACTGCGCTGAAAACCGTTAGTAATCGCCTGGCTTAAGGTATATCCGCCGCGCCGCCTGCCGCTTCCTGCCAGCCAAAAGCACTCTCAAAAATAGACCAAATGACGTGGGGAAGTTGCCCGATATTCATTACGCAAATTACCAGGCTGGTCAGTACCCAGATTATCGCCATCAACGGGACAAAGCCCTGCATGAGCCGGGCGACGCCATGAAGACGCGAGTGATTGCCAGCAGAGTAAAGACAGCGAGAATAATGCCTGTCACCAGCGGGGGAAAATCAAAAGAAAAACTCAGGGCGCGGGCAACGGCGTTCGCTTGAACTCCGCTGAAAATTATGCCATAGGCGATGAGCAAAAAGACGGCGAACAGAACGCCCATCCAGCGCATCCCCAGCCCGCGCGCCATATACCATGCCGGTCCGCCACGAAACTGCCCATTGACGTCACGTTCTTTATAAAGTTGTGCCAGAGAACATTCGGCAAACGAGGTCGCCATGCCGATAAACGCGGCAACCCACATCCAAAAGACGGCTCCAGGTCCACCGCGCGGTAATAGCCAGCGCAACGCCGGCCAGGTTGCCGCTACCCACGCGCGCCGCAAGACTGGTACACAATGACTGAAATGAGGTTAAACCGCCTGGCTGTGGATGAATGCTATTTTTAAGACTTTTGCCAAAC";
        let qual = vec![5; smc_seq.len()];

        let mut subreads_and_smc = SubreadsAndSmc::from_seq_qual(smc_seq, qual);
        // subreads_and_smc.add_subread_str("0".into(), "".to_string());
        // subreads_and_smc.add_subread_str("1".into(), "".to_string());
        // subreads_and_smc.add_subread_str("2".into(), "".to_string());
        // subreads_and_smc.add_subread_str("3".into(), "".to_string());
        // subreads_and_smc.add_subread_str("4".into(), "".to_string());
        subreads_and_smc.add_subread_str("ref/-1".into(), "ACAGATCTTTCCCTCGGGCATTCTCAAGACGGATCCCCATTTCCAGACGATAAGGCTGCATTAAATCGAGCGGGCGGAGTACGCCATACAAGCCGGAAGCATTCGCAAATGCTGTTGGCAAAATCGAAATCGTCTTCGCTGAAGGTTTCGGCCTGCAAGCCGGTGTAGACATCACCTTTAAACGCCAGAATCGCCTGGCGGGCATTCGCCGGCGTGAAATCTGGCTGCCAGTCATGAAAGCGAGCGGCGTTGATACCCGCCAGTTTGTCGCTGATGCGCATCAGCGTGCTAATCTGCGGAGGCGTCAGTTTCCGCGCCTCATGGATCAACTGCTGGGAATTGTCTAACAGCTCCGGCAGCGTATAGCGCGTGGTGGTCAACGGGCTTTGGTAATCAAGCGTTTTCGCAGGTGAAATAAGAATCAGCATATCCAGTCCTTGCAGGAAATTTATGCCGACTTTAGCAAAAAATGAGAATGAGTTGATCGATAGTTGTGATTACTCCTGCGAAACATCATCCCACGCGTCCGGAGAAAGCTGGCGACCGATATCCGGATAACGCAATGGATCAAACACCGGGCGCACGCCGAGTTTACGCTGGCGTAGATAATCACTGGCAATGGTATGAACCACAGGCGAGAGCAGTAAAATGGCGGTCAAATTGGTAATAGCCATGCAGGCCATTATGATATCTGCCAGTTGCCACATCAGCGGAAGGCTTAGCAAGGTGCCGCCGATGACCGTTGCGAAGGTGCAGATCCGCAAACACCAGATCGCTTTAGGGTTGTTCAGGCGTAAAAAGAAGAGATTGTTTTCGGCATAAATGTAGTTGGCAACGATGGAGCTGAAGGCAAACAGAATAACCCACAAGGGTAACAAACTCAGCACCCCAGGAACCCATTAGCACCCGCATCGCCTTCTGGATAAGCTGAATACCTTCCAGCGGCATGTAGGTTGTGCCGTTACCCGCCAGTAATATCAGCATGGCGCTTGCCGTACAGATGACCAGGGTGTCGATAAAAATGCCAATCATCTGGACAATCCCCTTGCGCTGCCGGATGCGGAGGCCAGGACGCCGCTGCCGCTGCCGCGTTTGGCGTCGAACCCATTCCCGCCTCATTGGAAAACATACTGCGCTGAAAACCGTTAGTAATCGCCTGGCTTAAGGTATATCCGCCGCGCCGCCTGCCGCTTCCTGCCAGCCAAAAGCACTCTCAAAAATAGACCAAATGACGTGGGGAAGTTGCCCGATATTCATTACGCAAATTACCAGGCTGGTCAGTACCCAGATTATCGCCATCAACGGGACAAAGCCCTGCATGAGCCGGGCGACGCCATGAAGACGCGAGTGATTGCCAGCAGAGTAAAGACAGCGAGAATAATGCCTGTCACCAGCGGGGGAAAATCAAAAGAAAAACTCAGGGCGCGGGCAACGGCGTTCGCTTGAACTCCGCTGAAAATTATGCCATAGGCGATGAGCAAAAAGACGGCGAACAGAACGCCCATCCAGCGCATCCCCAGCCCGCGCGCCATATACCATGCCGGTCCGCCACGAAACTGCCCATTGACGTCACGTTCTTTATAAAGTTGTGCCAGAGAACATTCGGCAAACGAGGTCGCCATGCCGATAAACGCGGCAACCCACATCCAAAAGACGGCTCCAGGTCCACCGCGCGGTAATAGCCAGCGCAACGCCGGCCAGGTTGCCGCTACCCACGCGCGCCGCAAGACTGGTACACAATGACTGAAATGAGGTTAAACCGCCTGGCTGTGGATGAATGCTATTTTTAAGACTTTTGCCAAAC".to_string());
        subreads_and_smc.add_subread_str("2026/sbr/5".into(), "GTTTGGCAAAAGTCTTAAAAATAGCATTCATCCACACCCAGCGGTTTACCCATTCAGTCATTGTGTACCAGTCTTGCGGCGCGGCTTGGGTAGCGGCAACCTGGCCGGCGTTGGCGCGGCTATTACCGCGGTGGCCTGTGAGCCGTCTTTTGATGTGGGATTGCCGGCGTTTATCGCATGGCGACCTCGTTTGCCAAGTTCTCTGGCACAACTTTTATAAAGAACGTGACGTCCAATGGGCAGTTTCGTGGGCGGACCGCATGGTATATGGCGCGCGGGCTGGGGGATGCGCTGGATGGGCGTTCATGTTCGCCGTCTTTTTGCTCATCGCCTATGGCATAATTTTCAGCGAGTTCAGCGAACGCCGTTGCCGCGCCCTGAGTTTCTTTTGATTTTCCCCGCTGGTGACAGCATTATTCTCGATGTCTTTACTCTGCTGAATCACATCGCGTCTTCATGTGGCGTCGCCCGAGTCCTCATGCAGGGCTTTGTCCCGTTGATGGAGCTAATCTGGGTACTGACCAGATCTGGTAATTTGCGTATAGAATATGGTGCAACTTCCCCACGCCATTTGGTCTATTTTTTGGATGTTTTTTGGGCTGGCAGGAAGCGGCAGGCGGCGCGGCGGGATATACCTTAAAGCAGGCGATTACTAACGGTTTTTTTCAGCCAGTATGTTTTTCCAATGAAGCGGGAGGGTTTCGAGCCAAACGCGCAGCGGCGTGGGGGGGGGCGGCTCCCTCGCCTCCGCATCCGCAGCGCAAAGGATTGTCCAGAAGATTGCATTTTATCGACACCCTGGTCATCTGTACAGGAAGCGCCGCTGATATTACGCGGTAACGGCACAACCTACAGCGCTGAGGTATTCAGCTTATCCAGAAGGCGATGCGGGTGCTAATGGGTTCCGTGTGGTGCTGAGTTTGTTACCCTTGTGGTTATCTTTTGCCTTCACCCATCGTGCCAACTTACATTGATGCGAAAAAATTCTTCTTTTTTACGCCTGAAACAACCTAAAGCGATCTGGTGTTTGCAAAGGATCTGCACCTTTCGCAACAGGTCATCGCGGCACCTTGCTAAGCCTTCCGCGATGTGGCAACTGGCAGATATCATAATGGCTGCATGGCTATACCGATTTGACGCCATTTTACTGCTCTCGCCTTGTGGTTCTACCATTGCCAAGTGATTATCTACAGCCAGCGTAAACTCGGCGTGCGCCGGTGTTTTGATCATTGCGTTATCCCGGAATCGGTCGCCCAGCTTTCTCCGGACGCGTGGGATGATGTTTCGCAGGTAATCACAATCAACTCCATTCTCATTTTTGCCTAAAGTCCGGCCATAATTTTCCTGCAAGGACTGGATATGCTGATTTTAATTTCACCTGCGAAACGCTTGATTACCAAAAGCCTTGACCATCCACGCCGTTATACGCTGCCGAGCTGTTAGAAATTCGCAGCAGTTGATCCCAGAGGTCCGGAAACTGACGCCTCCCGCAGATTAGCACGTTTTTGGGATGCGCTAGCGACAACTGGCGGGGTATCACGCCGCTCGCTTTCATTGGACTGGCAGCCAGATTTCACGCCGGCGATGCCCGCCGGCGATTCTGGCGTTAAAGGTGATGTCTACCACCAGGCTTGCAGCCAAACCTTCAGGAAGACGATTTCGATTTTGCCCAACAGCATTTGCGAATGCTTCCGCTTGGTAGGCTTACTCCGCCCTCCGATTTAATGCCTTTATCGTCCTGGAAATGGGATCCGGTCTTGAGAATGCCCGAGGGAAAGATCTGT".to_string());

        let mut target_idx = HashMap::new();
        target_idx.insert("smc".to_string(), (0, smc_seq.len()));
        let mut align_params = AlignParams::new()
            .set_m_score(4)
            .set_mm_score(10)
            .set_gap_open_penalty("4,48".to_string())
            .set_gap_extension_penalty("2,1".to_string())
            .set_poly_n_gap_left_align(false);

        let oup_params = OupParams::new();

        let (align_infos, no_hit_indices) = align_sbr_to_smc(
            &subreads_and_smc,
            (0..subreads_and_smc.subreads.len()).collect(),
            &target_idx,
            &MapParams {},
            &align_params,
            &oup_params,
            false,
        );

        println!("{}", align_infos.len());

        let align_res = build_msa_result_from_records(
            align_infos,
            &subreads_and_smc.smc.seq,
            &subreads_and_smc.smc.name,
            subreads_and_smc.smc.qual.as_deref(),
            Some(10),
            None,
        );

        let mut align_res = align_res.unwrap();
        align_res.identity = 1.0;
        align_res.query_aln_len = 100;
        // align_res = align_res.extract_error_region();

        println!("{}", align_res.msa_seqs.join("\n"));
    }
}
