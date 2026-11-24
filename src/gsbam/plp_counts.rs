//! `compute_max_ins_of_each_ref_position` 移植（原
//! `gskits::gsbam::plp_counts_from_records`，逐行对齐，blacklist 机制不再保留——
//! asts 调用时恒传 `None`）。

use std::{cmp, collections::HashMap};

use super::{aligned_pairs::aligned_pairs_full, bam_record_ext::{BamRecord, BamRecordExt}};

/// 每个被覆盖的参考位置上，紧随其后的最大插入碱基数（gap-compressed）。
/// 返回 map 包含所有被覆盖的 ref pos（值可能为 0）。
pub fn compute_max_ins_of_each_ref_position(
    records: &[BamRecord],
    rstart: Option<usize>,
    rend: Option<usize>,
) -> HashMap<i64, i32> {
    let mut pos2ins = HashMap::new();

    let rstart = rstart.map(|v| v as i64);
    let rend = rend.map(|v| v as i64);
    for record in records {
        let record_ext = BamRecordExt::new(record);
        let mut start = rstart.unwrap_or(record_ext.reference_start() as i64);
        let mut end = rend.unwrap_or(record_ext.reference_end() as i64);
        start = cmp::max(start, record_ext.reference_start() as i64);
        end = cmp::min(end, record_ext.reference_end() as i64);

        let mut rpos_cursor = None;
        let mut cur_ins = 0i32;
        let query_end = record_ext.query_alignment_end();

        let mut aligned_pair_full = aligned_pairs_full(record);
        if aligned_pair_full.is_empty() {
            continue;
        }

        // this make the following for loop correct.
        // if the query match to the last base of the ref seqence. the following for loop won't give the right result
        // but add this , it will get the right result.
        if let Some(last_q_pos) = aligned_pair_full.last().unwrap()[0] {
            aligned_pair_full.push([Some(last_q_pos + 1), None]);
        }

        for [qpos, rpos] in aligned_pair_full.into_iter() {
            if rpos.is_some() {
                rpos_cursor = rpos;
            }
            if rpos_cursor.is_none() {
                continue;
            }

            if rpos_cursor.unwrap() < start {
                continue;
            }

            if rpos_cursor.unwrap() >= end {
                // set the last rpos max ins and then break
                let rpos_ = rpos_cursor.unwrap();
                pos2ins.entry(rpos_ - 1).or_insert(0);
                *pos2ins.get_mut(&(rpos_ - 1)).unwrap() =
                    cmp::max(*pos2ins.get(&(rpos_ - 1)).unwrap(), cur_ins);
                break;
            }

            if let Some(qpos_) = qpos {
                if qpos_ as usize >= query_end {
                    // set the last rpos max ins and then break
                    let rpos_ = rpos_cursor.unwrap();

                    pos2ins.entry(rpos_).or_insert(0);
                    *pos2ins.get_mut(&rpos_).unwrap() =
                        cmp::max(*pos2ins.get(&rpos_).unwrap(), cur_ins);
                    break;
                }
            }

            if let Some(rpos_) = rpos {
                if rpos_ > start {
                    pos2ins.entry(rpos_ - 1).or_insert(0);
                    *pos2ins.get_mut(&(rpos_ - 1)).unwrap() =
                        cmp::max(*pos2ins.get(&(rpos_ - 1)).unwrap(), cur_ins);
                }

                cur_ins = 0;
            } else {
                cur_ins += 1;
            }
        }
    }

    pos2ins
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gsbam::aligned_pairs;
    use noodles_sam::alignment::record_buf::RecordBuf;
    use noodles_sam::alignment::record::cigar::op::{Kind, Op};
    use noodles_sam::alignment::record_buf::Cigar;
    use noodles_core::Position;

    fn record_with(pos: i64, cigar: &[(&Kind, usize)]) -> RecordBuf {
        let mut r = RecordBuf::builder()
            .set_alignment_start(Position::try_from(pos as usize + 1).unwrap())
            .set_cigar(
                cigar
                    .iter()
                    .map(|&(kind, len)| Op::new(*kind, len))
                    .collect::<Cigar>(),
            )
            .build();
        *r.flags_mut() = noodles_sam::alignment::record::Flags::empty();
        r
    }

    /// 对齐 gskits plp_counts_from_records.rs 的语义测试：
    /// `"2S1=1I2=1D1=3S"` @ pos 0 → pos0 的 max_ins=1，其余覆盖位为 0，共 5 个 key。
    #[test]
    fn test_max_ins() {
        let rec = record_with(
            0,
            &[
                (&Kind::SoftClip, 2),
                (&Kind::SequenceMatch, 1),
                (&Kind::Insertion, 1),
                (&Kind::SequenceMatch, 2),
                (&Kind::Deletion, 1),
                (&Kind::SequenceMatch, 1),
                (&Kind::SoftClip, 3),
            ],
        );
        let res = compute_max_ins_of_each_ref_position(&vec![rec], None, None);
        assert_eq!(res.len(), 5);
        assert_eq!(res[&0], 1);
        for k in 1..5 {
            assert_eq!(res[&k], 0);
        }
    }

    /// `"4=3S"` @ pos 0 → 4 个位置，max_ins 全 0（右侧 soft clip 不计）。
    #[test]
    fn test_max_ins_right_softclip() {
        let rec = record_with(
            0,
            &[(&Kind::SequenceMatch, 4), (&Kind::SoftClip, 3)],
        );
        let res = compute_max_ins_of_each_ref_position(&vec![rec], None, None);
        assert_eq!(res.len(), 4);
        assert!(res.values().all(|&v| v == 0));
    }

    /// aligned_pairs_full 与 rust-htslib `IterAlignedPairsFull` 的语义对齐：
    /// soft clip 产出 `[Some(q), None]` 并推进 qpos；hard clip 不产出、不推进。
    #[test]
    fn test_aligned_pairs() {
        let rec = record_with(
            2,
            &[
                (&Kind::HardClip, 5),
                (&Kind::SoftClip, 2),
                (&Kind::SequenceMatch, 2),
                (&Kind::Insertion, 1),
                (&Kind::Deletion, 2),
                (&Kind::SequenceMatch, 1),
                (&Kind::SoftClip, 1),
            ],
        );
        let pairs = aligned_pairs::aligned_pairs_full(&rec);
        assert_eq!(
            pairs,
            vec![
                // 2S
                [Some(0), None],
                [Some(1), None],
                // 2=
                [Some(2), Some(2)],
                [Some(3), Some(3)],
                // 1I
                [Some(4), None],
                // 2D
                [None, Some(4)],
                [None, Some(5)],
                // 1=
                [Some(5), Some(6)],
                // trailing 1S
                [Some(6), None],
            ]
        );
    }
}
