//! `aligned_pairs_full` 移植（替代 `rust_htslib::bam::ext::BamRecordExtensions`）。
//!
//! 与 rust-htslib 0.49 的 `IterAlignedPairsFull` 逐位对齐（这是关键：ats 的 MSA
//! 构建、cs 校验、plp 统计都直接消费这个 pair 序列，语义必须完全一致）：
//!
//! - 每个 CIGAR op 逐碱基产出 `[qpos, rpos]`（0-based，qpos 从 0 开始计）
//! - `=`/`X`/`M` → `[Some(q), Some(r)]`
//! - `I`/**`S`** → `[Some(q), None]`（soft clip 与 insertion 同等对待，会推进 qpos）
//! - `D`/`N` → `[None, Some(r)]`
//! - `H` → 不产出，也不推进 qpos
//! - `P` → panic（与原实现一致；padding 只用于 MSA，不影响 BAM 输出）
//!
//! 与 rust-htslib 唯一的差异：未设置比对起点（pos < 0）的记录返回空序列，
//! 原实现会从 -1 开始产出负坐标 pair（调用方在此之前就会因 `reference_start`
//! 断言失败，故此差异不可达）。

use noodles_sam as sam;
use sam::alignment::record::cigar::op::Kind;

use super::bam_record_ext::BamRecord;

/// rust-htslib 兼容的对类型：`[qpos, rpos]`（0-based）。
pub type AlignedPair = [Option<i64>; 2];

pub fn aligned_pairs_full(record: &BamRecord) -> Vec<AlignedPair> {
    let mut result = Vec::new();

    let mut rpos_cursor = match record.alignment_start() {
        Some(pos) => (usize::from(pos) - 1) as i64,
        None => return result,
    };
    let mut qpos_cursor: i64 = 0;

    for op in record.cigar().as_ref().iter() {
        let len = op.len();
        match op.kind() {
            Kind::SequenceMatch | Kind::SequenceMismatch | Kind::Match => {
                for _ in 0..len {
                    result.push([Some(qpos_cursor), Some(rpos_cursor)]);
                    qpos_cursor += 1;
                    rpos_cursor += 1;
                }
            }
            Kind::Insertion | Kind::SoftClip => {
                for _ in 0..len {
                    result.push([Some(qpos_cursor), None]);
                    qpos_cursor += 1;
                }
            }
            Kind::Deletion | Kind::Skip => {
                for _ in 0..len {
                    result.push([None, Some(rpos_cursor)]);
                    rpos_cursor += 1;
                }
            }
            Kind::HardClip => {}
            Kind::Pad => panic!("Padding (Cigar::Pad) is not supported."),
        }
    }

    result
}
