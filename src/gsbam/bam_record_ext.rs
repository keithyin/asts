//! `BamRecordExt`：对 noodles `RecordBuf` 的借用式扩展包装。
//!
//! 语义与 `gskits::gsbam::bam_record_ext::BamRecordExt` 逐一对齐：
//! - `get_int`：任意整型 aux → u32（有符号按 `as u32` 位截断），类型不符 → None
//! - `get_float`：Float → f32，其余 → None
//! - `get_uint_list`：UInt8/UInt16/UInt32 数组 → Vec<u32>，其余 → None
//! - `compute_identity`：gap-compressed，eqx CIGAR 专用，遇 `M` panic
//! - `compute_query_coverage`：(=,X,M,I) / 全 query 长度（含 soft clip）
//! - `reference_start/end`：CIGAR 推导的 0-based 半开区间（与原实现一致，
//!   不依赖 noodles 的 `alignment_end`，便于在写出前使用）

use noodles_sam as sam;
use sam::alignment::record::cigar::op::Kind;
use sam::alignment::record::data::field::Tag;
use sam::alignment::record_buf::data::field::value::Array;
use sam::alignment::record_buf::{data::field::Value, RecordBuf};

pub type BamRecord = RecordBuf;

pub struct BamRecordExt<'a> {
    bam_record: &'a BamRecord,
}

impl<'a> BamRecordExt<'a> {
    pub fn new(bam_record: &'a BamRecord) -> Self {
        Self { bam_record }
    }

    /// 存储的 query 序列（原始字节，不做碱基校验）。
    pub fn get_seq(&self) -> String {
        unsafe {
            String::from_utf8_unchecked(self.bam_record.sequence().as_ref().to_vec())
        }
    }

    /// query name（缺失时返回空串，与 htslib 行为一致）。
    pub fn get_qname(&self) -> String {
        unsafe {
            String::from_utf8_unchecked(
                self.bam_record
                    .name()
                    .unwrap_or_default()
                    .to_vec(),
            )
        }
    }

    /// phred 质量值（无 +33 偏移）。
    pub fn get_qual(&self) -> Vec<u8> {
        self.bam_record.quality_scores().as_ref().to_vec()
    }

    fn get_value(&self, tag: Tag) -> Option<&Value> {
        self.bam_record.data().get(&tag)
    }

    /// 任意整型 aux → u32。有符号值按 `as u32` 位转换（与原实现一致）。
    pub fn get_int(&self, tag: &[u8]) -> Option<u32> {
        let tag = Tag::new(tag[0], tag[1]);
        self.get_value(tag).and_then(|value| match value {
            Value::Int8(v) => Some(*v as u32),
            Value::UInt8(v) => Some(*v as u32),
            Value::Int16(v) => Some(*v as u32),
            Value::UInt16(v) => Some(*v as u32),
            Value::Int32(v) => Some(*v as u32),
            Value::UInt32(v) => Some(*v as u32),
            _ => None,
        })
    }

    /// Float aux → f32。
    pub fn get_float(&self, tag: &[u8]) -> Option<f32> {
        let tag = Tag::new(tag[0], tag[1]);
        self.get_value(tag).and_then(|value| match value {
            Value::Float(v) => Some(*v),
            _ => None,
        })
    }

    /// UInt8/UInt16/UInt32 数组 → Vec<u32>。
    pub fn get_uint_list(&self, tag: &[u8]) -> Option<Vec<u32>> {
        let tag = Tag::new(tag[0], tag[1]);
        self.get_value(tag).and_then(|value| match value {
            Value::Array(Array::UInt8(values)) => {
                Some(values.iter().map(|&v| v as u32).collect())
            }
            Value::Array(Array::UInt16(values)) => {
                Some(values.iter().map(|&v| v as u32).collect())
            }
            Value::Array(Array::UInt32(values)) => {
                Some(values.iter().copied().collect())
            }
            _ => None,
        })
    }

    pub fn get_cx(&self) -> Option<u8> {
        self.get_int(b"cx").map(|v| v as u8)
    }

    pub fn get_ch(&self) -> Option<u32> {
        self.get_int(b"ch")
    }

    pub fn get_np(&self) -> Option<u32> {
        self.get_int(b"np")
    }

    pub fn get_rq(&self) -> Option<f32> {
        self.get_float(b"rq")
    }

    #[allow(dead_code)]
    pub fn get_coverage(&self) -> Option<f32> {
        self.get_float(b"ec")
    }

    #[allow(dead_code)]
    pub fn get_identity(&self) -> Option<f32> {
        self.get_float(b"iy")
    }

    pub fn get_dw(&self) -> Option<Vec<u32>> {
        self.get_uint_list(b"dw")
    }

    pub fn get_ar(&self) -> Option<Vec<u32>> {
        self.get_uint_list(b"ar")
    }

    pub fn get_cr(&self) -> Option<Vec<u32>> {
        self.get_uint_list(b"cr")
    }

    pub fn get_nn(&self) -> Option<Vec<u32>> {
        self.get_uint_list(b"nn")
    }

    pub fn get_be(&self) -> Option<Vec<u32>> {
        self.get_uint_list(b"be")
    }

    /// gap-compressed identity：eq / (eq + X + I + D)。soft clip 与 RefSkip 不计入。
    /// 只支持 eqx CIGAR（`M` 直接 panic，与原实现一致）。
    pub fn compute_identity(&self) -> f32 {
        let mut aligned_span = 0usize;
        let mut matched = 0usize;
        self.bam_record.cigar().as_ref().iter().for_each(|op| match op.kind() {
            Kind::SequenceMatch => {
                matched += op.len();
                aligned_span += op.len();
            }
            Kind::SequenceMismatch | Kind::Insertion | Kind::Deletion => {
                aligned_span += op.len();
            }
            Kind::Match => panic!("must use eqx cigar"),
            _ => {}
        });
        let aligned_span = if aligned_span > 0 { aligned_span } else { 1 };
        matched as f32 / aligned_span as f32
    }

    /// (=,X,M,I) / 全 query 长度。分母与 rust-htslib 的
    /// `seq_len_from_cigar(include_hard_clip=true)` 一致：M/I/S/=/X 恒计入，
    /// H 也计入。
    pub fn compute_query_coverage(&self) -> f32 {
        let seq_len = self.seq_len_from_cigar(true);
        let aligned_span = self
            .bam_record
            .cigar()
            .as_ref()
            .iter()
            .map(|op| match op.kind() {
                Kind::SequenceMatch | Kind::SequenceMismatch | Kind::Match | Kind::Insertion => {
                    op.len()
                }
                _ => 0,
            })
            .sum::<usize>();
        let seq_len = if seq_len > 0 { seq_len } else { 1 };
        aligned_span as f32 / seq_len as f32
    }

    /// query 长度（CIGAR 推导），语义对齐 rust-htslib
    /// `BamRecordExtensions::seq_len_from_cigar(include_hard_clip)`：
    /// M/I/S/=/X 恒计入；H 仅在 `include_hard_clip` 为 true 时计入。
    fn seq_len_from_cigar(&self, include_hard_clip: bool) -> usize {
        self.bam_record
            .cigar()
            .as_ref()
            .iter()
            .map(|op| match op.kind() {
                Kind::SequenceMatch | Kind::SequenceMismatch | Kind::Match | Kind::Insertion | Kind::SoftClip => {
                    op.len()
                }
                Kind::HardClip if include_hard_clip => op.len(),
                _ => 0,
            })
            .sum::<usize>()
    }

    /// query 坐标下比对区间起点（leading soft clip 长度，否则 0）。
    pub fn query_alignment_start(&self) -> usize {
        match self.bam_record.cigar().as_ref().iter().next() {
            Some(op) if op.kind() == Kind::SoftClip => op.len(),
            _ => 0,
        }
    }

    /// query 坐标下比对区间终点（exclusive）：`query_alignment_start + (=,X,M,I) 长度`。
    pub fn query_alignment_end(&self) -> usize {
        let aligned_qlen = self
            .bam_record
            .cigar()
            .as_ref()
            .iter()
            .map(|op| match op.kind() {
                Kind::SequenceMatch | Kind::SequenceMismatch | Kind::Insertion | Kind::Match => {
                    op.len()
                }
                _ => 0,
            })
            .sum::<usize>();
        self.query_alignment_start() + aligned_qlen
    }

    /// 0-based 比对起点（未 set pos 时 panic，与原实现一致）。
    pub fn reference_start(&self) -> usize {
        let pos = self
            .bam_record
            .alignment_start()
            .expect("set_pos first");
        usize::from(pos) - 1
    }

    /// 0-based 比对终点（exclusive）= pos + (=,X,M,D,N) 长度。
    pub fn reference_end(&self) -> usize {
        let aligned_rlen = self
            .bam_record
            .cigar()
            .as_ref()
            .iter()
            .map(|op| match op.kind() {
                Kind::SequenceMatch
                | Kind::SequenceMismatch
                | Kind::Deletion
                | Kind::Match
                | Kind::Skip => op.len(),
                _ => 0,
            })
            .sum::<usize>();
        self.reference_start() + aligned_rlen
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles_core::Position;
    use noodles_sam::alignment::record::cigar::op::Op;
    use noodles_sam::alignment::record_buf::Cigar;

    fn record_with(pos: usize, ops: &[(Kind, usize)]) -> BamRecord {
        RecordBuf::builder()
            .set_alignment_start(Position::try_from(pos).unwrap())
            .set_cigar(
                ops.iter()
                    .map(|&(kind, len)| Op::new(kind, len))
                    .collect::<Cigar>(),
            )
            .build()
    }

    /// coverage 分母与 rust-htslib `seq_len_from_cigar(true)` 一致：S 与 H 都计入。
    /// "5H2S3=1X2I1D2S" → 分母 5+2+3+1+2+2=15，分子 3+1+2=6。
    #[test]
    fn test_compute_query_coverage_with_hard_clip() {
        let rec = record_with(
            1,
            &[
                (Kind::HardClip, 5),
                (Kind::SoftClip, 2),
                (Kind::SequenceMatch, 3),
                (Kind::SequenceMismatch, 1),
                (Kind::Insertion, 2),
                (Kind::Deletion, 1),
                (Kind::SoftClip, 2),
            ],
        );
        let ext = BamRecordExt::new(&rec);
        // 与 gskits 原版一致：query_alignment_start 只看第一个 CIGAR op，
        // 首个 op 是 H 时返回 0（即使后面还有 S）
        assert_eq!(ext.query_alignment_start(), 0);
        assert_eq!(ext.query_alignment_end(), 6);
        assert_eq!(ext.reference_start(), 0);
        assert_eq!(ext.reference_end(), 5);
        let cov = ext.compute_query_coverage();
        assert_eq!(cov, 6.0 / 15.0);
    }

    /// identity：eq/(eq+X+I+D)，S/H/N 不计入分母。
    #[test]
    fn test_compute_identity() {
        let rec = record_with(
            1,
            &[
                (Kind::SoftClip, 3),
                (Kind::SequenceMatch, 4),
                (Kind::SequenceMismatch, 1),
                (Kind::Insertion, 2),
                (Kind::Deletion, 1),
                (Kind::Skip, 10),
            ],
        );
        let ext = BamRecordExt::new(&rec);
        // 分母 4+1+2+1=8
        assert_eq!(ext.compute_identity(), 4.0 / 8.0);
    }
}
