//! gsbam 兼容层（原 `gskits::gsbam` 的 noodles 版移植）。
//!
//! 覆盖 asts 实际使用的 gsbam API：
//! - [`bam_record_ext::BamRecordExt`]：tag getter、identity/coverage、比对边界
//! - [`aligned_pairs::aligned_pairs_full`]：等价于 rust-htslib 的
//!   `bam::ext::BamRecordExtensions::aligned_pairs_full`
//! - [`plp_counts::compute_max_ins_of_each_ref_position`]
//! - [`fastx2bam::{fasta2bam, fastq2bam}`]
//!
//! 原实现中 asts 未使用的部分（`BamHeaderExt`、`PlpCnts`、`cigar_ext`、
//! blacklist 机制、`draw_aligned_seq` 等）不再移植。

pub mod aligned_pairs;
pub mod bam_record_ext;
pub mod fastx2bam;
pub mod plp_counts;

pub use aligned_pairs::{aligned_pairs_full, AlignedPair};
pub use bam_record_ext::{BamRecord, BamRecordExt};
