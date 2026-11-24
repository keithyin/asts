use std::{borrow::Cow, collections::HashMap, hash::Hash};

use noodles_core::Position;
use noodles_sam as sam;
use sam::alignment::record::{Flags, MappingQuality};
use sam::alignment::record_buf::{data::field::Value, QualityScores, RecordBuf, Sequence};
use sam::alignment::record::data::field::Tag;
use sam::alignment::record_buf::data::field::value::Array;

use crate::dna::reverse_complement;

use crate::align::cigar::ops_to_record_cigar;
use crate::align::Hit;
use crate::read_info::ReadInfo;

pub type BamRecord = sam::alignment::RecordBuf;

pub fn build_bam_record_from_mapping<T>(
    hit: &Hit,
    query_record: &ReadInfo,
    target_idx: &HashMap<T, (usize, usize)>,
) -> BamRecord
where
    T: std::borrow::Borrow<str> + Eq + Hash,
{
    let mut seq = Cow::Borrowed(query_record.seq.as_str());
    let is_rev = hit.is_reverse;

    if is_rev {
        seq = Cow::Owned(String::from_utf8(reverse_complement(seq.as_bytes())).unwrap());
    }

    // hit.query_start, hit.query_end 是相对于原始 query 而言的(即 未 reverse 的 query 而言)
    let cigar_str = ops_to_record_cigar(
        &hit.cigar_ops,
        hit.query_start,
        hit.query_end,
        seq.len(),
        is_rev,
    );

    let qual = if let Some(ref qual_) = query_record.qual {
        // 全 255 是 BAM 的 qual 缺失标记（如 fasta2bam 的占位值）：留空交给
        // noodles encoder 自动填 0xFF×len，与 rust-htslib 写 255 的输出一致
        // （noodles 校验 qual ≤ 93，255 不能直接写入）
        if qual_.iter().all(|&v| v == 255) {
            vec![]
        } else if is_rev {
            qual_.iter().copied().rev().collect()
        } else {
            qual_.clone()
        }
    } else {
        vec![]
    };

    let mut flags = Flags::empty();
    if is_rev {
        flags |= Flags::REVERSE_COMPLEMENTED;
    }

    // flag 组合与迁移前逐位一致：
    // primary → 清 secondary/supplementary；非 primary → 置 secondary+supplementary；
    // 再按 is_supplementary 收敛
    if hit.is_primary {
        flags &= !Flags::SECONDARY;
        flags &= !Flags::SUPPLEMENTARY;
    } else {
        flags |= Flags::SECONDARY;
        flags |= Flags::SUPPLEMENTARY;
    }

    if hit.is_supplementary {
        flags &= !Flags::SECONDARY;
        flags |= Flags::SUPPLEMENTARY;
    } else {
        flags &= !Flags::SUPPLEMENTARY;
    }
    flags &= !Flags::UNMAPPED;

    let mut bam_record = RecordBuf::builder()
        .set_name(query_record.name.as_bytes())
        .set_flags(flags)
        .set_sequence(Sequence::from(seq.as_bytes().to_vec()))
        .set_quality_scores(QualityScores::from(qual))
        .set_cigar(cigar_str)
        .set_alignment_start(Position::try_from(hit.target_start as usize + 1).unwrap())
        .set_reference_sequence_id(
            target_idx
                .get(hit.target_name.as_str())
                .unwrap()
                .0,
        )
        .build();

    // 与 htslib 的 set_mapq 逐字节一致：255 表示缺失（noodles 用 None 表达，
    // encoder 写出 0xFF）；0..=254 直接透传
    *bam_record.mapping_quality_mut() = MappingQuality::new(hit.mapq);

    let data = bam_record.data_mut();

    if let Some(cs) = hit.cs.as_ref() {
        data.insert(Tag::new(b'c', b's'), Value::String(cs.as_bytes().into()));
    }

    if let Some(md) = hit.md.as_ref() {
        data.insert(Tag::new(b'm', b'd'), Value::String(md.as_bytes().into()));
    }

    if let Some(np_) = query_record.np {
        data.insert(Tag::new(b'n', b'p'), Value::UInt16(np_ as u16));
    }

    if let Some(ch_) = query_record.ch {
        data.insert(Tag::new(b'c', b'h'), Value::UInt32(ch_ as u32));
    }

    if let Some(rq_) = query_record.rq {
        data.insert(Tag::new(b'r', b'q'), Value::Float(rq_));
    }

    if let Some(be_) = &query_record.be {
        data.insert(
            Tag::new(b'b', b'e'),
            Value::Array(Array::UInt32(be_.clone())),
        );
    }

    if let Some(dw_) = &query_record.dw {
        let dw_ = if is_rev {
            dw_.iter().copied().rev().collect::<Vec<_>>()
        } else {
            dw_.clone()
        };
        data.insert(Tag::new(b'd', b'w'), Value::Array(Array::UInt8(dw_)));
    }

    if let Some(ar_) = &query_record.ar {
        let ar_ = if is_rev {
            ar_.iter().copied().rev().collect::<Vec<_>>()
        } else {
            ar_.clone()
        };
        data.insert(Tag::new(b'a', b'r'), Value::Array(Array::UInt8(ar_)));
    }

    if let Some(cr_) = &query_record.cr {
        let cr_ = if is_rev {
            cr_.iter().copied().rev().collect::<Vec<_>>()
        } else {
            cr_.clone()
        };
        data.insert(Tag::new(b'c', b'r'), Value::Array(Array::UInt8(cr_)));
    }

    if let Some(nn_) = &query_record.nn {
        let nn_ = if is_rev {
            nn_.iter().copied().rev().collect::<Vec<_>>()
        } else {
            nn_.clone()
        };
        data.insert(Tag::new(b'n', b'n'), Value::Array(Array::UInt8(nn_)));
    }

    if let Some(wd_) = &query_record.wd {
        let wd_ = if is_rev {
            wd_.iter().copied().rev().collect::<Vec<_>>()
        } else {
            wd_.clone()
        };
        data.insert(Tag::new(b'w', b'd'), Value::Array(Array::UInt8(wd_)));
    }

    if let Some(sd_) = &query_record.sd {
        let sd_ = if is_rev {
            sd_.iter().copied().rev().collect::<Vec<_>>()
        } else {
            sd_.clone()
        };
        data.insert(Tag::new(b's', b'd'), Value::Array(Array::UInt8(sd_)));
    }

    if let Some(sp_) = &query_record.sp {
        let sp_ = if is_rev {
            sp_.iter().copied().rev().collect::<Vec<_>>()
        } else {
            sp_.clone()
        };
        data.insert(Tag::new(b's', b'p'), Value::Array(Array::UInt8(sp_)));
    }

    bam_record
}
