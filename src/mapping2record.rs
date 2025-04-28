use std::{borrow::Cow, collections::HashMap, hash::Hash};

use gskits::dna::reverse_complement;

use rust_htslib::bam::record::{Aux, AuxArray};

pub use gskits;

use crate::read_info::ReadInfo;

pub type BamRecord = rust_htslib::bam::record::Record;
pub type BamWriter = rust_htslib::bam::Writer;
pub type BamReader = rust_htslib::bam::Reader;

pub fn build_bam_record_from_mapping<T>(
    hit: &minimap2::Mapping,
    query_record: &ReadInfo,
    target_idx: &HashMap<T, (usize, usize)>,
) -> BamRecord
where
    T: std::borrow::Borrow<str> + Eq + Hash,
{
    // println!("{:?}", hit);

    let mut bam_record = BamRecord::new();

    let mut seq = Cow::Borrowed(query_record.seq.as_str());
    let mut is_rev = false;

    match hit.strand {
        minimap2::Strand::Reverse => {
            is_rev = true;
            seq = Cow::Owned(String::from_utf8(reverse_complement(seq.as_bytes())).unwrap());
        }
        _ => {}
    };

    if is_rev {
        bam_record.set_reverse();
    }

    // hit.query_start, hit.query_end 是相对于原始 query 而言的(即 未 reverse 的 query 而言)
    let aln_info = hit.alignment.as_ref().unwrap();
    let cigar_str = mm2::convert_mapping_cigar_to_record_cigar(
        aln_info.cigar.as_ref().unwrap(),
        hit.query_start as usize,
        hit.query_end as usize,
        seq.len(),
        is_rev,
    );

    let qual = if let Some(ref qual_) = query_record.qual {
        if is_rev {
            qual_.iter().copied().rev().collect()
        } else {
            qual_.clone()
        }
    } else {
        vec![255; seq.len()]
    };

    bam_record.set(
        query_record.name.as_bytes(),
        Some(&cigar_str),
        seq.as_bytes(),
        &qual,
    );
    if is_rev {
        bam_record.set_reverse();
    }

    // reference start
    bam_record.set_pos(hit.target_start as i64);
    bam_record.set_mpos(-1);
    // bam_record.set_mpos(mpos);
    // bam_record.reference_end()
    bam_record.set_mapq(hit.mapq as u8);

    bam_record.set_tid(
        target_idx
            .get(hit.target_name.as_ref().unwrap().as_str())
            .unwrap()
            .0 as i32,
    );
    bam_record.set_mtid(-1);

    if hit.is_primary {
        bam_record.unset_secondary();
        bam_record.unset_supplementary();
    } else {
        bam_record.set_secondary();
        bam_record.set_supplementary();
    }

    if hit.is_supplementary {
        bam_record.unset_secondary();
        bam_record.set_supplementary();
    } else {
        bam_record.unset_supplementary();
    }
    bam_record.unset_unmapped();

    if let Some(cs) = aln_info.cs.as_ref() {
        bam_record
            .push_aux(b"cs", rust_htslib::bam::record::Aux::String(cs))
            .unwrap();
    }

    if let Some(md) = aln_info.md.as_ref() {
        bam_record
            .push_aux(b"md", rust_htslib::bam::record::Aux::String(md))
            .unwrap();
    }

    if let Some(np_) = query_record.np {
        bam_record.push_aux(b"np", Aux::U16(np_ as u16)).unwrap();
    }

    if let Some(ch_) = query_record.ch {
        bam_record.push_aux(b"ch", Aux::U32(ch_ as u32)).unwrap();
    }

    if let Some(rq_) = query_record.rq {
        bam_record.push_aux(b"rq", Aux::Float(rq_)).unwrap();
    }

    if let Some(be_) = &query_record.be {
        bam_record
            .push_aux(b"be", Aux::ArrayU32(AuxArray::from(be_)))
            .unwrap();
    }

    if let Some(dw_) = &query_record.dw {
        if is_rev {
            let dw_ = dw_.iter().copied().rev().collect::<Vec<_>>();
            bam_record
                .push_aux(b"dw", Aux::ArrayU8(AuxArray::from(&dw_)))
                .unwrap();
        } else {
            bam_record
                .push_aux(b"dw", Aux::ArrayU8(AuxArray::from(dw_)))
                .unwrap();
        }
    }

    if let Some(ar_) = &query_record.ar {
        if is_rev {
            let ar_ = ar_.iter().copied().rev().collect::<Vec<_>>();
            bam_record
                .push_aux(b"ar", Aux::ArrayU8(AuxArray::from(&ar_)))
                .unwrap();
        } else {
            bam_record
                .push_aux(b"ar", Aux::ArrayU8(AuxArray::from(ar_)))
                .unwrap();
        }
    }

    if let Some(cr_) = &query_record.cr {
        if is_rev {
            let cr_ = cr_.iter().copied().rev().collect::<Vec<_>>();
            bam_record
                .push_aux(b"cr", Aux::ArrayU8(AuxArray::from(&cr_)))
                .unwrap();
        } else {
            bam_record
                .push_aux(b"cr", Aux::ArrayU8(AuxArray::from(cr_)))
                .unwrap();
        }
    }

    if let Some(nn_) = &query_record.nn {
        if is_rev {
            let nn_ = nn_.iter().copied().rev().collect::<Vec<_>>();
            bam_record
                .push_aux(b"nn", Aux::ArrayU8(AuxArray::from(&nn_)))
                .unwrap();
        } else {
            bam_record
                .push_aux(b"nn", Aux::ArrayU8(AuxArray::from(nn_)))
                .unwrap();
        }
    }

    if let Some(wd_) = &query_record.wd {
        let mut wd_ = Cow::Borrowed(wd_);
        if is_rev {
            wd_ = Cow::Owned(wd_.iter().copied().rev().collect::<Vec<_>>());
        }
        bam_record
            .push_aux(b"wd", Aux::ArrayU8(AuxArray::from(wd_.as_ref())))
            .unwrap();
    }

    if let Some(sd_) = &query_record.sd {
        let mut sd_ = Cow::Borrowed(sd_);
        if is_rev {
            sd_ = Cow::Owned(sd_.iter().copied().rev().collect::<Vec<_>>());
        }
        bam_record
            .push_aux(b"sd", Aux::ArrayU8(AuxArray::from(sd_.as_ref())))
            .unwrap();
    }

    if let Some(sp_) = &query_record.sp {
        let mut sp_ = Cow::Borrowed(sp_);
        if is_rev {
            sp_ = Cow::Owned(sp_.iter().copied().rev().collect::<Vec<_>>());
        }
        bam_record
            .push_aux(b"sp", Aux::ArrayU8(AuxArray::from(sp_.as_ref())))
            .unwrap();
    }

    bam_record
}
