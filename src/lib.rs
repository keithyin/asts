use std::collections::HashMap;

use crossbeam::channel::Sender;
use gskits::ds::ReadInfo;
use minimap2::{Aligner, PresetSet};
use mm2::params::InputFilterParams;
use rust_htslib::bam::{ext::BamRecordExtensions, Read};

// cCsSiIf int8, uint8, int16, uint16, int32, uint32, float
type BamRecord = rust_htslib::bam::record::Record;

pub struct SubreadsAndSmc {
    pub smc: ReadInfo,
    pub subreads: Vec<ReadInfo>,
}

pub struct SingleChannelAlignRes {
    pub records: Vec<BamRecord>,
}

impl SubreadsAndSmc {
    pub fn new(smc_record: &rust_htslib::bam::Record) -> Self {
        Self {
            smc: ReadInfo::from_bam_record(smc_record, None),
            subreads: vec![],
        }
    }

    pub fn add_subread(&mut self, record: &rust_htslib::bam::Record) {
        self.subreads.push(ReadInfo::from_bam_record(record, None));
    }
}

pub fn subreads_and_smc_generator(
    sorted_sbr_bam: &str,
    sorted_smc_bam: &str,
    input_filter_params: &InputFilterParams,
    sender: crossbeam::channel::Sender<SubreadsAndSmc>,
) {
    let mut smc_bam_reader = rust_htslib::bam::Reader::from_path(sorted_smc_bam).unwrap();
    smc_bam_reader.set_threads(5).unwrap();
    let mut subreads_bam_reader = rust_htslib::bam::Reader::from_path(sorted_sbr_bam).unwrap();
    subreads_bam_reader.set_threads(5).unwrap();

    let mut smc_records = smc_bam_reader.records();
    let mut subreads_records = subreads_bam_reader.records();

    let mut smc_record_opt = smc_records.next();
    let mut sbr_record_opt = subreads_records.next();

    loop {
        if smc_record_opt.is_none() {
            break;
        }

        let smc_record = smc_record_opt.unwrap().unwrap();

        if !input_filter_params.valid(&smc_record) {
            smc_record_opt = smc_records.next();
            continue;
        }

        let mut subreads_and_smc = SubreadsAndSmc::new(&smc_record);
        let smc = gskits::gsbam::bam_record_ext::BamRecordExt::new(&smc_record);

        let mut found = false;
        loop {
            if sbr_record_opt.is_none() {
                break;
            }

            let sbr = sbr_record_opt.as_ref().unwrap().as_ref().unwrap();
            let sbr = gskits::gsbam::bam_record_ext::BamRecordExt::new(sbr);
            if sbr.get_ch().unwrap() != smc.get_ch().unwrap() {
                if found {
                    break;
                }
                sbr_record_opt = subreads_records.next();
                continue;
            }

            found = true;
            let sbr = sbr_record_opt.as_ref().unwrap().as_ref().unwrap();
            subreads_and_smc.add_subread(sbr);

            sbr_record_opt = subreads_records.next();
        }

        sender.send(subreads_and_smc).unwrap();
        if sbr_record_opt.is_none() {
            break;
        }
        smc_record_opt = smc_records.next();
    }
}

pub fn align_worker(
    recv: crossbeam::channel::Receiver<SubreadsAndSmc>,
    sender: Sender<mm2::AlignResult>,
    target_idx: &HashMap<String, (usize, usize)>,
) {
    for subreads_and_smc in recv {
        // tracing::info!("sbr_cnt:{}-{}", subreads_and_smc.smc.name, subreads_and_smc.subreads.len());
        let align_infos = align(&subreads_and_smc, target_idx);
        let align_res = mm2::AlignResult {
            records: align_infos
                .into_iter()
                .filter(|rec| rec.is_some())
                .map(|v| v.unwrap())
                .collect::<Vec<_>>(),
        };
        sender.send(align_res).unwrap();
    }
}

pub fn align(
    subreads_and_smc: &SubreadsAndSmc,
    target_idx: &HashMap<String, (usize, usize)>,
) -> Vec<Option<BamRecord>> {
    let aligner = build_asts_aligner(subreads_and_smc.smc.seq.len() < 200);

    let aligner = aligner
        .with_seq_and_id(
            subreads_and_smc.smc.seq.as_bytes(),
            subreads_and_smc.smc.name.as_bytes(),
        )
        // .with_seq(subreads_and_smc.smc.seq.as_bytes())
        .unwrap();

    let mut bam_records = vec![];

    for subread in subreads_and_smc.subreads.iter() {
        let hits = aligner
            .map(
                subread.seq.as_bytes(),
                true,
                true,
                None,
                Some(&[0x4000000, 0x40000000]),
                Some(subread.name.as_bytes())
            )
            .unwrap();
        let mut bam_record = None;

        for hit in hits {
            // no supp is needed !!
            if hit.is_primary && !hit.is_supplementary {
                bam_record = Some(mm2::build_bam_record_from_mapping(
                    &hit, subread, target_idx,
                ));
            }
        }
        bam_records.push(bam_record);
    }
    bam_records
}

/// https://github.com/PacificBiosciences/actc/blob/main/src/PancakeAligner.cpp#L128
/// https://github.com/lh3/minimap2/blob/master/minimap.h
fn build_asts_aligner(short_insert: bool) -> Aligner<PresetSet> {
    let mut aligner = Aligner::builder()
        .map_ont()
        .with_cigar() // cigar_str has bug in minimap2="0.1.20+minimap2.2.28"
        .with_sam_out()
        .with_index_threads(1);

    aligner.idxopt.k = 7;
    aligner.idxopt.w = 5;

    // this is all map-ont default
    // aligner.mapopt.zdrop = 400;
    // aligner.mapopt.zdrop_inv = 200;
    // aligner.mapopt.bw = 500;
    // aligner.mapopt.end_bonus = 0;

    aligner.mapopt.best_n = 1;
    aligner.mapopt.q_occ_frac = 0.0;

    if short_insert {
        aligner.idxopt.k = 4;
        aligner.idxopt.w = 1;
        aligner.mapopt.min_cnt = 2;
        aligner.mapopt.min_dp_max = 10; // min dp score
        aligner.mapopt.min_chain_score = 10; // this is important for short insert
        aligner.mapopt.min_ksw_len = 0;
    }

    aligner

}


pub fn draw_aligned_seq(record: &BamRecord, ref_seq: &[u8], r_start: Option<usize>, r_end: Option<usize>) -> (String, String) {

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
            unsafe {(*query_seq.get_unchecked(qpos_ as usize)) as char}
        } else {
            '-'
        };

        let r_char = if let Some(rpos_) = rpos{
            unsafe {(*ref_seq.get_unchecked(rpos_ as usize)) as char}

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
    use gskits::{fastx_reader::fasta_reader::FastaFileReader, gsbam::bam_record_ext::{BamReader, BamRecordExt}};
    use mm2::build_bam_record_from_mapping;

    use super::*;

    #[test]
    fn test_align() {
        let subreads_and_smc = SubreadsAndSmc {
            smc: ReadInfo::new_fa_record(
                "hello".to_string(),
                "CACTCTTTAAAAGGGGGGTTGAGG".to_string(),
            ),
            subreads: vec![
                ReadInfo::new_fa_record("".to_string(), "CACTCTTTAAAAGGGGGGGGTTGAGG".to_string()),
                ReadInfo::new_fa_record("".to_string(), "CACTATTTAAAAGGGGGGTTGAGG".to_string()),
                ReadInfo::new_fa_record("".to_string(), "CACTCTTTAAAAGGGGGGTTGAGA".to_string()),
                ReadInfo::new_fa_record("".to_string(), "CACTCTTTAAAAGGGGGGAAGAGG".to_string()),
            ],
        };

        let mut target2idx = HashMap::new();
        target2idx.insert("hello".to_string(), (0, subreads_and_smc.smc.seq.len()));
        let records = align(&subreads_and_smc, &target2idx);
        for record in records {
            println!("{:?}", record.unwrap().cigar());
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

        let mut aligner = minimap2::Aligner::builder()
            .map_ont()
            .with_cigar() // cigar_str has bug in minimap2="0.1.20+minimap2.2.28"
            .with_sam_out()
            .with_index_threads(1);

        aligner.idxopt.k = 7;
        aligner.idxopt.w = 5;

        let aligner = aligner
            .with_seq_and_id(ref_fa.seq.as_bytes(), ref_fa.name.as_bytes())
            // .with_seq(subreads_and_smc.smc.seq.as_bytes())
            .unwrap();

        let hits = aligner
            .map(
                query_fa.seq.as_bytes(),
                true,
                true,
                None,
                Some(&[0x4000000, 0x40000000]),
                Some(query_fa.name.as_bytes())
            )
            .unwrap();

        for hit in hits {
            println!("query2ref {:?}", hit)
        }

        let mut aligner = minimap2::Aligner::builder()
            .map_ont()
            .with_cigar() // cigar_str has bug in minimap2="0.1.20+minimap2.2.28"
            .with_sam_out()
            .with_index_threads(1);

        aligner.idxopt.k = 7;
        aligner.idxopt.w = 5;

        let aligner = aligner
            .with_seq_and_id(query_fa.seq.as_bytes(), query_fa.name.as_bytes())
            // .with_seq(subreads_and_smc.smc.seq.as_bytes())
            .unwrap();

        let hits = aligner
            .map(
                ref_fa.seq.as_bytes(),
                true,
                true,
                None,
                Some(&[0x4000000, 0x40000000]),
                Some(ref_fa.name.as_bytes())
            )
            .unwrap();

        for hit in hits {
            println!("ref2query {:?}", hit)
        }
    }

    #[test]
    fn test_build_aligner() {
        let aligner = build_asts_aligner(true);
        println!("{:?}", aligner.idxopt);
        println!("{:?}", aligner.mapopt);

        let aligner = build_asts_aligner(false);
        println!("{:?}", aligner.idxopt);
        println!("{:?}", aligner.mapopt);
    }

    #[test]
    fn test_short_insert() {
        let aligner = build_asts_aligner(true);
        println!("{:?}", aligner.idxopt);
        println!("{:?}", aligner.mapopt);
        println!("\n\n\n");

        let ref_seq = "CCCCTTTTAAAAGGGGGAGAGA";
        let q_seq = "CCCCTTTTAAAAGGGGGAGAGA";
        let aligner = aligner.with_seq(ref_seq.as_bytes()).unwrap();

        for hit in aligner
            .map(q_seq.as_bytes(), false, false, None, None, Some(b"query"))
            .unwrap()
        {
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
        let aligner = build_asts_aligner(false);
        let aligner = aligner.with_seq_and_id(icing_reads.as_bytes(), b"icing").unwrap();
        let mut target_idx = HashMap::new();
        target_idx.insert("icing".to_string(), (0, icing_reads.len()));
        for hit in aligner.map(sbr.as_bytes(), false, false, None, None, Some(b"sbr")).unwrap() {
            if hit.is_primary {
                let query_record = ReadInfo::new_fa_record("name".to_string(), sbr.to_string());
                let record = build_bam_record_from_mapping(&hit, &query_record, &target_idx);

                let (ref_aligned, query_aligned) = draw_aligned_seq(&record, icing_reads.as_bytes(), Some(590), Some(620));
                println!("\n{}\n{}", ref_aligned, query_aligned)

            }

        }


    }

}
