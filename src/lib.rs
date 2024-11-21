use std::collections::HashMap;

use crossbeam::channel::Sender;
use rust_htslib::bam::Read;

pub mod cli;
use gskits::gsbam::bam_record_ext::BamRecordExt;
use gskits::dna::reverse_complement;
type BamRecord = rust_htslib::bam::record::Record;

pub fn add(left: u64, right: u64) -> u64 {
    left + right
}

pub struct Smc {
    pub name: String,
    pub seq: String,
}

impl Smc {
    pub fn from_bam_record(record: &rust_htslib::bam::Record) -> Self {
        let record = BamRecordExt::new(record);

        Self {
            name: record.get_qname(),
            seq: record.get_seq(),
        }
    }
}

pub struct Subread {
    pub name: String,
    pub seq: String,
}

impl Subread {
    pub fn from_bam_record(record: &rust_htslib::bam::Record) -> Self {
        let record = BamRecordExt::new(record);

        Self {
            name: record.get_qname(),
            seq: record.get_seq(),
        }
    }

    pub fn reverse(&self) -> Self {
        Self {
            name: self.name.clone(),
            seq: reverse_complement(&self.seq),
        }
    }
}

pub struct SubreadsAndSmc {
    pub smc: Smc,
    pub subreads: Vec<Subread>,
}

pub struct SingleChannelAlignRes {
    pub records: Vec<BamRecord>,
}

impl SubreadsAndSmc {
    pub fn new(smc_record: &rust_htslib::bam::Record) -> Self {
        Self {
            smc: Smc::from_bam_record(smc_record),
            subreads: vec![],
        }
    }

    pub fn add_subread(&mut self, record: &rust_htslib::bam::Record) {
        self.subreads.push(Subread::from_bam_record(record));
    }
}

pub fn subreads_and_smc_generator(
    sorted_sbr_bam: &str,
    sorted_smc_bam: &str,
    sender: crossbeam::channel::Sender<SubreadsAndSmc>,
) {
    let mut smc_bam_reader = rust_htslib::bam::Reader::from_path(sorted_smc_bam).unwrap();
    smc_bam_reader.set_threads(5).unwrap();
    let mut subreads_bam_reader = rust_htslib::bam::Reader::from_path(sorted_sbr_bam).unwrap();
    subreads_bam_reader.set_threads(5).unwrap();

    let mut smc_records = smc_bam_reader.records();
    let mut subreads_records = subreads_bam_reader.records();

    let mut smc_record_opt = smc_records.next();
    let mut sbr_record_opt  = subreads_records.next();

    loop {
        if smc_record_opt.is_none() {
            break;
        }

        let smc_record = smc_record_opt.unwrap().unwrap();
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
    sender: Sender<mm2::SingleQueryAlignResult>,
    target_idx: &HashMap<String, (usize, usize)>,
) {
    for subreads_and_smc in recv {
        let align_infos = align(&subreads_and_smc, target_idx);
        let align_res = mm2::SingleQueryAlignResult {
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
    let mut aligner = minimap2::Aligner::builder()
        .map_ont()
        .with_cigar() // cigar_str has bug in minimap2="0.1.20+minimap2.2.28"
        .with_sam_out()
        .with_index_threads(1);

    aligner.idxopt.k = 7;
    aligner.idxopt.w = 5;

    aligner = aligner
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
            )
            .unwrap();
        let mut bam_record = None;

        for hit in hits {
            if hit.is_primary && !hit.is_supplementary {
                let q_record = mm2::fille_reader::QueryRecord {
                    qname: subread.name.clone(),
                    sequence: subread.seq.clone(),
                    ch: None,
                    np: None
                };
                bam_record = Some(mm2::build_bam_record_from_mapping(
                    &hit, &q_record, target_idx,
                ));
            }
        }
        bam_records.push(bam_record);
    }
    bam_records
}

#[cfg(test)]
mod tests {
    use gskits::fastx_reader::fasta_reader::FastaFileReader;
    use minimap2::Aligner;

    use super::*;

    #[test]
    fn it_works() {
        let result = add(2, 2);
        assert_eq!(result, 4);
    }
    

    #[test]
    fn test_align() {
        let subreads_and_smc = SubreadsAndSmc{
            smc: Smc { name: "hello".to_string(), seq: "CACTCTTTAAAAGGGGGAGAGAGATCTCAGGATGATTCAAAGATTGAAGAAGTTAAGAAAAAACGGGGATAAGGCATGAAAGATTGCGAGAATTGTTGATGTAGAGAAGGTTGAG
AAAAAGTTTCTCGGCAAGCCTATTACCGTGTGGAAACTTATTGGAACATCCCCAAGATGTTCCCACTATTAGAGAAAAAAGTTAGAGAACATCCAGCAGTGTGGACATCTTCGAATACGATATTCCATTTGCAAAGAGATACCTCATCGACAAAGGCCTAATACCAATGGA
GGGGGAAGAAGAGCTAAAGATTCTTGCCTTCGATATAGAAACCCTCTATCACGAAGGAGAAGAGTTTGGAAAGGCCCAATTATAATGATTAGTTATGCAGATGAGAATGAAGCAAAGGTGATTACTTGGAAAAACATAGATCTTCCATACGTTGAGGTGTATCAAGCGAGA
GAGAGATGATAAAGAGATTTCTCAGGATTATCAGGGAGAGGATCCTGACATTATATGTTACTTATAATGGAGACTCATTCGACTTCCATATTTAGCGAAAAGGGCAGAAAAACTTGGGGATTAAATTAACCATTGGAAGAGATGGAAGCGAGCCCAAGATGCAGAGAATAG
GCGATATGACGGCTGTAGAAGTCAAGGAAGAATACATTTCGACTGTATCATGTAATAACAAGGACATAAATCTCCCAACATACACACTAGAGGCTGTATATGAAGCAATTTTTGGAAAGCCAAAGGAGAAGGTATACGCCGACGAGATAGCAAAAGCCTGGGAAAGTGGAG
AGAACCTTGAGAGAGTTGCCAAATACTCGATGGAAGATGCAAAGGCAACTTATGAACTCGGGAAAGAATTCCTTCCAATGGAAATTCAGCTTTCAAGATTAGTTGGACAACCTTTATGGGATGTTTCAAGGTCAAGCACAGGGAAGACCTGTAGAGTGGTTCTTACTTAGG
AAGCCTACGAAAGAAACGAAGTAGCATCTCTCTCCCCCTTTTAAAAGGG".to_string()},
            subreads: vec![
                Subread{name: "".to_string(), seq: "AAAAACCCTTTTAAAAGGGGGAGAGAGATTCAGGATGATTCAAAGATTGAAGAAGTTACCAGAAAATAACGGGATAAAGGCATGAAAGATTGCGAGAATGTTGATGTAGAGAAGGTTGA
GAAAAAGTTTCTCGGCAAGCCTATTACCGTGTGTGAAACTTTATTTGGAACTCCCCAAGATGTTCCCACTATTAGAGAAAAAGTTAGAGAACACCAAGCACAGTTTGTGGACATCTTACCCGATACATATTCCATTTGCAAGAGATACCTCATCACAAAGGCCTAATACCA
ATGAGGGCGGAAGAAGAGAAAGATTCTTGCCTTCGTAAGAAACCCTCTATCACGAGGAGAAGAGTTTGGAAAAGGCCCAATTATAATGATTAGTTTATCCAGATGATGAATGAAGCAAAGGTGATTACATTGGGAAAAACATAGACTCCATAACGTTGAGGTGTATCAAGC
GAGAGAGAGATGATAAAGAGATTTCTCAGGATTATCAGGGAGAGGATCCGACATTATAGTTACTTATAATGAGACTCATTCGACTTCCATATTTAAGCGAAAGGGCAGAAAAACTTTGGGATAAATTAACCATTGGAAGGATGGAGCGAGCCAAAGATCCAGAGAATAGGG
ATATGACGGCGTAGAAGTCAGGGAAGAAAACTCACTGTATCATGTAATAACAAGGACAATAAATCTCCCAACATACACACTAGGGCTGTATATGAAGCAATTTTTGGAAAGCAAAGGAGGAAGGTATACGCCGACGAGATAGCAAAAGCCTGGAAAGTGGAGAAATTGAGA
GAGTTGCCAAATACTCGATGGAAGATGCAAAAGGCAACTTAGAACTCGGGAAGAATTCCTTCCAATGGAAATTCAGCTTTCAGGGATTAGTGGACACCTTTTATGGGATGTTCAGGTCAAGCACAGGGACCTGTAGAGTGGTTCTTACTTAGGAAGCACGAAAGAAACGAA
GTAGCATCCTCTCCCCCTTTTAAAAGGG".to_string()},
            
                Subread{name: "".to_string(), seq: "ACCCTTTTAAAAGGGGGAGAGAGATGCACGTTCTTTCGTAGGCTTCCTAAGTAAGAACCATCACAGGTTTCCCTGTGCTTTGACCTTGAAACTCCCCATAAGGTTGTCCAACATCTTGG
CTGTAATTTTTCCATTGGAAGGAATTCTTTCC
GAGTACATAAGTTGCCTTTGCATCTTCCATCGAGTATTTGGCAACTCTCTCAAGGTTCTCCCTCCACTTTCCCAGGCTTTTCCTATCTCGTCGGCGTATACCTTCTCCTTTGGCTTTCAAAATTGCTTCATATAAAGAGCCTCTAGTGTGTATGTTGGGAGATTAGTCCTT
GTTTATTACATGATACGAAAGTCGAAATGTATTTCTTCCTGACTTCTACAGCCGTCATATCGCCTATTCTCTGCAT
CTTGGGCTCGCTTCCATCTCTCCAATGGTTAATTTAATCCCCAAGTTTTTCTGCCCTTTTCGGCTAAATAGGAAGTGAATGAGTCTCCCCCCCCCAAATTATAAGTAACATATAATGTCAGAATCCTTCTCCCTGATAATCCTGAGAAATCTCGTTTATCATCTCTTCTCG
CTTGATACAACCTCAACGTAGGAAGATCTAATGTTTTTCCAAGTAATCACCTTTTGCTTCATTTTCATCTCATACT
AATCAATTATAATTGGGCCTTCCAAACTCTTCTCCTTCCGTGATAGAGGGTTCATATCGAAGCAAGAATCTTTAGCTCTTCTTCCCCCTCCATGGTAATTTAGGGCCTTTGTCGAATGAGGTATCTCTTTGCCTCTTGCAAATGGAATATTCGTATATCAAGATGGTCCAC
ACTGCTGGATGTTCTCTAACTTTTTTCTCTAATAGTGGGAACATCTTGGGGATGTTTCCAATAAGTTTTCCACACG
TAATAGGCTTGCCGAGAAATTTTTCCAACCTTCTCTAACATCAACAATCTCGCAATCTTTCAATGCCTTTCCCCCGTTTTTTCTTAACTTTCTTCAATCTTTGAATCCCCTGAGACTCCTTCCCCCTTTTAAAGAGGGG".to_string()},
                Subread{name: "".to_string(), seq:  "CCCCCTTTAAAAGGGGGAGAGAGATCTCAGGATGATTCAAAGATTGAAGAAGTTAAGAAAAACGGGGAAGGCATGAAGATTGCGAGGAATTGTTGATGGTAGAGAAGGTTGAGAAAAGTTTCTCGGCAAGCTATTACCCGTGTAGAACTTA
TTGGAACATCCCCAAGATGTTCCCACTATTAGAAAAAAGTTAGAGAACATCCAGCAGTGTGGACATCTTCGAATACGATATTCCATTTGCAAAGAGATACCTCATCGAGACAAACCTAATACCAATGGAGGGGAAGAAGAGCTAAAGATCTGCCTCGATATAGAAACCCTCTATCACGAAGGAGAAGATTTGGAAAGGCCCAATTATAATGATTAGTTATGCAGAGAGAATGAAGCAAAGGTGATTA
CTTGGAAAACATAAGATCTTCCATACGTTGAGATGTATCAAGCGAGAGAGAGTGTGATCAAGAGATTCTCAGGATTATCAGGGAGAGGATCCTGACATTAATAATAGTTACTTATAATGGAGAACTCATTCGACTTCCCCATATTTAGCGAAAGGGGCAGAAAAACTTGGGATTAATTTAACCATTGGAAGAATGGAAGCGAGCCAGATGCAGAGAAAGCATATGACAGGCTGTAGAAGTCAAGGAA
GAATACATTTCGACTGTATCATGTAATAACAGGACTAAATCTCCACATACAACTTAAGAGGCTGTAGCATTTTTGGAAAGCCAAGGAGAGGTATACGCCGACGAAGAATAGCAACAAGCCTGGGAAGTGGAGAAACCTTGAGAGAGGTTGCCAAATACTCGATGGAAGATGCATTGGCACTTATGAACTCGGGAAAGAATTCCTCCATGGAAATTCAGCTTTCAAGATAGTTGGACAACCTTTATGG
GATGTTTCCAAGGCAAGCACAGGGAACCTGTAGAGTGGTTCTTACTTAGGAAAGCCACGAAAGAAACGAAGTAGCATGCTCTCTCCCCTTTTAAACAGGGAAACC".to_string()},
                Subread{name: "".to_string(), seq:"CCCCTTTTAAAAGGGGGAGAGAGATGCTACTTCGTTTCTTTCGTAGGCTTTCTAAGTAAAACCACTCTAACAGGTTCCACCCGTTCCCTGTGCTTGACCTTGAAACACCATAAAAAGGTTGGTCCAACTAATCTTTGAAAAGCTGAAATTCCATATGGAAGGATTCTTTCCCGAGTACATAAGTGCCCTTTGCATCTTCCACGAGTATTTGGCAACTCTCTCAAGGTTCCTCTCCACTTTCCCAGGCTTTTCTATCTCGTCGCGTAATACCTTTCTCCTTTGGCTTTTCCGAAAAAATTGCTTCAATACAGCCCTAGGTGTGTAGTTGGGAGATTTATTGTCCTTGTTATACATGATACAAGTCAGAAATGTAATTCTCCCTTGACTCTACAGCGGTCATATCGCCTATCTCTGCATCTTGGGCTCCCGCTTCCATCTCTTCCAATGGTTAATTTAAATCCCCAATTTTCTGCCCTTTTTCGCTAAAATATGGGAGTCGAATGAGTCTA".to_string()}
]
    };

    let mut target2idx = HashMap::new();
    target2idx.insert("hello".to_string(), (0, subreads_and_smc.smc.seq.len()));
    align(&subreads_and_smc, &target2idx);

    }


    #[test]
    fn test_align2() {
        let query_fa = FastaFileReader::new("../test_data/query.fa".to_string());
        let query_fa = gskits::fastx_reader::read_fastx(query_fa);
        let query_fa = &query_fa[0];

        let ref_fa = FastaFileReader::new("../test_data/ref.fa".to_string());
        let ref_fa = gskits::fastx_reader::read_fastx(ref_fa);
        let ref_fa = &ref_fa[0];



        let mut aligner = minimap2::Aligner::builder()
        .map_ont()
        .with_cigar() // cigar_str has bug in minimap2="0.1.20+minimap2.2.28"
        .with_sam_out()
        .with_index_threads(1);

        aligner.idxopt.k = 7;
        aligner.idxopt.w = 5;

        aligner = aligner
            .with_seq_and_id(
                ref_fa.sequence.as_bytes(),
                ref_fa.name.as_bytes(),
            )
            // .with_seq(subreads_and_smc.smc.seq.as_bytes())
            .unwrap();


        let hits = aligner
            .map(
                query_fa.sequence.as_bytes(),
                true,
                true,
                None,
                Some(&[0x4000000, 0x40000000]),
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

        aligner = aligner
            .with_seq_and_id(
                query_fa.sequence.as_bytes(),
                query_fa.name.as_bytes(),
            )
            // .with_seq(subreads_and_smc.smc.seq.as_bytes())
            .unwrap();


        let hits = aligner
            .map(
                ref_fa.sequence.as_bytes(),
                true,
                true,
                None,
                Some(&[0x4000000, 0x40000000]),
            )
            .unwrap();

        for hit in hits {
            println!("ref2query {:?}", hit)
        }


    }

}
