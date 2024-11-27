use std::{collections::HashMap, fs, path, thread};

use asts::{align_worker, subreads_and_smc_generator};
use gskits::{
    fastx_reader::fastx2bam::{fasta2bam, fastq2bam},
    samtools::{samtools_bai, sort_by_coordinates, sort_by_tag},
};
use mm2::params::{InputFilterParams, OupParams};
use rust_htslib::bam::Read;

use time;
use tracing_subscriber;

// use std::str::FromStr;

use clap::{self, Args, Parser};
#[derive(Debug, Parser, Clone)]
#[command(version, about, long_about=None)]
pub struct Cli {
    #[arg(long = "threads")]
    pub threads: Option<usize>,

    // #[arg(long="preset", default_value_t=String::from_str("map-ont").unwrap())]
    // pub preset: String,
    #[command(flatten)]
    pub io_args: IoArgs,

    // #[command(flatten)]
    // pub index_args: cli::IndexArgs,

    // #[command(flatten)]
    // pub map_args: cli::MapArgs,

    // #[command(flatten)]
    // pub align_args: cli::AlignArgs,
    #[command(flatten)]
    pub oup_args: OupArgs,
}

// #[derive(Debug, Args, Clone, Copy)]
// pub struct IndexArgs {
//     #[arg(long, help = "minimizer kmer")]
//     kmer: Option<usize>,

//     #[arg(long, help = "minimizer window size")]
//     wins: Option<usize>,
// }

#[derive(Debug, Args, Clone)]
pub struct IoArgs {
    #[arg(short = 'q', help = "subreads.bam")]
    pub sbr: String,

    #[arg(
        short = 't',
        help = "smc.bam/.fasta/.fa/.fna/.fq/.fastq, if fasta/fastq provided, please set tn-delim and ch-idx"
    )]
    pub smc: String,

    #[arg(short = 'p', help = "output a file named ${p}.bam")]
    pub prefix: String,

    #[arg(
        long = "tn-delim",
        help = "tn-delim and ch-idx is used for extract channel from fastq/fasta smc file"
    )]
    pub target_name_delim: Option<String>,
    #[arg(long = "ch-idx")]
    pub channel_idx: Option<usize>,

    #[arg(
        long = "np-range",
        help = "1-3,5,7-9 means [[1, 3], [5, 5], [7, 9]]. smc np_range. only valid for smc input that contains np field"
    )]
    pub np_range: Option<String>,

    #[arg(
        long = "rq-range",
        help = "0.9~1.1 means 0.9<=rq<=1.1. smc rq_range. only valid for smc input that contains rq field"
    )]
    pub rq_range: Option<String>,
}

impl IoArgs {
    pub fn to_input_filter_params(&self) -> InputFilterParams {
        let mut param = InputFilterParams::new();
        param = if let Some(ref np_range_str) = self.np_range {
            param.set_np_range(np_range_str)
        } else {
            param
        };

        param = if let Some(ref rq_range_str) = self.rq_range {
            param.set_rq_range(rq_range_str)
        } else {
            param
        };

        param
    }
}

#[derive(Debug, Args, Clone, Default)]
pub struct OupArgs {
    #[arg(long="oupIyT", default_value_t=-1.0, help="remove the record from the result bam file when the identity < identity_threshold")]
    pub oup_identity_threshold: f32,

    #[arg(long="oupCovT", default_value_t=-1.0, help="remove the record from the result bam file when the coverage < coverage_threshold")]
    pub oup_coverage_threshold: f32,
}

impl OupArgs {
    pub fn to_oup_params(&self) -> OupParams {
        let mut param = OupParams::new();
        param = param
            .set_discard_secondary(true)
            .set_discard_supplementary(true)
            .set_oup_identity_threshold(self.oup_identity_threshold)
            .set_oup_coverage_threshold(self.oup_coverage_threshold);
        param
    }
}

fn build_target_to_idx(smc_bam: &str) -> HashMap<String, (usize, usize)> {
    let mut reader = rust_htslib::bam::Reader::from_path(smc_bam).unwrap();
    reader.set_threads(10).unwrap();

    let mut target2idx = HashMap::new();

    for (idx, record) in reader.records().enumerate() {
        let record = record.unwrap();

        let qname = unsafe { String::from_utf8_unchecked(record.qname().to_owned()) };

        let qlen = record.seq_len();
        target2idx.insert(qname, (idx, qlen));
    }
    target2idx
}

fn main() {
    let time_fmt = time::format_description::parse(
        "[year]-[month padding:zero]-[day padding:zero] [hour]:[minute]:[second]",
    )
    .unwrap();

    // let timer = tracing_subscriber::fmt::time::OffsetTime::new(time_offset, time_fmt);
    // let timer = tracing_subscriber::fmt::time::LocalTime::new(time_fmt);
    let time_offset =
        time::UtcOffset::current_local_offset().unwrap_or_else(|_| time::UtcOffset::UTC);
    let timer = tracing_subscriber::fmt::time::OffsetTime::new(time_offset, time_fmt);

    tracing_subscriber::fmt::fmt().with_timer(timer).init();

    let mut tmp_files = vec![];

    let args = Cli::parse();
    let smc_fname = if args.io_args.smc.ends_with(".bam") {
        args.io_args.smc.clone()
    } else {
        let delim = args
            .io_args
            .target_name_delim
            .as_ref()
            .expect(&format!("--tn-delim & --ch-idx need to be provided"));
        let channel_idx = args
            .io_args
            .channel_idx
            .expect(&format!("--tn-delim & --ch-idx need to be provided"));
        if args.io_args.smc.ends_with("fq") || args.io_args.smc.ends_with("fastq") {
            tracing::info!("fastq2bam, {} -> .bam file", args.io_args.smc);
            fastq2bam(&args.io_args.smc, delim, channel_idx)
        } else if args.io_args.smc.ends_with("fa")
            || args.io_args.smc.ends_with("fasta")
            || args.io_args.smc.ends_with("fna")
        {
            tracing::info!("fasta2bam, {} -> .bam file", args.io_args.smc);
            fasta2bam(&args.io_args.smc, delim, channel_idx)
        } else {
            panic!("not a valid smc file format");
        }
    };

    if !smc_fname.eq(&args.io_args.smc) {
        tmp_files.push(smc_fname.clone());
    }

    tracing::info!("sorting sbr.bam {}", args.io_args.sbr);
    let sorted_sbr = sort_by_tag(&args.io_args.sbr, "ch", None);

    tracing::info!("sorting smc.bam {}", smc_fname);
    let sorted_smc = sort_by_tag(&smc_fname, "ch", None);

    tmp_files.push(sorted_sbr.clone());
    tmp_files.push(sorted_smc.clone());

    tracing::info!("building target to idx map");
    let target2idx = build_target_to_idx(&sorted_smc);

    let oup_params = args.oup_args.to_oup_params();
    let input_filter_params = args.io_args.to_input_filter_params();

    let o_path = thread::scope(|s| {
        let args = &args;
        let sorted_sbr = &sorted_sbr;
        let sorted_smc = &sorted_smc;
        let target2idx = &target2idx;
        let oup_params = &oup_params;
        let input_filter_params = &input_filter_params;

        let (sbr_and_smc_sender, sbr_and_smc_recv) = crossbeam::channel::bounded(1000);
        s.spawn(move || {
            subreads_and_smc_generator(
                sorted_sbr,
                sorted_smc,
                input_filter_params,
                sbr_and_smc_sender,
            );
        });

        let threads = args.threads.unwrap_or(num_cpus::get_physical());
        let (align_res_sender, align_res_recv) = crossbeam::channel::bounded(1000);
        for _ in 0..threads {
            let sbr_and_smc_recv_ = sbr_and_smc_recv.clone();
            let align_res_sender_ = align_res_sender.clone();
            s.spawn(move || {
                align_worker(sbr_and_smc_recv_, align_res_sender_, target2idx);
            });
        }
        drop(sbr_and_smc_recv);
        drop(align_res_sender);

        let o_path = format!("{}.bam", args.io_args.prefix);
        mm2::bam_writer::write_bam_worker(
            align_res_recv,
            target2idx,
            &o_path,
            oup_params,
            "asts",
            env!("CARGO_PKG_VERSION"),
            true,
        );
        o_path
    });

    tracing::info!("sorting result bam");
    sort_by_coordinates(&o_path, None);

    tracing::info!("indexing result bam");
    samtools_bai(&o_path, true, None).unwrap();

    for tmp_file in tmp_files {
        if path::Path::new(&tmp_file).exists() {
            fs::remove_file(&tmp_file).unwrap();
            tracing::info!("removed tmp file {}", tmp_file);
        }
    }
}
