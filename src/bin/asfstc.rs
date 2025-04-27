use std::{
    collections::HashMap,
    fs,
    io::{BufWriter, Write},
    path,
    str::FromStr,
    sync::{Arc, Mutex},
    thread,
};

use asts::{
    reporter::Reporter, sbr_and_cs_to_cs::align_sbr_and_fake_cs_to_cs_worker,
    sbr_and_ref_to_cs::MsaResult, subreads_and_smc_generator,
};
use crossbeam::channel::Receiver;
use mm2::gskits::{
    fastx_reader::fastx2bam::{fasta2bam, fastq2bam},
    pbar::{self, DEFAULT_INTERVAL},
    samtools::sort_by_tag,
};
use asts::params::{InputFilterParams, OupParams, AlignParams, MapParams};
use rust_htslib::bam::Read;

use time;
use tracing_subscriber;

// use std::str::FromStr;

use clap::{self, Args, Parser};
#[derive(Debug, Parser, Clone)]
#[command(
    version,
    about,
    long_about = "align subreads to cs, then output msa result"
)]
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
    #[command(flatten)]
    pub align_args: AlignArgs,

    #[command(flatten)]
    pub oup_args: OupArgs,
}

#[derive(Debug, Args, Clone)]
pub struct IoArgs {
    #[arg(short = 'q', help = "subreads.bam")]
    pub sbr: String,

    #[arg(
        short = 't',
        help = "smc.bam/.fasta/.fa/.fna/.fq/.fastq, if fasta/fastq provided, please set tn-delim and ch-idx"
    )]
    pub smc: String,

    #[arg(short = 'p', help = "output a file named ${p}.asrtc.txt")]
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
        help = "1:3,5,7:9 means [[1, 3], [5, 5], [7, 9]]. target np_range. only valid for target that contains np field"
    )]
    pub np_range: Option<String>,

    #[arg(
        long = "rq-range",
        default_value_t=String::from_str("0.0:0.999").unwrap(),
        help = "0.0:0.999 means 0.0<=rq<=0.999. target rq_range. only valid for target that contains rq field"
    )]
    pub rq_range: String,
}

impl IoArgs {
    pub fn to_input_filter_params(&self) -> InputFilterParams {
        let mut param = InputFilterParams::new();
        param = if let Some(ref np_range_str) = self.np_range {
            param.set_np_range(np_range_str)
        } else {
            param
        };

        // param = if let Some(ref rq_range_str) = self.rq_range {
        //     param.set_rq_range(rq_range_str)
        // } else {
        //     param
        // };

        param = param.set_rq_range(&self.rq_range);

        param
    }
}

#[derive(Debug, Args, Clone, Default)]
pub struct AlignArgs {
    #[arg(
        short = 'm',
        default_value_t = 2,
        help = "matching_score>=0"
    )]
    matching_score: i32,

    #[arg(
        short = 'M',
        default_value_t = 5,
        help = "mismatch_penalty >=0"
    )]
    mismatch_penalty: i32,

    #[arg(short = 'o', default_value_t=String::from_str("2,24").unwrap() ,help = "gap_open_penalty >=0")]
    gap_open_penalty: String,

    #[arg(short = 'e', default_value_t=String::from_str("1,0").unwrap(), help = "gap_extension_penalty >=0")]
    gap_extension_penalty: String,
}

impl AlignArgs {
    pub fn to_align_params(&self) -> AlignParams {
        let mut param = AlignParams::new();
        param = param
            .set_m_score(self.matching_score)
            .set_mm_score(self.mismatch_penalty)
            .set_gap_open_penalty(self.gap_open_penalty.clone())
            .set_gap_extension_penalty(self.gap_extension_penalty.clone());

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
            .set_oup_coverage_threshold(self.oup_coverage_threshold)
            .set_pass_through_tags(None);
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

    let time_offset =
        time::UtcOffset::current_local_offset().unwrap_or_else(|_| time::UtcOffset::UTC);
    let timer = tracing_subscriber::fmt::time::OffsetTime::new(time_offset, time_fmt);

    let mut tmp_files = vec![];
    let args = Cli::parse();
    let o_path = format!("{}.asrtc.txt", args.io_args.prefix);
    let log_path = format!("{}.asrtc.log", args.io_args.prefix);
    let log_file = std::fs::File::create(&log_path).unwrap();
    let (non_blocking, _guard) = tracing_appender::non_blocking(log_file);
    tracing_subscriber::fmt::fmt()
        .with_timer(timer)
        .with_ansi(false)
        .with_writer(non_blocking)
        .init();

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
            tracing::error!("not a valid smc file format. valid file format : .bam/.fq/.fa/.fasta/.fna, but got:{}", args.io_args.smc);
            panic!("exit. read log file for more information. {}", log_path);
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

    let map_params = MapParams::default();
    let align_params = args.align_args.to_align_params();
    let reporter = Arc::new(Mutex::new(Reporter::default()));

    thread::scope(|s| {
        let tot_threads = args.threads.unwrap_or(num_cpus::get());
        assert!(tot_threads >= 10, "at least 10 threads");

        let args = &args;
        let sorted_sbr = &sorted_sbr;
        let sorted_smc = &sorted_smc;
        let target2idx = &target2idx;
        let oup_params = &oup_params;
        let input_filter_params = &input_filter_params;
        let map_params = &map_params;
        let align_params = &align_params;

        let (sbr_and_smc_sender, sbr_and_smc_recv) = crossbeam::channel::bounded(1000);
        let reporter_ = reporter.clone();
        s.spawn(move || {
            subreads_and_smc_generator(
                sorted_sbr,
                sorted_smc,
                input_filter_params,
                oup_params,
                sbr_and_smc_sender,
                reporter_,
            );
        });

        let align_threads = args.threads.unwrap_or(num_cpus::get()) - 4;
        let (align_res_sender, align_res_recv) = crossbeam::channel::bounded(1000);
        for idx in 0..align_threads {
            let sbr_and_smc_recv_ = sbr_and_smc_recv.clone();
            let align_res_sender_ = align_res_sender.clone();
            let reporter_ = reporter.clone();

            thread::Builder::new()
                .name(format!("align_sbr_to_smc_worker-{}", idx))
                .spawn_scoped(s, move || {
                    align_sbr_and_fake_cs_to_cs_worker(
                        sbr_and_smc_recv_,
                        align_res_sender_,
                        target2idx,
                        map_params,
                        align_params,
                        oup_params,
                        reporter_,
                    )
                })
                .unwrap();
        }
        drop(sbr_and_smc_recv);
        drop(align_res_sender);

        write_msa_result(align_res_recv, &o_path);
    });

    tracing::info!(
        "\n--------Reporter-----------\n{}\n---------------------------------",
        reporter.lock().unwrap()
    );

    for tmp_file in tmp_files {
        if path::Path::new(&tmp_file).exists() {
            fs::remove_file(&tmp_file).unwrap();
            tracing::info!("removed tmp file {}", tmp_file);
        }
    }
}

fn write_msa_result(recv: Receiver<MsaResult>, o_path: &str) {
    let mut buf_writer = BufWriter::new(fs::File::create(o_path).unwrap());

    let pb = pbar::get_spin_pb(
        format!("asctr: writing msa result to {}", o_path),
        DEFAULT_INTERVAL,
    );

    for msa in recv {
        let json_str = serde_json::to_string(&msa).unwrap();
        pb.inc(1);
        writeln!(&mut buf_writer, "{json_str}").unwrap();
    }
    pb.finish();
}
