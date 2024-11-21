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

#[derive(Debug, Args, Clone, Copy)]
pub struct IndexArgs {
    #[arg(long, help = "minimizer kmer")]
    kmer: Option<usize>,

    #[arg(long, help = "minimizer window size")]
    wins: Option<usize>,
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

    #[arg(short = 'p', help = "output a file named ${p}.bam")]
    pub prefix: String,

    #[arg(
        long = "tn-delim",
        help = "tn-delim and ch-idx is used for extract channel from fastq/fasta smc file"
    )]
    pub target_name_delim: Option<String>,

    #[arg(long = "ch-idx")]
    pub channel_idx: Option<usize>,
}

#[derive(Debug, Args, Clone, Default)]
pub struct OupArgs {
    #[arg(long="oupIyT", default_value_t=-1.0, help="remove the record from the result bam file when the identity < identity_threshold")]
    pub oup_identity_threshold: f32,

    #[arg(long="oupCovT", default_value_t=-1.0, help="remove the record from the result bam file when the coverage < coverage_threshold")]
    pub oup_coverage_threshold: f32,
}
