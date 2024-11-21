use std::str::FromStr;

use clap::{self, Args, Parser};
use mm2::cli;
#[derive(Debug, Parser, Clone)]
#[command(version, about, long_about=None)]
pub struct Cli {
    #[arg(long = "threads")]
    pub threads: Option<usize>,

    #[arg(long="preset", default_value_t=String::from_str("map-ont").unwrap())]
    pub preset: String,

    #[command(flatten)]
    pub io_args: IoArgs,

    #[command(flatten)]
    pub index_args: cli::IndexArgs,

    #[command(flatten)]
    pub map_args: cli::MapArgs,

    #[command(flatten)]
    pub align_args: cli::AlignArgs,

    #[command(flatten)]
    pub oup_args: cli::OupArgs,
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
    #[arg(
        long = "sbr",
        help = "query file paths, 
    if multiple query_filepath are provided, 
    their query_names will be rewriten to ___0, ___1, and so on, 
    based on the order of the filenames"
    )]
    pub sbr: String,

    #[arg(long = "smc", help=".bam/.fasta/.fa/.fna/.fq/.fastq, if fasta/fastq provided, please set tn-delim and ch-idx")]
    pub smc: String,

    #[arg(short = 'p', help = "output a file named ${p}.bam")]
    pub prefix: String,

    #[arg(long="tn-delim", help="tn-delim and ch-idx is used for extract channel from fastq/fasta smc file")]
    pub target_name_delim: Option<String>,

    #[arg(long="ch-idx")]
    pub channel_idx: Option<usize>
}
