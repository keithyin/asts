use std::{collections::HashMap, fs, path, thread};

use asts::{align_worker, cli::Cli, subreads_and_smc_generator};
use clap::Parser;
use gskits::{
    fastx_reader::fastx2bam::{fasta2bam, fastq2bam},
    samtools::{samtools_bai, sort_by_coordinates, sort_by_tag},
};
use mm2::bam_writer::BamOupArgs;
use rust_htslib::bam::Read;

use time;
use tracing_subscriber;

fn build_target_to_idx(smc_bam: &str) -> HashMap<String, (usize, usize)> {
    let mut reader = rust_htslib::bam::Reader::from_path(smc_bam).unwrap();
    reader.set_threads(4).unwrap();

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
        let delim = args.io_args.target_name_delim.as_ref().unwrap();
        let channel_idx = args.io_args.channel_idx.unwrap();
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

    let o_path = thread::scope(|s| {
        let args = &args;
        let sorted_sbr = &sorted_sbr;
        let sorted_smc = &sorted_smc;
        let target2idx = &target2idx;

        let (sbr_and_smc_sender, sbr_and_smc_recv) = crossbeam::channel::bounded(1000);
        s.spawn(move || {
            subreads_and_smc_generator(sorted_sbr, sorted_smc, sbr_and_smc_sender);
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
            &BamOupArgs {
                iy_threshold: args.oup_args.oup_identity_threshold,
                ec_threshold: args.oup_args.oup_coverage_threshold,
                no_sencondary: true,
                no_supplementry: true,
            },
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
