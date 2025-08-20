use std::sync::{Arc, Mutex};

use asts::params::{InputFilterParams, OupParams};
use asts::{reporter::Reporter, subreads_and_smc_generator};
use criterion::{black_box, criterion_group, criterion_main, Criterion};
use mm2::gskits::samtools::sort_by_tag;

fn subreads_and_smc_generator_benchmark(c: &mut Criterion) {
    let sbr_bam = "/data/ccs_data/speed-test/13k-5h/4cc/20241220_Sync_Y0701_01_H01_Run0001_called.adapter.filtered.bam";
    let smc_bam = "/data/ccs_data/speed-test/13k-5h/4cc/20241220_Sync_Y0701_01_H01_Run0001_called.adapter.filtered.smc_all_reads.bam";
    let input_filter_param = InputFilterParams::new();

    let sorted_sbr = sort_by_tag(sbr_bam, "ch", None);
    let sorted_smc = sort_by_tag(smc_bam, "ch", None);
    c.bench_function("subreads_and_smc_generator_benchmark", |b| {
        b.iter(|| {
            let (sender, recv) = crossbeam::channel::unbounded();
            let recv = black_box(recv);
            let reporter = Arc::new(Mutex::new(Reporter::default()));

            subreads_and_smc_generator(
                black_box(&sorted_sbr),
                black_box(&sorted_smc),
                black_box(&input_filter_param),
                black_box(&OupParams::default()),
                black_box(sender),
                reporter,
                None,
            );
        })
    });
}

criterion_group!(benches, subreads_and_smc_generator_benchmark);
criterion_main!(benches);
