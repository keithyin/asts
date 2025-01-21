use asts::subreads_and_smc_generator;
use criterion::{black_box, criterion_group, criterion_main, Criterion};
use gskits::samtools::sort_by_tag;
use mm2::params::InputFilterParams;

fn subreads_and_smc_generator_benchmark(c: &mut Criterion) {
    let sbr_bam = "";
    let smc_bam = "";
    let input_filter_param = InputFilterParams::new();

    let sorted_sbr = sort_by_tag(sbr_bam, "ch", None);
    let sorted_smc = sort_by_tag(smc_bam, "ch", None);
    c.bench_function("subreads_and_smc_generator_benchmark", |b| {
        b.iter(|| {
            let (sender, _) = crossbeam::channel::unbounded();
            subreads_and_smc_generator(
                black_box(&sorted_sbr),
                black_box(&sorted_smc),
                black_box(&input_filter_param),
                black_box(sender),
            );
        })
    });
}

criterion_group!(benches, subreads_and_smc_generator_benchmark);
criterion_main!(benches);
