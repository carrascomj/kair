use criterion::{criterion_group, criterion_main, Criterion};

extern crate kair;

use good_lp::coin_cbc;
use kair::{flux_analysis::fba, ModelLP};

use std::str::FromStr;

fn read_ecoli() {
    let file_str = include_str!("../tests/EcoliCore.xml");
    ModelLP::from_str(file_str).unwrap();
}

fn optimize_ecoli() {
    let file_str = include_str!("../tests/EcoliCore.xml");
    let model = ModelLP::from_str(file_str).unwrap();
    fba(model, coin_cbc).unwrap();
}

fn create_lp_benchmark(c: &mut Criterion) {
    c.bench_function("Populate E. coli core", |b| b.iter(read_ecoli));
}

fn optimize_benchmark(c: &mut Criterion) {
    c.bench_function("Optimize E. coli core", |b| b.iter(optimize_ecoli));
}

criterion_group!(benches, optimize_benchmark, create_lp_benchmark);
criterion_main!(benches);
