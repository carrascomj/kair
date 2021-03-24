use criterion::{criterion_group, criterion_main, Criterion};

extern crate good_lp;
extern crate kair;

use good_lp::{coin_cbc, constraint, Expression, ProblemVariables, SolverModel};
use kair::{flux_analysis::fba, ModelLP};
use std::str::FromStr;

fn read_ecoli() {
    let file_str = include_str!("../tests/EcoliCore.xml");
    ModelLP::from_str(file_str).unwrap();
}

fn read_recon1() {
    let file_str = include_str!("../tests/RECON1.xml");
    let mut model = ModelLP::from_str(file_str).unwrap();
    let mut problem = ProblemVariables::new();
    model.populate_model(&mut problem);
    let mut problem = problem.maximise(model.get_objective()).using(coin_cbc);
    for (_, cons) in model.stoichiometry.iter() {
        problem.add_constraint(constraint::eq(cons.iter().sum::<Expression>(), 0f32));
    }
}

fn optimize_ecoli() {
    let file_str = include_str!("../tests/EcoliCore.xml");
    let mut model = ModelLP::from_str(file_str).unwrap();
    fba(&mut model, coin_cbc).unwrap();
}

fn create_big_lp_benchmark(c: &mut Criterion) {
    c.bench_function("Populate Recon1", |b| b.iter(read_recon1));
}

fn create_lp_benchmark(c: &mut Criterion) {
    c.bench_function("Populate E. coli core", |b| b.iter(read_ecoli));
}

fn optimize_benchmark(c: &mut Criterion) {
    c.bench_function("Optimize E. coli core", |b| b.iter(optimize_ecoli));
}

criterion_group!(
    benches,
    create_lp_benchmark,
    optimize_benchmark,
    create_big_lp_benchmark,
);

criterion_main!(benches);
