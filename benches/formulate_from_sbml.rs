use criterion::{criterion_group, criterion_main, Criterion};

extern crate good_lp;
extern crate kair;

use good_lp::{coin_cbc, constraint, Expression, ProblemVariables, SolverModel};
use kair::{flux_analysis::fba, flux_analysis::fva, ModelLp};
use std::str::FromStr;

fn read_ecoli() {
    let file_str = include_str!("../tests/EcoliCore.xml");
    ModelLp::from_str(file_str).unwrap();
}

fn read_recon1() {
    let file_str = include_str!("../tests/RECON1.xml");
    let mut model = ModelLp::from_str(file_str).unwrap();
    let mut problem = ProblemVariables::new();
    model.populate_model(&mut problem);
    let mut problem = problem.maximise(model.get_objective()).using(coin_cbc);
    for (_, cons) in model.stoichiometry.iter() {
        problem.add_constraint(constraint::eq(cons.iter().sum::<Expression>(), 0f32));
    }
}

fn optimize_ecoli() {
    let file_str = include_str!("../tests/EcoliCore.xml");
    let mut model = ModelLp::from_str(file_str).unwrap();
    fba(&mut model, coin_cbc).unwrap();
}

fn fva_ecoli_small() {
    let file_str = include_str!("../tests/EcoliCore.xml");
    let mut model = ModelLp::from_str(file_str).unwrap();
    let reactions: Vec<String> = model.reactions.iter().map(|(k, _v)| k.clone()).collect();
    fva(&mut model, coin_cbc, &reactions).unwrap();
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

fn fva_benchmark(c: &mut Criterion) {
    c.bench_function("FVA of E. coli core", |b| b.iter(fva_ecoli_small));
}

criterion_group!(
    benches,
    create_lp_benchmark,
    optimize_benchmark,
    create_big_lp_benchmark,
    fva_benchmark
);

criterion_main!(benches);
