extern crate kair;

use good_lp::default_solver;
use kair::{
    flux_analysis::{fba, fva},
    ModelLP,
};
use std::str::FromStr;

const EXAMPLE: &str = include_str!("../tests/EcoliCore.xml");

#[test]
fn read_ecoli() {
    ModelLP::from_str(&EXAMPLE).unwrap();
}

#[test]
fn verify_bound() {
    let model = ModelLP::from_str(&EXAMPLE).unwrap();
    assert_eq!((model.reactions["R_ATPM"].lb * 100.).round() as i32, 839);
}

#[test]
fn verify_neg_bound() {
    let model = ModelLP::from_str(&EXAMPLE).unwrap();
    println!(
        "{:?}",
        &model
            .reactions
            .iter()
            .map(|(id, _)| id.to_string())
            .filter(|id| id.starts_with("R_EX"))
            .collect::<Vec::<String>>()
    );
    assert_eq!(model.reactions["R_EX_glc__D_e"].lb.round() as i32, -10);
}

#[test]
fn optimize_ecoli() {
    let mut model = ModelLP::from_str(&EXAMPLE).unwrap();
    assert_eq!(
        (fba(&mut model, default_solver).unwrap()["R_BIOMASS_Ecoli_core_w_GAM"] * 10000.).round()
            as i32,
        8739
    )
}
#[test]
fn flux_variability_analysis_looks_fine() {
    let mut model = ModelLP::from_str(&EXAMPLE).unwrap();
    let reactions: Vec<String> = model.reactions.iter().map(|(k, _v)| k.clone()).collect();
    let sol = fva(&mut model, default_solver, &reactions).unwrap();
    let total_flux: f64 = sol.values().map(|(low, up)| low + up).sum();
    println!("{:?}", sol);
    assert_eq!(
        (sol["R_BIOMASS_Ecoli_core_w_GAM"].0 * 10000.).round() as i32,
        8739
    );
    assert!(total_flux > 0f64);
}
