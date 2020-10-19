extern crate kair;

use kair::{Fbc, ModelLP};
use std::str::FromStr;

const EXAMPLE: &'static str = include_str!("../tests/EcoliCore.xml");

#[test]
fn read_ecoli() {
    ModelLP::from_str(&EXAMPLE).unwrap();
}

#[test]
fn verify_bound() {
    let model = ModelLP::from_str(&EXAMPLE).unwrap();
    assert_eq!(model.reactions["R_ATPM"].lb(&model.config), 8.39);
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
    assert_eq!(model.reactions["R_EX_glc__D_e"].lb(&model.config), -10.);
}

#[test]
fn optimize_ecoli() {
    let model = ModelLP::from_str(&EXAMPLE).unwrap();
    assert_eq!(
        (model.optimize().unwrap()["R_BIOMASS_Ecoli_core_w_GAM_"] * 10000.).round(),
        8739.
    )
}
