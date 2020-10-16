extern crate kair;

use kair::{Fbc, ModelLP};
use std::str::FromStr;

#[test]
fn read_ecoli() {
    let file_str = std::fs::read_to_string("tests/EcoliCore.xml").unwrap();
    ModelLP::from_str(&file_str).unwrap();
}

#[test]
fn verify_bound() {
    let file_str = std::fs::read_to_string("tests/EcoliCore.xml").unwrap();
    let model = ModelLP::from_str(&file_str).unwrap();
    assert_eq!(model.reactions["R_ATPM"].lb(&model.config), 8.39);
}

#[test]
fn verify_neg_bound() {
    let file_str = std::fs::read_to_string("tests/EcoliCore.xml").unwrap();
    let model = ModelLP::from_str(&file_str).unwrap();
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
    let file_str = std::fs::read_to_string("tests/EcoliCore.xml").unwrap();
    let model = ModelLP::from_str(&file_str).unwrap();
    assert_eq!(
        (model.optimize().unwrap()["R_BIOMASS_Ecoli_core_w_GAM_"] * 10000.).round(),
        8739.
    )
}
