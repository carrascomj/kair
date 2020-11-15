extern crate kair;
extern crate lp_modeler;

use kair::ModelLP;
use std::str::FromStr;

fn main() {
    let file_str = std::fs::read_to_string("examples/EcoliCore.xml").unwrap();
    let model = ModelLP::from_str(&file_str).unwrap();
    assert_eq!(&model.constraints.len(), &144);
}
