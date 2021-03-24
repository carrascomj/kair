extern crate good_lp;
extern crate kair;

use good_lp::coin_cbc;
use kair::{flux_analysis::fba, ModelLP};
use std::str::FromStr;

fn main() {
    let file_str = std::fs::read_to_string("examples/EcoliCore.xml").unwrap();
    let mut model = ModelLP::from_str(&file_str).unwrap();
    // println!("{:?}", model.metabolites_lp);
    for (name, val) in fba(&mut model, coin_cbc).unwrap().iter() {
        println!("{} = {}", name, val)
    }
}
