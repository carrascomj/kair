extern crate kair;
extern crate good_lp;

use kair::{ModelLP, flux_analysis::fba};
use good_lp::coin_cbc;
use std::str::FromStr;

fn main() {
    let file_str = std::fs::read_to_string("examples/EcoliCore.xml").unwrap();
    let model = ModelLP::from_str(&file_str).unwrap();
    // println!("{:?}", model.metabolites_lp);
    for (name, val) in fba(model, coin_cbc).unwrap().iter() {
        println!("{} = {}", name, val)
    }
}
