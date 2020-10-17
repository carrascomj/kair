extern crate kair;
extern crate lp_modeler;

use kair::ModelLP;
use std::str::FromStr;

fn main() {
    let file_str = std::fs::read_to_string("examples/EcoliCore.xml").unwrap();
    let model = ModelLP::from_str(&file_str).unwrap();
    // println!("{:?}", model.metabolites_lp);
    println!(
        "Model has {:?} constraints",
        &model.constraints.len()
    );
    println!("Model has {:?} variables", &model.variables.len());
    for (name, val) in model.optimize().unwrap().iter() {
        println!("{} = {}", name, val)
    }
}
