extern crate kai;
extern crate lp_modeler;

use kai::ModelLP;
use std::str::FromStr;

fn main() {
    let file_str = std::fs::read_to_string("examples/EcoliCore.xml").unwrap();
    let model = ModelLP::from_str(&file_str).unwrap();
    // println!("{:?}", model.metabolites_lp);
    println!(
        "Model has {:?} constraints",
        &model.problem.constraints.len()
    );
    println!("Model has {:?} variables", &model.problem.variables().len());
    for (name, val) in model.optimize().unwrap().iter() {
        println!("{} = {}", name, val)
    }
}
