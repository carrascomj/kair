extern crate kair;
extern crate lp_modeler;
extern crate log4rs;

use kair::ModelLP;
use std::str::FromStr;


fn main() {
    let file_str = std::fs::read_to_string("examples/EcoliCore.xml").unwrap();
    let model = ModelLP::from_str(&file_str).unwrap();
    log4rs::init_file("log4rs.yml", Default::default()).unwrap();
    model.optimize().unwrap();
}
