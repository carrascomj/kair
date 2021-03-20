//! COBRA methods that take an LpProblem and a Solver
use crate::ModelLP;
use good_lp::{constraint, solvers::Solver, Expression, ProblemVariables, Solution, SolverModel};

use std::collections::HashMap;

/// Optimize the model according to Flux Balance Analysis (FBA).
/// FBA: [https://pubmed.ncbi.nlm.nih.gov/20212490/](https://pubmed.ncbi.nlm.nih.gov/20212490/)
///
/// # Example
/// ```
/// use kair::{ModelLP, fba};
/// use std::{str::FromStr, convert::Into};
/// use good_lp::default_solver;
///
/// # use std::{fs::File, io::{BufReader, prelude::*}};
///
/// # let file = std::fs::File::open("examples/EcoliCore.xml").unwrap();
/// # let mut buf_reader = BufReader::new(file);
/// # let mut contents = String::new();
/// # buf_reader.read_to_string(&mut contents).unwrap();
/// // contents is a &str containing a SBML document
/// let model = ModelLP::from_str(&contents).unwrap();
/// println!("{:?}", fba(model, default_solver).unwrap())
/// ```
pub fn fba<S: Solver>(
    mut model: ModelLP,
    solver: S,
) -> Result<HashMap<String, f64>, Box<dyn std::error::Error>> {
    let mut problem = ProblemVariables::new();
    model.populate_model(&mut problem);
    let mut problem = problem.maximise(model.get_objective()).using(solver);
    for (_, cons) in model.stoichiometry.iter() {
        problem.add_constraint(constraint::eq(cons.iter().sum::<Expression>(), 0f32));
    }
    let solution = problem.solve().unwrap();
    Ok(model
        .variables
        .iter()
        .map(|(id, var)| (id.clone(), solution.value(*var)))
        .collect())
}
