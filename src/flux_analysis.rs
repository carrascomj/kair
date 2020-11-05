//! COBRA methods that take an LpProblem and a Solver
use lp_modeler::solvers::SolverTrait;

use std::collections::HashMap;

/// Optimize the model according to Flux Balance Analysis (FBA).
/// FBA: [https://pubmed.ncbi.nlm.nih.gov/20212490/](https://pubmed.ncbi.nlm.nih.gov/20212490/)
///
/// # Example
/// ```
/// use kair::{ModelLP, fba};
/// use std::{str::FromStr, convert::Into};
/// use lp_modeler::solvers::CbcSolver;
///
/// # use std::{fs::File, io::{BufReader, prelude::*}};
///
/// # let file = std::fs::File::open("examples/EcoliCore.xml").unwrap();
/// # let mut buf_reader = BufReader::new(file);
/// # let mut contents = String::new();
/// # buf_reader.read_to_string(&mut contents).unwrap();
/// // contents is a &str containing a SBML document
/// let model = &ModelLP::from_str(&contents).unwrap();
/// println!("{:?}", fba(&model.into(), CbcSolver::new()).unwrap())
/// ```
pub fn fba<T: SolverTrait>(
    problem: &T::P,
    solver: T,
) -> Result<HashMap<String, f32>, Box<dyn std::error::Error>> {
    let solution = solver.run(problem)?;
    Ok(solution.results)
}
