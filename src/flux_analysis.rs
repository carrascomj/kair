//! COBRA methods that take an LpProblem and a Solver
use lp_modeler::solvers::SolverTrait;

use std::collections::HashMap;

/// Optimize the model according to Flux Balance Analysis (FBA).
/// FBA: [https://pubmed.ncbi.nlm.nih.gov/20212490/](https://pubmed.ncbi.nlm.nih.gov/20212490/)
///
/// # Example
/// ```ignore
/// use kair::{ModelLP, fba};
/// use std::{str::FromStr, convert::Into};
/// use lp_modeler::solvers::CbcSolver;
///
/// let model = &ModelLP::from_str(&include_str!("../tests/EcoliCore.xml")).unwrap();
/// println!("{:?}", fba(&model.into(), CbcSolver::new()).unwrap())
/// ```
pub fn fba<T: SolverTrait>(
    problem: &T::P,
    solver: T,
) -> Result<HashMap<String, f32>, Box<dyn std::error::Error>> {
    let (_, solution) = solver.run(problem)?;
    Ok(solution)
}

#[cfg(test)]
mod tests {
    use super::super::*;
    use std::convert::Into;

    #[test]
    fn direct_fba() {
        let file_str = include_str!("../examples/EcoliCore.xml");
        let model = &ModelLP::from_str(&file_str).unwrap();
        println!("{:?}", fba(&model.into(), CbcSolver::new()).unwrap())
    }
}
