//! COBRA methods that take an LpProblem and a Solver
use crate::model::ModelLP;
use good_lp::{
    solvers::ObjectiveDirection, solvers::StaticSolver, ProblemVariables, Solution, Solver,
    SolverModel,
};

use std::collections::HashMap;
use std::sync::mpsc::channel;
use std::thread;

// generic type for the Errors implemented by different solver interfaces
type SolverError<S> = <<S as Solver>::Model as good_lp::SolverModel>::Error;

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
/// let mut model = ModelLP::from_str(&contents).unwrap();
/// println!("{:?}", fba(&mut model, default_solver).unwrap())
/// ```
pub fn fba<S: Solver>(
    model: &mut ModelLP,
    solver: S,
) -> Result<HashMap<String, f64>, SolverError<S>> {
    _fva_step(model, solver, ObjectiveDirection::Maximisation)
}

fn _fva_step<S: Solver>(
    model: &mut ModelLP,
    solver: S,
    direction: ObjectiveDirection,
) -> Result<HashMap<String, f64>, SolverError<S>> {
    let mut problem = ProblemVariables::new();
    model.populate_model(&mut problem);
    let mut problem = problem
        .optimise(direction, model.get_objective())
        .using(solver);
    model.add_constraints::<S>(&mut problem);
    let solution = problem.solve()?;
    Ok(model
        .variables
        .iter()
        .map(|(id, var)| (id.clone(), solution.value(*var)))
        .collect())
}

/// Perform [Flux Variability Analysis](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2963619//).
///
/// 1. FBA for the default objective and fix its flux to that solution.
/// 2. For each reaction:
///     1. Minimise a FBA with the reaction as objective.
///     2. Maximize a FBA with reaction as objective.
/// 3. Report solution.
///
/// The returned `HashMap<String, (f64, f64)` contains the reaction id as key
/// and a tuple of the lower possible flux and the upper possible flux, respectively.
///
/// # Example
/// ```
/// use kair::{ModelLP, flux_analysis::fva};
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
/// let mut model = ModelLP::from_str(&contents).unwrap();
/// let reactions = &model.  reactions.iter().map(|(k, _v)| k.clone()).collect::<Vec<String>>();
/// println!("Reaction  LowerFlux  UpperFlux\n{:?}", fva(
///     &mut model,
///     default_solver,
///     reactions,
/// ).unwrap())
/// ```
pub fn fva<S>(
    model: &mut ModelLP,
    solver: S,
    reactions: &[String],
) -> Result<HashMap<String, (f64, f64)>, Box<dyn std::error::Error>>
where
    S: StaticSolver + Clone + Send + Sync,
{
    let original_solution = fba(model, solver.clone())?;
    let fix_to = original_solution[&model.objective];
    let objective = model.get_objective_reaction()?;
    objective.lb = fix_to;
    objective.ub = fix_to;
    let cpus = num_cpus::get();
    let (tx, rx) = channel();
    let reacs_per_job = reactions.len() / num_cpus::get();

    for i in 0..cpus {
        let mut model = model.clone();
        let tx = tx.clone();
        let solver = solver.clone();
        let (lower, mut upper) = (i * reacs_per_job, reacs_per_job * (i + 1));
        if (cpus - 1) == i {
            upper = reactions.len()
        }
        let reactions = reactions[lower..upper].to_vec();
        thread::spawn(move || {
            for reaction in reactions {
                model.objective = reaction.clone();
                let upper_value = match fba(&mut model, solver.clone()) {
                    Ok(sol) => sol[&model.objective],
                    _ => std::f64::NAN,
                };
                let lower_value =
                    match _fva_step(&mut model, solver.clone(), ObjectiveDirection::Minimisation) {
                        Ok(sol) => sol[&model.objective],
                        _ => std::f64::NAN,
                    };
                tx.send((reaction.clone(), (lower_value, upper_value)))
                    .expect("Could not send data!");
            }
        });
    }
    let mut result = HashMap::new();
    for _ in 0..(reactions.len()) {
        let (reac_id, bounds) = rx.recv()?;
        result.insert(reac_id, bounds);
    }
    Ok(result)
}
