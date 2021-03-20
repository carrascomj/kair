#![deny(unsafe_code)]
#![deny(bare_trait_objects)]
#![deny(unconditional_recursion)]
#![warn(missing_docs)]
#![warn(trivial_casts)]
#![warn(trivial_numeric_casts)]
#![warn(unreachable_pub)]
#![warn(unused_qualifications)]

//! Constraint-Based Reconstruction and Analysis.
//! It uses the [rust_sbml](https://docs.rs/rust_sbml/0.3.0/rust_sbml/) to read a
//! [SBML](http://sbml.org/Special/specifications/sbml-level-3/version-2/core/release-2/sbml-level-3-version-2-release-2-core.pdf)
//! document and translates it into a LP formulation using [lp_modeler](https://jcavat.github.io/rust-lp-modeler/lp_modeler/index.html)
//!
//! # COBRA methods available
//!
//! * [Flux Balance Analysis](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3108565/).
//!
//! # Examples
//!
//! Read and optimize [e_coli_core model](http://bigg.ucsd.edu/models/e_coli_core)
//! ```
//! use kair::{ModelLP, flux_analysis::fba};
//! use std::{str::FromStr, fs::File, io::{BufReader, prelude::*}};
//! use good_lp::default_solver;
//!
//! let file = std::fs::File::open("examples/EcoliCore.xml").unwrap();
//! let mut buf_reader = BufReader::new(file);
//! let mut contents = String::new();
//! buf_reader.read_to_string(&mut contents).unwrap();
//! let model = ModelLP::from_str(&contents).unwrap();
//! for (name, val) in fba(model, default_solver).unwrap().iter() {
//!     println!("{} = {}", name, val)
//! }
//! ```
//!
//! # Additional links
//!
//! * [Github repository](https://github.com/carrascomj/kair), open to PRs!
//! * [rust_sbml](https://docs.rs/rust_sbml/0.3.0/rust_sbml/): SBML parser in rust.
//! * [cobrapy](https://github.com/opencobra/cobrapy/): fully featured COBRA package written in Python.
pub mod flux_analysis;

pub use flux_analysis::fba;

use good_lp::{
    constraint, variable, Constraint, Expression, ProblemVariables, Solver, SolverModel, Variable,
};
use rust_sbml::{Model, Parameter, Reaction, Species, SpeciesReference};

use std::collections::HashMap;
use std::str::FromStr;

/// LP problem as a Flux Balance Analysis formulation.
///
/// See: [What is flux balance analysis?, Orth et al., 2010](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3108565/)
///
/// Being $f(\overrightarrow{z})$ a function to optimize (historically, biomass or ATPase), $S$ and
/// stoichimetry matrix and $\overrightarrow{v}$ the flux vector representing the reactions in the reconstruction:
///
/// $$ \text{Max.} f(\overrightarrow{z}) \newline \text{subject to}\medspace S\overrightarrow{v} = 0 \newline \text{where}\medspace lb_j \le v_j \le ub_j $$
pub struct ModelLP {
    /// Id from SBML document
    pub id: String,
    /// Name from SBML document
    pub name: String,
    /// Metabolites from the SBML document
    pub metabolites: HashMap<String, Species>,
    /// Reactions from the SBML document
    pub reactions: HashMap<String, Reaction>,
    /// Parsed from reactions, variables of LP problem
    variables: HashMap<String, Variable>,
    /// Parameters from the SBML document
    pub config: HashMap<String, Parameter>,
    /// Reaction id to be used as the objective in the LP problem
    pub objective: String,
    stoichiometry: HashMap<String, Vec<Expression>>,
    /// Parsed from the stoichiometry matrix
    pub constraints: Vec<Constraint>,
}

/// Indicate that the struct can be translated to a flux balance variable (generally, reactions)
pub trait Fbc {
    /// Lower bound for the variable
    fn lb(&self, _parameters: &HashMap<String, Parameter>) -> f32 {
        0.
    }
    /// Upper bound for the variable
    fn ub(&self, _parameters: &HashMap<String, Parameter>) -> f32 {
        0.
    }
    /// Name of the variable to be used in the LP formulation
    fn name_var(&self) -> String {
        "".to_string()
    }
}

impl Fbc for Reaction {
    fn lb(&self, parameters: &HashMap<String, Parameter>) -> f32 {
        match self.lower_bound.as_ref() {
            // a parameter in reaction is guaranteed to be on the list of parameters of SBML
            Some(s) => parameters[s].value.unwrap() as f32,
            _ => parameters["cobra_default_lb"].value.unwrap() as f32,
        }
    }
    fn ub(&self, parameters: &HashMap<String, Parameter>) -> f32 {
        match self.upper_bound.as_ref() {
            // a parameter in reaction is guaranteed to be on the list of parameters of SBML
            Some(s) => parameters[s].value.unwrap() as f32,
            _ => parameters["cobra_default_ub"].value.unwrap() as f32,
        }
    }
    fn name_var(&self) -> String {
        format!(
            "{}_{}",
            self.id,
            match self.compartment.as_ref() {
                Some(s) => s.to_owned(),
                _ => String::from(""),
            }
        )
    }
}

impl ModelLP {
    /// Read and call the LP builder.
    ///
    /// # Example
    /// ```
    /// use kair::ModelLP;
    /// use std::str::FromStr;
    /// # use std::{fs::File, io::{BufReader, prelude::*}};
    ///
    /// # let file = std::fs::File::open("examples/EcoliCore.xml").unwrap();
    /// # let mut buf_reader = BufReader::new(file);
    /// # let mut contents = String::new();
    /// # buf_reader.read_to_string(&mut contents).unwrap();
    /// // contents is a &str containing a SBML document
    /// ModelLP::from_str(&contents).unwrap();
    /// ```
    pub fn new(input_sbml: Model) -> Self {
        Self::from(input_sbml)
    }

    fn reac_expr(&self, met: &SpeciesReference, reac: &str, com: f64) -> Expression {
        Expression::from(self.variables[reac].to_owned())
            * match met.stoichiometry {
                Some(val) => com * val,
                None => com * 1f64,
            }
    }
    /// Build LP problem as an FBA formulation.
    fn populate_model(&mut self, problem: &mut ProblemVariables) {
        self.add_vars(problem);
        let mut stoichiometry = HashMap::<String, Vec<Expression>>::new();
        // Build a constraint (stoichiometry) table metabolites x reactions.
        for (reac_id, reaction) in self.reactions.iter() {
            reaction
                .list_of_reactants
                .species_references
                .iter()
                .for_each(|sref| {
                    let cons = &mut stoichiometry
                        .entry(sref.species.to_owned())
                        .or_insert_with(Vec::new);
                    cons.push(self.reac_expr(sref, reac_id, -1.))
                });
            reaction
                .list_of_products
                .species_references
                .iter()
                .for_each(|sref| {
                    let cons = &mut stoichiometry
                        .entry(sref.species.to_owned())
                        .or_insert_with(Vec::new);
                    cons.push(self.reac_expr(sref, reac_id, 1.));
                });
        }
        self.stoichiometry = stoichiometry;
    }
    /// Get objective variable given the objective identifier
    pub fn get_objective(&self) -> Variable {
        self.variables[&self.objective]
    }
    /// Add the constraints to th problem
    pub fn add_constraints<S: Solver>(&self, model: &mut S::Model) {
        for (_, cons) in self.stoichiometry.iter() {
            model.add_constraint(constraint::eq(cons.iter().sum::<Expression>(), 0f32));
        }
    }
    fn add_vars(&mut self, problem: &mut ProblemVariables) {
        self.variables = self
            .reactions
            .iter()
            .map(|(id, reac)| {
                (
                    id.to_owned(),
                    problem.add(
                        variable()
                            .min(reac.lb(&self.config))
                            .max(reac.ub(&self.config)),
                    ),
                )
            })
            .collect();
    }
}

impl FromStr for ModelLP {
    type Err = Box<dyn std::error::Error>;

    fn from_str(input_sbml: &str) -> Result<Self, Box<dyn std::error::Error>> {
        Ok(Self::new(Model::parse(input_sbml)?))
    }
}

impl From<Model> for ModelLP {
    fn from(model: Model) -> ModelLP {
        let metabolites = model.species;
        let config = model.parameters;
        let reactions = model.reactions;
        let objective = model.objectives.unwrap()[0].to_owned();
        let id = match model.id {
            Some(s) => s,
            _ => "".to_string(),
        };
        let name = match model.name {
            Some(s) => s,
            _ => "".to_string(),
        };

        ModelLP {
            id,
            name,
            metabolites,
            reactions,
            variables: HashMap::<_, _>::new(),
            config,
            objective,
            stoichiometry: HashMap::<_, _>::new(),
            constraints: Vec::<_>::new(),
        }
    }
}
