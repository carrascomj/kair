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
//! use kair::ModelLP;
//! use std::str::FromStr;
//!
//! let file_str = std::fs::read_to_string("examples/EcoliCore.xml").unwrap();
//! let model = ModelLP::from_str(&file_str).unwrap();
//! println!(
//!     "Model has {:?} constraints and {:?} variables",
//!     &model.constraints.len(),
//!     &model.variables.len()
//! );
//! for (name, val) in model.optimize().unwrap().iter() {
//!     println!("{} = {}", name, val)
//! }
//! ```
//!
//! # Additional links
//!
//! * [Github repository](https://github.com/carrascomj/kair), open to PRs!
//! * [rust_sbml](https://docs.rs/rust_sbml/0.3.0/rust_sbml/): SBML parser in rust.
//! * [cobrapy](https://github.com/opencobra/cobrapy/): fully featured COBRA package written in Python.
mod flux_analysis;

pub use flux_analysis::fba;

use lp_modeler::dsl::*;
use lp_modeler::solvers::CbcSolver;
use rust_sbml::{Model, Parameter, Reaction, Specie, SpeciesReference};
use uuid::Uuid;

use std::collections::HashMap;
use std::ops::AddAssign;
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
    /// Metabolites from the SBML document
    pub metabolites: HashMap<String, Specie>,
    /// Reactions from the SBML document
    pub reactions: HashMap<String, Reaction>,
    /// Parsed from reactions, variables of LP problem
    pub variables: HashMap<String, LpExpression>,
    /// Parameters from the SBML document
    pub config: HashMap<String, Parameter>,
    /// Reaction id to be used as the objective in the LP problem
    pub objective: String,
    obj_expr: Option<LpExpression>,
    stoichiometry: HashMap<String, Vec<LpExpression>>,
    /// Parsed from the stoichiometry matrix
    pub constraints: Vec<LpConstraint>,
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

impl Problem for ModelLP {
    fn add_objective_expression(&mut self, expr: &LpExpression) {
        let (_, simpl_expr) = split_constant_and_expr(&simplify(expr));
        self.obj_expr = Some(simpl_expr);
    }

    fn add_constraints(&mut self, expr: &LpConstraint) {
        self.constraints.push(expr.clone());
    }
}

impl AddAssign<LpConstraint> for ModelLP {
    fn add_assign(&mut self, _rhs: LpConstraint) {
        self.add_constraints(&_rhs);
    }
}

impl AddAssign<LpExpression> for ModelLP {
    fn add_assign(&mut self, _rhs: LpExpression) {
        self.add_objective_expression(&_rhs);
    }
}

impl ModelLP {
    /// Read and call the LP builder.
    ///
    /// # Example
    /// ```
    /// use kair::ModelLP;
    /// use std::str::FromStr;
    ///
    /// let file_str = std::fs::read_to_string("examples/EcoliCore.xml").unwrap();
    /// ModelLP::from_str(&file_str).unwrap();
    /// ```
    pub fn new(input_sbml: Model) -> Self {
        let mut model = Self::from(input_sbml);
        model.populate_model();
        model
    }
    /// Covenience method to solve Flux Balance Analysis (FBA) with the CbcSolver.
    ///
    /// FBA: [https://pubmed.ncbi.nlm.nih.gov/20212490/](https://pubmed.ncbi.nlm.nih.gov/20212490/)
    ///
    /// # Example
    /// ```
    /// use kair::ModelLP;
    /// use std::str::FromStr;
    ///
    /// let file_str = std::fs::read_to_string("examples/EcoliCore.xml").unwrap();
    /// let model = ModelLP::from_str(&file_str).unwrap();
    /// println!("{:?}", model.optimize().unwrap())
    /// ```
    pub fn optimize(&self) -> Result<HashMap<String, f32>, Box<dyn std::error::Error>> {
        let solver = CbcSolver::new();

        fba(&self.into(), solver)
    }
    fn reac_expr(&self, met: &SpeciesReference, reac: &str, com: f32) -> Option<LpExpression> {
        Some(
            match met.stoichiometry {
                Some(val) => com * (val as f32),
                None => com * 1f32,
            } * self.variables[reac].to_owned(),
        )
    }
    /// Build LP problem as an FBA formulation.
    fn populate_model(&mut self) {
        let mut stoichiometry = HashMap::<String, Vec<LpExpression>>::new();
        // Build a constraint (stoichiometry) table metabolites x reactions.
        for (reac_id, reaction) in self.reactions.iter() {
            reaction.list_of_reactants.0.iter().for_each(|sref| {
                let cons = &mut stoichiometry
                    .entry(sref.species.to_owned())
                    .or_insert_with(Vec::new);
                cons.push(self.reac_expr(sref, reac_id, -1.).unwrap());
            });
            reaction.list_of_products.0.iter().for_each(|sref| {
                let cons = &mut stoichiometry
                    .entry(sref.species.to_owned())
                    .or_insert_with(Vec::new);
                cons.push(self.reac_expr(sref, reac_id, 1.).unwrap());
            });
        }
        // Then, add each metabolite column as a constraint.
        for (_, cons) in stoichiometry.iter() {
            *self += cons.sum().ge(0.);
            *self += cons.sum().le(0.);
        }
        self.stoichiometry = stoichiometry;
        // Add the problem as the objective defined in the SBML document
        let obj = self.variables[&self.objective].clone();
        *self += obj;
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
        let reactions_lp: &HashMap<String, LpExpression> = &reactions
            .iter()
            .map(|(id, reac)| {
                (
                    id.to_owned(),
                    LpExpression::ConsCont(
                        LpContinuous::new(reac.name_var().as_ref())
                            .lower_bound(reac.lb(&config))
                            .upper_bound(reac.ub(&config)),
                    ),
                )
            })
            .collect();
        let objective = model.objectives[0].to_owned();

        ModelLP {
            metabolites,
            reactions,
            variables: reactions_lp.to_owned(),
            config,
            objective,
            obj_expr: None,
            stoichiometry: HashMap::<_, _>::new(),
            constraints: Vec::<_>::new(),
        }
    }
}

impl From<&ModelLP> for LpProblem {
    fn from(model: &ModelLP) -> LpProblem {
        LpProblem {
            // TODO: from sbml
            name: "COBRA model",
            unique_name: format!("COBRA model_{}", Uuid::new_v4()),
            objective_type: LpObjective::Maximize,
            obj_expr: model.obj_expr.clone(),
            constraints: model.constraints.clone(),
        }
    }
}
