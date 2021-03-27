//! Structs for the formulation of the LP problem from the SBML model
use custom_error::custom_error;
use good_lp::{constraint, variable, Expression, ProblemVariables, Solver, SolverModel, Variable};
use rust_sbml::{Model, Parameter, Reaction, Species, SpeciesReference};

use std::collections::HashMap;
use std::str::FromStr;

custom_error! {
    /// Error for inconsistencies on the SBML document
    pub SBMLError
    /// When a reaction uses an unknown parameter
    InconsistentModel{/// parameter name
        param: String} = "reaction points to {param} but it does not exist in model.parameters",
    /// When a parameter.value is accessed but None
    EmptyParameter{/// parameter name
        param: String} = "the parameter {param} exists but it holds no value",
    /// When the model.objective is not in model.reactions
    InconsistentObjective{/// objective name
        obj: String} = "model.objective points to {obj}, which could not be found in the model."
}

/// LP problem as a Flux Balance Analysis formulation.
///
/// See: [What is flux balance analysis?, Orth et al., 2010](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3108565/)
///
/// Being $f(\overrightarrow{z})$ a function to optimize (historically, biomass or ATPase), $S$ and
/// stoichimetry matrix and $\overrightarrow{v}$ the flux vector representing the reactions in the reconstruction:
///
/// $$ \text{Max.} f(\overrightarrow{z}) \newline \text{subject to}\medspace S\overrightarrow{v} = 0 \newline \text{where}\medspace lb_j \le v_j \le ub_j $$
#[derive(Clone)]
pub struct ModelLp {
    /// Id from SBML document
    pub id: String,
    /// Name from SBML document
    pub name: String,
    /// Metabolites from the SBML document
    pub metabolites: HashMap<String, Species>,
    /// Reactions from the SBML document
    pub reactions: HashMap<String, ReactionLp>,
    /// Parsed from reactions, variables of LP problem
    pub variables: HashMap<String, Variable>,
    /// Parameters from the SBML document
    pub config: HashMap<String, Parameter>,
    /// Reaction id to be used as the objective in the LP problem
    pub objective: String,
    /// Parsed stoichiometry matrix
    pub stoichiometry: HashMap<String, Vec<Expression>>,
}

/// Reaction struct translated from a SBML Reaction for ease of use.
#[derive(Clone)]
pub struct ReactionLp {
    /// lower bound of the reaction
    pub lb: f64,
    /// upper bound of the reaction
    pub ub: f64,
    id: String,
    reactants: Vec<SpeciesReference>,
    products: Vec<SpeciesReference>,
}

impl ReactionLp {
    fn from_reaction(
        reaction: Reaction,
        parameters: &HashMap<String, Parameter>,
    ) -> Result<ReactionLp, SBMLError> {
        Ok(ReactionLp {
            id: format!(
                "{}_{}",
                reaction.id,
                match reaction.compartment.as_ref() {
                    Some(s) => s.to_owned(),
                    _ => String::from(""),
                }
            ),
            lb: match reaction.lower_bound.as_ref() {
                // a parameter in reaction is guaranteed to be on the list of parameters of SBML
                Some(s) => parameters
                    .get(s)
                    .ok_or(SBMLError::InconsistentModel {
                        param: s.to_owned(),
                    })?
                    .value
                    .ok_or(SBMLError::EmptyParameter {
                        param: s.to_owned(),
                    })?,
                _ => match parameters.get("cobra_default_lb") {
                    Some(param) => param.value.ok_or(SBMLError::EmptyParameter {
                        param: String::from("cobra_default_lb"),
                    })?,
                    _ => -1000.,
                },
            },
            ub: match reaction.upper_bound.as_ref() {
                // a parameter in reaction is guaranteed to be on the list of parameters of SBML
                Some(s) => parameters
                    .get(s)
                    .ok_or(SBMLError::InconsistentModel {
                        param: s.to_owned(),
                    })?
                    .value
                    .ok_or(SBMLError::EmptyParameter {
                        param: s.to_owned(),
                    })?,
                _ => match parameters.get("cobra_default_ub") {
                    Some(param) => param.value.ok_or(SBMLError::EmptyParameter {
                        param: String::from("cobra_default_ub"),
                    })?,
                    _ => 1000.,
                },
            },
            reactants: reaction.list_of_reactants.species_references,
            products: reaction.list_of_products.species_references,
        })
    }
}

impl ModelLp {
    /// Read and call the LP builder.
    ///
    /// # Example
    /// ```
    /// use kair::ModelLp;
    /// use std::str::FromStr;
    /// # use std::{fs::File, io::{BufReader, prelude::*}};
    ///
    /// # let file = std::fs::File::open("examples/EcoliCore.xml").unwrap();
    /// # let mut buf_reader = BufReader::new(file);
    /// # let mut contents = String::new();
    /// # buf_reader.read_to_string(&mut contents).unwrap();
    /// // contents is a &str containing a SBML document
    /// ModelLp::from_str(&contents).unwrap();
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
    pub fn populate_model(&mut self, problem: &mut ProblemVariables) {
        self.add_vars(problem);
        let mut stoichiometry = HashMap::<String, Vec<Expression>>::new();
        // Build a constraint (stoichiometry) table metabolites x reactions.
        for (reac_id, reaction) in self.reactions.iter() {
            reaction.reactants.iter().for_each(|sref| {
                let cons = &mut stoichiometry
                    .entry(sref.species.to_owned())
                    .or_insert_with(Vec::new);
                cons.push(self.reac_expr(sref, reac_id, -1.))
            });
            reaction.products.iter().for_each(|sref| {
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
    /// Get objective variable given the objective identifier
    pub fn get_objective_reaction(&mut self) -> Result<&mut ReactionLp, SBMLError> {
        self.reactions
            .get_mut(&self.objective)
            .ok_or(SBMLError::InconsistentObjective {
                obj: self.objective.to_owned(),
            })
    }
    /// Add the constraints to th problem
    pub fn add_constraints<S: Solver>(&self, model: &mut S::Model) {
        for (_, cons) in self.stoichiometry.iter() {
            model.add_constraint(constraint::eq(cons.iter().sum::<Expression>(), 0.));
        }
    }
    fn add_vars(&mut self, problem: &mut ProblemVariables) {
        self.variables = self
            .reactions
            .iter()
            .map(|(id, reac)| {
                (
                    id.to_owned(),
                    problem.add(variable().min(reac.lb).max(reac.ub)),
                )
            })
            .collect();
    }
}

impl FromStr for ModelLp {
    type Err = Box<dyn std::error::Error>;

    fn from_str(input_sbml: &str) -> Result<Self, Box<dyn std::error::Error>> {
        Ok(Self::new(Model::parse(input_sbml)?))
    }
}

impl From<Model> for ModelLp {
    fn from(mut model: Model) -> ModelLp {
        let metabolites = model.species;
        let config = model.parameters;
        let mut reactions = HashMap::new();
        let reac_ids: Vec<String> = model.reactions.iter().map(|(k, _)| k.to_owned()).collect();
        for key in reac_ids.iter() {
            reactions.insert(
                key.to_owned(),
                // key comes from model.reactions.keys, which is iterated just once
                ReactionLp::from_reaction(model.reactions.remove(key).unwrap(), &config).unwrap(),
            );
        }
        let objective = model.objectives.unwrap()[0].to_owned();
        let id = match model.id {
            Some(s) => s,
            _ => "".to_string(),
        };
        let name = match model.name {
            Some(s) => s,
            _ => "".to_string(),
        };

        ModelLp {
            id,
            name,
            metabolites,
            reactions,
            variables: HashMap::<_, _>::new(),
            config,
            objective,
            stoichiometry: HashMap::<_, _>::new(),
        }
    }
}
