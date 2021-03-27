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
//! * [Flux Variablity Analysis](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2963619/).
//!
//! # Examples
//!
//! Read and optimize [e_coli_core model](http://bigg.ucsd.edu/models/e_coli_core)
//! ```
//! use kair::{ModelLp, flux_analysis::fba};
//! use std::{str::FromStr, fs::File, io::{BufReader, prelude::*}};
//! use good_lp::default_solver;
//!
//! let file = std::fs::File::open("examples/EcoliCore.xml").unwrap();
//! let mut buf_reader = BufReader::new(file);
//! let mut contents = String::new();
//! buf_reader.read_to_string(&mut contents).unwrap();
//! let mut model = ModelLp::from_str(&contents).unwrap();
//! for (name, val) in fba(&mut model, default_solver).unwrap().iter() {
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
pub mod model;

pub use flux_analysis::fba;
pub use model::ModelLp;
