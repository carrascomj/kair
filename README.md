# KAIr (COBRA Alternative In rust)

[![Crates.io](https://img.shields.io/crates/v/kair.svg)](https://crates.io/crates/kair)
[![Documentation](https://docs.rs/kair/badge.svg)](https://docs.rs/kair/)
[![Build](https://github.com/carrascomj/kair/workflows/build/badge.svg)](https://github.com/carrascomj/kair)
[![Codecov](https://codecov.io/github/carrascomj/kair/coverage.svg?branch=trunk)](https://codecov.io/gh/carrascomj/kair)

*COnstraint-Based Reconstruction and Analysis* (COBRA) methods
enable the use of knowledge-based reconstructions of the metabolism of a
particular organism to simulate its metabolic network.

**kair** provides the translation from a [SBML](http://sbml.org/Special/specifications/sbml-level-3/version-2/core/release-2/sbml-level-3-version-2-release-2-core.pdf) (using [rust_sbml](https://github.com/carrascomj/rust_sbml/)) document to the most basic
Linear Programming formulation of COBRA: Flux Balance Analysis (FBA). Being
`f(z)` a function to optimize (historically, the biomass pseudoreaction or the ATPase),
`S` and stoichimetry matrix; and `v` the flux vector representing
the reactions in the reconstruction:

<a href="https://www.codecogs.com/eqnedit.php?latex=\setlength\arraycolsep{1.5pt}&space;\begin{array}{l@{\quad}&space;r&space;c&space;r&space;c&space;r}&space;\max&space;&&space;f(\vec{z})&space;\\&space;\mathrm{s.t.}&space;&&space;S\vec(v)&space;\\&space;\mathrm{where}&space;&&space;lb_j&space;\le&space;v_j&space;\le&space;ub_j&space;\end{array}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\setlength\arraycolsep{1.5pt}&space;\begin{array}{l@{\quad}&space;r&space;c&space;r&space;c&space;r}&space;\max&space;&&space;f(\vec{z})&space;\\&space;\mathrm{s.t.}&space;&&space;S\vec(v)&space;\\&space;\mathrm{where}&space;&&space;lb_j&space;\le&space;v_j&space;\le&space;ub_j&space;\end{array}" title="\setlength\arraycolsep{1.5pt} \begin{array}{l@{\quad} r c r c r} \max & f(\vec{z}) \\ \mathrm{s.t.} & S\vec(v) \\ \mathrm{where} & lb_j \le v_j \le ub_j \end{array}" /></a>

The FBA problem can then be optimized thanks to [lp_modeler](https://github.com/jcavat/rust-lp-modeler).

See [What is flux balance analysis?, Orth et al., 2010](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3108565/)
for a brief description of FBA.

## Example
Add **kair** it to your Cargo.toml:
```toml
[dependencies]
kair = "0.1.0"
```
Some `use` statements to get started.
```rust
use kair::ModelLP;
use std::str::FromStr;
```
First, read the SBML document, we will be using the [e_coli_core model](http://bigg.ucsd.edu/models/e_coli_core).
```rust
let file_str = std::fs::read_to_string("examples/EcoliCore.xml").unwrap();
let model = ModelLP::from_str(&file_str).unwrap();
```
Having read the document, the LP problem is already formulated. We can print
some information about it:
```rust
println!(
    "Model has {:?} constraints",
    &model.problem.constraints.len()
);
println!("Model has {:?} variables", &model.problem.variables().len());
```
_Output_
```
Model has 144 constraints
Model has 95 variables
```
Finally, we can optimize it and print the solution, which is just a
[HashMap](https://doc.rust-lang.org/std/collections/struct.HashMap.html) of
pairs _variable name_ -> _solution value_.
```rust
for (name, val) in model.optimize().unwrap().iter() {
    println!("{} = {}", name, val)
}
```
_Output_
```
R_EX_co2_e_ = 22.809834
R_ATPM_ = 8.39
R_H2Ot_ = -29.175827
R_GLNS_ = 0.22346173
...
R_BIOMASS_Ecoli_core_w_GAM_ = 0.8739215
...
R_EX_pi_e_ = -3.214895
R_SUCOAS_ = -5.064376
R_PGL_ = 4.959985
R_TKT1_ = 1.4969838
```

To run this example, on the root of this repository, run
```shell
cargo run --example ecoli
```