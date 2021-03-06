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

<a href="https://www.codecogs.com/eqnedit.php?latex=\begin{array}{l@{\quad}&space;r&space;c&space;r}&space;\max&space;&&space;f(\overrightarrow&space;z)&space;\\&space;\mathrm{s.t.}&space;&&space;S\overrightarrow{v}&space;=&space;0&space;\\&space;\mathrm{where}&&space;lb_j&space;\le&space;v_j&space;\le&space;ub_j&space;\end{array}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\begin{array}{l@{\quad}&space;r&space;c&space;r}&space;\max&space;&&space;f(\overrightarrow&space;z)&space;\\&space;\mathrm{s.t.}&space;&&space;S\overrightarrow{v}&space;=&space;0&space;\\&space;\mathrm{where}&&space;lb_j&space;\le&space;v_j&space;\le&space;ub_j&space;\end{array}" title="\begin{array}{l@{\quad} r c r} \max & f(\overrightarrow z) \\ \mathrm{s.t.} & S\overrightarrow{v} = 0 \\ \mathrm{where}& lb_j \le v_j \le ub_j \end{array}" /></a>

The FBA problem can then be optimized thanks to [lp_modeler](https://github.com/jcavat/rust-lp-modeler).

See [What is flux balance analysis?, Orth et al., 2010](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3108565/)
for a brief description of FBA.

## Installation
Add **kair** it to your Cargo.toml:
```toml
[dependencies]
kair = "0.5.0"
```

In addition, add [`good_lp`](https://github.com/rust-or/good_lp) with the solver of choice, for instance `coin_cbc` (default):
```toml
[dependencies]
good_lp = { version="1.1.0", default_features=true }
```

Make sure you have installed the [Cbc solver](https://github.com/coin-or/Cbc#binaries)
(other solvers do not require installation).
```shell
# Debian
sudo apt install coinor-cbc
# Arch
sudo pacman -S coin-or
# Mac OS
brew tap coin-or-tools/coinor && brew install coin-or-tools/coinor/cbc
```

## Example
Some `use` statements to get started.
```rust
use kair::{ModelLP, fba, flux_analysis::fva};
use good_lp::default_solver;
use std::str::FromStr;
```
First, read the SBML document, we will be using the [e_coli_core model](http://bigg.ucsd.edu/models/e_coli_core).
```rust
let file_str = std::fs::read_to_string("examples/EcoliCore.xml").unwrap();
let model = ModelLP::from_str(&file_str).unwrap();
```
Now, we can optimize it and print the solution, which is just a
[HashMap](https://doc.rust-lang.org/std/collections/struct.HashMap.html) of
pairs _variable name_ -> _solution value_.
```rust
for (name, val) in fba(&mut model, default_solver).unwrap().iter() {
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

Flux variability analysis is also implemented:
```rust
let reactions: Vec<String> = model.reactions.iter().map(|(k, _v)| k.clone()).collect();
for (name, val) in fva(&mut model, default_solver, reactions).unwrap().iter() {
    println!("{} = {:?}", name, val)
}
```
_Output (you would need to use a bigger model to see the difference)_
```
R_ACONTa = (6.007249575350586, 6.007249575350007)
R_ACALD = (0.0, 0.0)
R_ACKr = (-0.0, -0.0)
R_ICDHyr = (6.007249575351851, 6.007249575350007)
R_CO2t = (-22.80983331020489, -22.809833310205118)
R_RPI = (-2.2815030940668573, -2.2815030940674283)
R_ADK1 = (-0.0, -0.0000000000003395200787181807)
R_PGK = (-16.0235261431673, -16.02352614316787)
R_SUCCt3 = (0.0, -0.0000000000004168517383125921)
R_EX_pyr_e = (0.0, 0.0)
```
