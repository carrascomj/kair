0.2.0
------------
* Change in API: call `model.constraints.len()` and `model.variables.len()` instead
of accessing the now removed underlying `LpProblem`.
* Benchmarks: within noise range from previous baseline.
* Model implements the trait
[`Problem`](https://github.com/jcavat/rust-lp-modeler/blob/master/src/dsl/problem.rs#L26).
* Model implements `Into<LpProblem>`.
