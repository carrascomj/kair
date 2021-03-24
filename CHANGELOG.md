Next release
------------

0.4.0
-----
* Replaced `lp-modeler` with `good_lp` ([#5](https://github.com/carrascomj/kair/pull/5))
* Implement [Flux Variability Analysis](https://github.com/carrascomj/kair/pull/6).
* Fix bug were doctests were racing for the file handle and sometimes failed because of that.

0.3.0
-----
* Refactored flux analysis methods to separated module.
* Flux analysis methods accept a `LpProblem` and a struct implementing
[`SolverTrain`](https://github.com/jcavat/rust-lp-modeler/blob/master/src/solvers/mod.rs)
so that the user can decide the type of solver.
* `ModelLP::optimizes(&self)` is kept as a convenience function.

0.2.0
-----
* Change in API: call `model.constraints.len()` and `model.variables.len()` instead
of accessing the now removed underlying `LpProblem`.
* Benchmarks: within noise range from previous baseline.
* Model implements the trait
[`Problem`](https://github.com/jcavat/rust-lp-modeler/blob/master/src/dsl/problem.rs#L26).
* Model implements `Into<LpProblem>`.
