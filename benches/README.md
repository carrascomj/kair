# Benchmarks
* **Optimize E. coli core**: reads the model from SBML, generates the problem and
optimizes the FBA formulation.
* **Populate E. coli core**: reads the model from SBML and generates the problem.

## Run the benchmarks yourself
```
cargo bench
```

## Results
```
Optimize E. coli core   time:   [22.771 ms 22.860 ms 22.962 ms]
Found 9 outliers among 100 measurements (9.00%)
  7 (7.00%) high mild
  2 (2.00%) high severe

Populate E. coli core   time:   [5.2549 ms 5.2796 ms 5.3092 ms]
Found 7 outliers among 100 measurements (7.00%)
  2 (2.00%) high mild
  5 (5.00%) high severe
```

## Plots and statistics
Criterion generates an informative HTML report after running `cargo bench` at
`target/criterion/report`.
