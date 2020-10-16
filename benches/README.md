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
Optimize E. coli core   time:   [34.284 ms 34.515 ms 34.772 ms]
Found 12 outliers among 100 measurements (12.00%)
  9 (9.00%) high mild
  3 (3.00%) high severe

Populate E. coli core   time:   [15.368 ms 15.391 ms 15.416 ms]
Found 5 outliers among 100 measurements (5.00%)
  4 (4.00%) high mild
  1 (1.00%) high severe
```

## Plots and statistics
Criterion generates an informative HTML report after running `cargo bench` at
`target/criterion/report`.
