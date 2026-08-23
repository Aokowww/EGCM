# Reproducibility

## Suggested workflow (R)

1. Create an RStudio Project in the repository root.
2. Use `renv` to lock package versions:

```r
install.packages("renv")
renv::init()           # first time
renv::snapshot()       # after installing required packages
# Later, on a fresh machine:
# renv::restore()
```

3. Run simulations / analysis scripts from the repo root so relative paths resolve correctly.

## Outputs

Write all derived outputs (figures, tables, intermediate `.RData`) to `results/`.

## Public validation commands

Run the core checks from the repository root:

```bash
Rscript tests/test_cluster_gcm_core.R
Rscript tests/test_adni_stress_dgp.R
Rscript code/simulation/simulate_non_iid_bc.R \
  configs/non_iid_bc_simulation_smoke.yaml
```

The frozen 1,000-repetition Case B/C extension uses
`configs/non_iid_bc_simulation.yaml`. Its reviewed aggregate output is stored
in `results_public/simulation/non_iid_bc_summary.csv`; replicate-level and
checkpoint outputs remain under the ignored `results/` directory.

The ADNI scripts require independent ADNI approval and locally downloaded
source tables. They never fetch participant rows from the public repository.
Only the reviewed aggregate tables under `results_public/adni/` are versioned.
