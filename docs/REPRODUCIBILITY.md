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
