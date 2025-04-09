# moiraine_publication_code

This repository contains the code used to generate the results presented in Angelin-Bonnet et al, 2025 "moiraine: an R package to construct reproducible pipelines for the application and comparison of multi-omics integration methods".

## Repository content

Key files:

- `renv.lock`: record of the packages and their version used in the repository.

- `_targets.R`: code used for the analysis, in the form of a `targets` pipeline (see the [targets user manual](https://books.ropensci.org/targets/) for more information)

- `helper_functions.R`: script, collection of custom R functions used in the analysis (mostly to generate the plots and tables for the manuscript).

## Reproduce the analysis

In order to reproduce the analysis, run the following commands in R (from the project's root directory):

```r
# Install necessary packages with correct version
renv::restore()

# Execute the targets pipeline
targets::tar_make()
```

If you are not familiar with the `targets` package, you can learn about it in the [targets user manual](https://books.ropensci.org/targets/).
