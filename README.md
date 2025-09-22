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

# Execute the targets pipeline (for main paper analysis)
Sys.setenv(TAR_PROJECT = "main")
targets::tar_make()

# Save the autometric log (to avoid overriding it when running the pipeline
# again) and pipeline metadata
fs::file_move("autometric.log", "output/autometric_main_targets.log")
targets::tar_meta() |> readr::write_csv("output/meta_main_targets.csv")

# Generate figures for the manuscript (will be saved into output/ folder)
source("make_manuscript_figures.R")

# Prepare datasets for additional example datasets
source("prepare_datasets.R")

# Execute the targets pipeline (for analysis of full vs reduced genomics dataset)
Sys.setenv(TAR_PROJECT = "li_genomics")
targets::tar_make()

# Save the autometric log and pipeline metadata
fs::file_move("autometric.log", "output/autometric_li_genomics_targets.log")
targets::tar_meta() |> readr::write_csv("output/meta_li_genomics_targets.csv")

# Execute the targets pipeline (for analysis of EATRIS-Plus dataset)
Sys.setenv(TAR_PROJECT = "eatris_plus")
targets::tar_make()

# Save the autometric log and pipeline metadata
fs::file_move("autometric.log", "output/autometric_eatris_plus_targets.log")
targets::tar_meta() |> readr::write_csv("output/meta_eatris_plus_targets.csv")

quarto::quarto_render("supplementary_material_figures.qmd")
```

If you are not familiar with the `targets` package, you can learn about it in the [targets user manual](https://books.ropensci.org/targets/).
