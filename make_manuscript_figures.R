library(targets)
library(autometric)
library(tidyverse)
library(here)
library(fs)

source(here("R/manuscript_figures_helper_functions.R"))

dir_create("output")

dpi <- 900

# Saving figures from pipeline --------------------------------------------

ggsave(
  here("output/figure_pca_rnaseq.pdf"),
  plot = tar_read(figure_pca_rnaseq),
  width = 7,
  height = 10,
  units = "in",
  dpi = dpi
)

ggsave(
  here("output/figure_samples_score_integration.pdf"),
  plot = tar_read(figure_samples_score_integration),
  width = 8,
  height = 8,
  units = "in",
  dpi = dpi
)

ggsave(
  here("output/figure_mofa_features.pdf"),
  plot = tar_read(figure_mofa_features),
  width = 10,
  height = 8,
  units = "in",
  dpi = dpi
)

ggsave(
  here("output/figure_methods_comparison.pdf"),
  plot = tar_read(figure_methods_comparison),
  width = 8.5,
  height = 10,
  units = "in",
  dpi = dpi
)

p <- moiraine::plot_running_time(
  # target_patterns = c(
  #   "sPLS-DA\n(pre-filtering)" = "^individual_splsda",
  #   "sPLS" = "^spls_", 
  #   "sO2PLS" = "^so2pls_", 
  #   "MOFA" = "^mofa_", 
  #   "DIABLO" = "^diablo_"
  # ),
  # patterns_to_methods = c("sPLS-DA", "sPLS", "sO2PLS", "MOFA", "DIABLO"),
  target_exclude_patterns = "^mofa_unsupervised"
)

ggsave(
  here("output/figure_running_times.pdf"),
  plot = p,
  width = 6,
  height = 4.5,
  units = "in",
  dpi = dpi
)

## Supplementary tables
write_suppl_tables(
  c(
    list("Prefiltering" = tar_read(suppl_table_prefiltering)),
    tar_read(suppl_tables_integration)
  ),
  output_file = here("output/supplementary_tables.xlsx")
)

