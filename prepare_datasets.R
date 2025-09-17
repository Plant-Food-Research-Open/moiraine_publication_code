library(MultiAssayExperiment)
library(tidyverse)
library(janitor)
library(fs)
library(here)
library(withr)

dir_temp <- tempdir()
fs::dir_create("data")

# TCGA dataset ------------------------------------------------------------

tcga_file <- fs::path(dir_temp, "brca.tcga.RData")
download.file(
  "https://github.com/xlucpu/MOVICS/raw/refs/heads/master/inst/extdata/brca.tcga.RData",
  tcga_file
)

load(tcga_file)

# mRNA dataset
tcga_list <- brca.tcga[c("count", "lncRNA.expr", "meth.beta", "mut.status")] |> 
  purrr::map(\(x) tibble::as_tibble(x, rownames = "feature_id")) |> 
  rlang::set_names("rnaseq", "lncrna", "meth", "mut")


# EATRIS-plus dataset -----------------------------------------------------

mae_file <- fs::path(dir_temp, "mae_mae.rds")
download.file(
  "https://zenodo.org/records/10782800/files/mae_mae.rds", 
  mae_file
)

options(timeout = 600)
mae_h5_file <- fs::path(dir_temp, "mae_experiments.h5")
download.file(
  "https://zenodo.org/records/10782800/files/mae_experiments.h5", 
  mae_h5_file,
  mode = "wb"
)

## needed to be able to read the data in memory
fs::file_copy(mae_h5_file, getwd())

eatrisplus_mae <- readRDS(mae_file)

omics_layers <- c(
  "Acylcarnitines | batch-adjusted | imputed missing values",
  "Amino acids | batch-adjusted | imputed missing values",
  "Lipidomics, positive | transformed",
  "Lipidomics, negative | transformed",
  "Proteomics | imputed missing values",
  "mRNA-seq-2 | batch-adjusted | corrected",
  "miRNA-seq, mature | batch-adjusted | cell-type adjusted",
  "EM-seq | cell-type adjusted | 100,000 most variable CpG sites"
)
omics_names <- stringr::str_extract(omics_layers, "[^\\|]+(?= \\|)") |> 
  stringr::str_remove_all("-") |> 
  stringr::str_replace_all("\\s|(, )", "_") |> 
  stringr::str_to_lower()

eatrisplus_mae <- eatrisplus_mae[, , omics_layers]

smeta <- MultiAssayExperiment::colData(eatrisplus_mae) |> 
  tibble::as_tibble() |> 
  janitor::clean_names()

data_list <- MultiAssayExperiment::assays(eatrisplus_mae) |> 
  purrr::map(\(x) tibble::as_tibble(x, rownames = "feature_id")) |> 
  rlang::set_names(omics_names)


fmeta_list <- eatrisplus_mae |> 
  purrr::map(\(x) {
    SummarizedExperiment::rowData(x) |>
      tibble::as_tibble(rownames = "feature_id") |> 
      janitor::clean_names()
  }) |> 
  rlang::set_names(omics_names)


readr::write_csv(
  smeta, 
  here::here(stringr::str_glue("data/eatrisplus_samples_metadata.csv"))
)
purrr::walk(
  omics_names,
  \(i) {
    readr::write_csv(
      data_list[[i]], 
      here::here(stringr::str_glue("data/eatrisplus_{i}_data.csv"))
    )
    write_csv(
      fmeta_list[[i]], 
      here::here(stringr::str_glue("data/eatrisplus_{i}_features_metadata.csv"))
    )
  }
)

fs::file_delete(basename(mae_h5_file))

# Li, 2022 full genomics dataset ------------------------------------------

geno_file <- fs::path(dir_temp, "genotype.csv")
download.file(
  "https://borealisdata.ca/api/access/datafile/411517",
  geno_file
)

snp_file <- fs::path(dir_temp, "snp_info.txt")
download.file(
  "https://borealisdata.ca/api/access/datafile/411514",
  snp_file
)

## Genotypes features metadata
geno_fmeta_df <- readr::read_tsv(
  snp_file, 
  name_repair = make_clean_names,
  show_col_types = FALSE
) |>
  dplyr::select(-index) |>
  dplyr::rename(marker = name) |>
  dplyr::mutate(
    snp = stringr::str_remove_all(snp, "\\[|\\]"),
    chromosome = as.character(chromosome)
  ) |>
  tidyr::separate_wider_delim(
    cols = snp,
    delim = "/",
    names = c("ref", "alt")
  )

dosage_values <- c("AA" = 0, "AB" = 1, "BB" = 2)

geno_df <- readr::read_csv(geno_file, show_col_types = FALSE) |>
  janitor::remove_empty("cols") |>
  dplyr::select(-Neogen.ID, -ID_RNA) |>
  dplyr::rename(sample_id = Submitted.ID) |>
  tidyr::pivot_longer(
    cols = -sample_id,
    names_to = "marker",
    values_to = "dosage"
  ) |>
  dplyr::mutate(
    dosage = dosage_values[dosage]
  ) |>
  tidyr::pivot_wider(
    names_from = sample_id,
    values_from = dosage
  )

## Filtering according to the paper; ie removing markers with:
## - > 10% more missing values
## - minor allele frequency < 5%
## - located on sex chromosomes
markers_filtering <- geno_df |>
  dplyr::left_join(
    dplyr::select(geno_fmeta_df, marker, chromosome),
    by = "marker"
  ) |>
  tidyr::pivot_longer(
    cols = -c(marker, chromosome),
    names_to = "sample_id",
    values_to = "dosage"
  ) |>
  dplyr::group_by(marker, chromosome) |>
  dplyr::summarise(
    n_tot = dplyr::n(),
    missing = sum(is.na(dosage)),
    mas = sum(dosage, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    frac_missing = missing / n_tot,
    maf = mas / (2 * (n_tot - missing)),
    maf = dplyr::case_when(
      maf > 0.5 ~ 1 - maf,
      TRUE ~ maf
    ),
    sex_chr = chromosome %in% c("MT", "X", "Y"),
    retained = (frac_missing <= 0.1) & (maf >= 0.05) & !sex_chr
  ) |>
  dplyr::filter(retained) |>
  dplyr::pull(marker)


geno_df <- geno_df |>
  dplyr::filter(marker %in% markers_filtering)
geno_fmeta_df <- geno_fmeta_df |>
  dplyr::filter(marker %in% markers_filtering)

readr::write_csv(
  geno_df, 
  here::here(stringr::str_glue("data/li_genomics_data.csv"))
)
readr::write_csv(
  geno_fmeta_df, 
  here::here(stringr::str_glue("data/li_genomics_features_metadata.csv"))
)
