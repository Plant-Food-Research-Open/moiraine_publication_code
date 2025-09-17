library(targets)
library(tarchetypes)
library(moiraine)
library(here)
library(autometric)

tar_option_set(
  packages = c(
    "here",
    "moiraine",
    "mixOmics"
  )
)

source(here("R/pipeline_helper_functions.R"))

if (tar_active()) {
  log_start(path = "autometric.log", seconds = 1)
}

list(
  ##==============##
  ## Data loading ----
  ##==============##
  
  ## Importing the full genomics dataset
  tar_target(
    dataset_file_geno_full,
    here::here("data/li_genomics_data.csv"),
    format = "file"
  ),
  tar_target(
    data_geno_full,
    import_dataset_csv(
      dataset_file_geno_full, 
      col_id = "marker", 
      features_as_rows = TRUE
    )
  ),
  
  tar_target(
    fmetadata_file_geno_full,
    here::here("data/li_genomics_features_metadata.csv"),
    format = "file"
  ),
  tar_target(
    fmetadata_geno_full,
    import_fmetadata_csv(
      fmetadata_file_geno_full,
      col_id = "marker",
      col_types = c("chromosome" = "c")
    )
  ),
  
  ## Importing the reduced genomics dataset
  tar_target(
    dataset_file_geno_reduced,
    system.file("extdata/genomics_dataset.csv", package = "moiraine"),
    format = "file"
  ),
  tar_target(
    data_geno_reduced,
    import_dataset_csv(
      dataset_file_geno_reduced, 
      col_id = "marker", 
      features_as_rows = TRUE
    )
  ),
  
  tar_target(
    fmetadata_file_geno_reduced,
    system.file("extdata/genomics_features_info.csv", package = "moiraine"),
    format = "file"
  ),
  tar_target(
    fmetadata_geno_reduced,
    import_fmetadata_csv(
      fmetadata_file_geno_reduced,
      col_id = "marker",
      col_types = c("chromosome" = "c")
    )
  ),
  
  ## Samples metadata import
  moiraine::import_smetadata_csv_factory(
    files = system.file("extdata/samples_info.csv", package = "moiraine"),
    col_ids = "animal_id",
    target_name_suffixes = "all"
  ),
  
  ## Creating omics sets for each dataset
  tar_target(
    set_geno_full,
    create_omics_set(
      data_geno_full,
      omics_type = "genomics",
      features_metadata = fmetadata_geno_full,
      samples_metadata = smetadata_all
    )
  ),
  tar_target(
    set_geno_reduced,
    create_omics_set(
      data_geno_reduced,
      omics_type = "genomics",
      features_metadata = fmetadata_geno_reduced,
      samples_metadata = smetadata_all
    )
  ),
  
  ## Creating the MultiDataSet object
  tar_target(
    mo_set_full,
    moiraine::create_multiomics_set(list(set_geno_full))
  ),
  tar_target(
    mo_set_reduced,
    moiraine::create_multiomics_set(list(set_geno_reduced))
  ),
  
  ##============================================##
  ## Individual PCA and missing data imputation ----
  ##============================================##
  
  ## Running a PCA on each dataset
  moiraine::pca_complete_data_factory(
    mo_set_full,
    complete_data_name = "mo_set_complete_full",
    target_name_prefix = "full_"
  ),
  moiraine::pca_complete_data_factory(
    mo_set_reduced,
    complete_data_name = "mo_set_complete_reduced",
    target_name_prefix = "reduced_"
  ),
  
  ##===============##
  ## Pre-filtering ----
  ##===============##
  
  ## Supervised feature selection based on disease status
  moiraine::feature_preselection_splsda_factory(
    mo_set_complete_full,
    group = "status",
    to_keep_ns = c("snps" = 1000),
    filtered_set_target_name = "mo_presel_supervised_full",
    seed_perf = c("snps" = -1591100874),
    target_name_prefix = "full_"
  ),
  moiraine::feature_preselection_splsda_factory(
    mo_set_complete_reduced,
    group = "status",
    to_keep_ns = c("snps" = 1000),
    filtered_set_target_name = "mo_presel_supervised_reduced",
    seed_perf = c("snps" = 5843266),
    target_name_prefix = "reduced_"
  )
)