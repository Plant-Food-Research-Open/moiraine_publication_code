library(targets)
library(tarchetypes)
library(moiraine)
library(here)
library(autometric)

tar_option_set(
  packages = c(
    "moiraine",
    "MOFA2",
    "mixOmics",
    "readr",
    "tibble",
    "tidyr",
    "dplyr",
    "ggplot2",
    "patchwork"
  )
)

if (tar_active()) {
  log_start(path = "autometric.log", seconds = 1)
}

list(
  ##==============##
  ## Data loading ----
  ##==============##
  
  ## Data import using a target factory
  moiraine::import_dataset_csv_factory(
    files = c(
      here::here("data/eatrisplus_mirnaseq_mature_data.csv"),
      here::here("data/eatrisplus_mrnaseq2_data.csv"),
      here::here("data/eatrisplus_proteomics_data.csv"),
      here::here("data/eatrisplus_emseq_data.csv"),
      here::here("data/eatrisplus_lipidomics_negative_data.csv"),
      here::here("data/eatrisplus_lipidomics_positive_data.csv"),
      here::here("data/eatrisplus_acylcarnitines_data.csv"),
      here::here("data/eatrisplus_amino_acids_data.csv")
    ),
    col_ids = rep("feature_id", 8),
    features_as_rowss = rep(TRUE, 8),
    target_name_suffixes = c("mirnaseq", "mrnaseq", "proteo", "emseq", 
                             "lipido_neg", "lipido_pos", "acylca", "aa")
  ),
  
  ## Features metadata import
  moiraine::import_fmetadata_csv_factory(
    files = c(
      here::here("data/eatrisplus_mirnaseq_mature_features_metadata.csv"),
      here::here("data/eatrisplus_mrnaseq2_features_metadata.csv"),
      here::here("data/eatrisplus_proteomics_features_metadata.csv"),
      here::here("data/eatrisplus_emseq_features_metadata.csv"),
      here::here("data/eatrisplus_lipidomics_negative_features_metadata.csv"),
      here::here("data/eatrisplus_lipidomics_positive_features_metadata.csv"),
      here::here("data/eatrisplus_acylcarnitines_features_metadata.csv"),
      here::here("data/eatrisplus_amino_acids_features_metadata.csv")
    ),
    col_ids = c("feature_id"),
    target_name_suffixes = c("mirnaseq", "mrnaseq", "proteo", "emseq", 
                             "lipido_neg", "lipido_pos", "acylca", "aa")
  ),
  
  ## Samples metadata import
  moiraine::import_smetadata_csv_factory(
    files = here::here("data/eatrisplus_samples_metadata.csv"),
    col_ids = "subject_id",
    target_name_suffixes = "all"
  ),
  
  ## Creating omics sets for each dataset
  moiraine::create_omics_set_factory(
    datasets = c(data_mirnaseq, data_mrnaseq, data_proteo, data_emseq, 
                 data_lipido_neg, data_lipido_pos, data_acylca, data_aa),
    omics_types = c("transcriptomics", "transcriptomics", "phenomics", "phenomics",
                    "metabolomics", "metabolomics", "metabolomics", "metabolomics"),
    features_metadatas = c(fmetadata_mirnaseq, fmetadata_mrnaseq, fmetadata_proteo, fmetadata_emseq, 
                           fmetadata_lipido_neg, fmetadata_lipido_pos, fmetadata_acylca, fmetadata_aa),
    samples_metadatas = c(smetadata_all, smetadata_all, smetadata_all, smetadata_all,
                          smetadata_all, smetadata_all, smetadata_all, smetadata_all)
  ),

  ## Creating the MultiDataSet object
  tar_target(
    mo_set,
    moiraine::create_multiomics_set(
      list(set_mirnaseq, set_mrnaseq, set_proteo, set_emseq, 
           set_lipido_neg, set_lipido_pos, set_acylca, set_aa),
      datasets_names = c("mirnaseq", "mrnaseq", "proteo", "emseq", 
                         "lipido_neg", "lipido_pos", "acylca", "aa")
    )
  ),
  
  ##=========================##
  ## Datasets transformation ----
  ##=========================##

  tar_target(
    density_plot,
    plot_density_data(mo_set, combined = FALSE, scales = "free")
  ),
  
  # Applying transformations to the datasets
  moiraine::transformation_datasets_factory(
    mo_set,
    c("phenotypes+proteo" = "best-normalize-manual",
      "metabolome+lipido_neg" = "best-normalize-manual",
      "metabolome+lipido_pos" = "best-normalize-manual"),
    methods = "center_scale",
    transformed_data_name = "mo_set_transformed"
  ),

  tar_target(
    density_plot_transformed,
    plot_density_data(mo_set_transformed, combined = FALSE, scales = "free")
  ),

  ##============================================##
  ## Individual PCA and missing data imputation ----
  ##============================================##

  ## Running a PCA on each dataset
  moiraine::pca_complete_data_factory(
    mo_set_transformed,
    complete_data_name = "mo_set_complete"
  ),
  
  ##===============##
  ## Pre-filtering ----
  ##===============##

  ## Supervised feature selection based on disease status
  moiraine::feature_preselection_splsda_factory(
    mo_set_complete,
    group = "sex",
    to_keep_ns = c("rnaseq+mirnaseq" = 1000,
                   "rnaseq+mrnaseq" = 1000,
                   "phenotypes+proteo" = 1000,
                   "phenotypes+emseq" = 1000),
    filtered_set_target_name = "mo_presel_supervised",
    seed_perf = c("rnaseq+mirnaseq" = 142,
                  "rnaseq+mrnaseq" = 143,
                  "phenotypes+proteo" = 144,
                  "phenotypes+emseq" = 145),
    ncomp_max = 2
  ),

  tar_target(
    individual_splsda_perf_plot,
    moiraine::plot_feature_preselection_splsda(individual_splsda_perf)
  ),

  ##=================##
  ## DIABLO pipeline ----
  ##=================##

  ## Creating the DIABLO input
  tar_target(
    diablo_input,
    moiraine::get_input_mixomics_supervised(
      mo_presel_supervised,
      group = "sex"
    )
  ),

  ## Running sPLS on each dataset to construct the design matrix
  moiraine::diablo_pairwise_pls_factory(diablo_input),

  ## Initial DIABLO run with no feature selection and large number of components
  tar_target(
    diablo_novarsel,
    moiraine::diablo_run(
      diablo_input,
      diablo_design_matrix,
      ncomp = 4
    )
  ),

  ## Cross-validation for number of components
  tar_target(
    diablo_perf_res,
    mixOmics::perf(
      diablo_novarsel,
      validation = "Mfold",
      folds = 10,
      nrepeat = 10,
      cpus = 3
    )
  ),

  ## Plotting cross-validation results (for number of components)
  tar_target(
    diablo_perf_plot,
    moiraine::diablo_plot_perf(diablo_perf_res)
  ),

  ## Selected value for ncomp
  tar_target(
    diablo_optim_ncomp,
    moiraine::diablo_get_optim_ncomp(diablo_perf_res)
  ),

  ## Cross-validation for number of features to retain
  tar_target(
    diablo_tune_res,
    moiraine::diablo_tune(
      diablo_input,
      diablo_design_matrix,
      ncomp = diablo_optim_ncomp,
      keepX_list = purrr::map(names(mo_set), \(x) c(10, 30)),
      validation = "Mfold",
      folds = 10,
      nrepeat = 5,
      dist = "centroids.dist",
      cpus = 3,
      seed = 9354
    )
  ),

  ## Plotting cross-validation results (for number of features)
  tar_target(
    diablo_tune_plot,
    moiraine::diablo_plot_tune(diablo_tune_res)
  ),

  tar_target(
    diablo_tune_table,
    diablo_table_optim_keepX(diablo_tune_res)
  ),

  ## Final DIABLO run
  tar_target(
    diablo_final_run,
    moiraine::diablo_run(
      diablo_input,
      diablo_design_matrix,
      ncomp = diablo_optim_ncomp,
      keepX = diablo_tune_res$choice.keepX
    )
  ),

  tar_target(
    diablo_output,
    moiraine::get_output(diablo_final_run)
  ),

  # ##===============##
  # ## sPLS pipeline ----
  # ##===============##
  # 
  # ## Creating sPLS input
  # tar_target(
  #   spls_input,
  #   moiraine::get_input_spls(
  #     mo_presel_supervised,
  #     mode = "canonical",
  #     datasets = c("rnaseq", "metabolome")
  #   )
  # ),
  # 
  # ## Initial PLS run with no feature selection and large number of components
  # tar_target(
  #   spls_novarsel,
  #   moiraine::spls_run(
  #     spls_input,
  #     ncomp = 4
  #   )
  # ),
  # 
  # ## Cross-validation for number of components
  # tar_target(
  #   spls_perf_res,
  #   mixOmics::perf(
  #     spls_novarsel,
  #     validation = "Mfold",
  #     folds = 10,
  #     nrepeat = 10,
  #     cpus = 3
  #   )
  # ),
  # 
  # ## Plotting cross-validation results (for number of components)
  # ## Can try criterion = 'Q2.total', 'cor.tpred', 'cor.upred', 'RSS.tpred',
  # ## 'RSS.upred' (but avoid 'RSS' and 'PRESS')
  # tar_target(
  #   spls_perf_plot,
  #   plot(spls_perf_res, criterion = "Q2.total")
  # ),
  # 
  # ## Selected value for ncomp
  # tar_target(
  #   spls_optim_ncomp,
  #   moiraine::spls_get_optim_ncomp(spls_perf_res, min_ncomp = 2)
  # ),
  # 
  # ## Cross-validation for number of features to retain
  # tar_target(
  #   spls_tune_res,
  #   moiraine::spls_tune(
  #     spls_input,
  #     ncomp = spls_optim_ncomp,
  #     keepX = seq(10, 100, 10),
  #     keepY = seq(10, 100, 10),
  #     validation = "Mfold",
  #     folds = 10,
  #     nrepeat = 5,
  #     measure = "cor",
  #     cpus = 3,
  #     seed = -584594170
  #   )
  # ),
  # 
  # ## Plotting cross-validation results (for number of features)
  # tar_target(
  #   spls_tune_plot,
  #   moiraine::spls_plot_tune(spls_tune_res)
  # ),
  # 
  # tar_target(
  #   spls_tune_table,
  #   spls_table_optim_keepX(spls_tune_res)
  # ),
  # 
  # ## Final sPLS run
  # tar_target(
  #   spls_final_run,
  #   moiraine::spls_run(
  #     spls_input,
  #     ncomp = spls_optim_ncomp,
  #     keepX = spls_tune_res$choice.keepX,
  #     keepY = spls_tune_res$choice.keepY
  #   )
  # ),
  # 
  # tar_target(
  #   spls_output,
  #   moiraine::get_output(spls_final_run)
  # ),
  # 
  # 
  # ##=================##
  # ## sO2PLS pipeline ----
  # ##=================##
  # 
  # ## Creating sO2PLS input
  # tar_target(
  #   omicspls_input,
  #   moiraine::get_input_omicspls(
  #     mo_presel_supervised,
  #     datasets = c("rnaseq", "metabolome")
  #   )
  # ),
  # 
  # ## Adjusted cross-validation for number of components
  # tar_target(
  #   so2pls_cv_adj,
  #   moiraine::so2pls_crossval_o2m_adjR2(
  #     omicspls_input,
  #     a = 1:5,
  #     ax = seq(0, 10, by = 2),
  #     ay = seq(0, 10, by = 2),
  #     nr_folds = 10,
  #     nr_cores = 6,
  #     seed = 127
  #   )
  # ),
  # tar_target(
  #   so2pls_cv_adj_res,
  #   moiraine::so2pls_get_optim_ncomp_adj(so2pls_cv_adj)
  # ),
  # 
  # ## Plotting adjusted cross-validation results
  # tar_target(
  #   so2pls_cv_adj_plot,
  #   moiraine::so2pls_plot_cv_adj(so2pls_cv_adj)
  # ),
  # 
  # ## Standard cross-validation for number of components
  # tar_target(
  #   so2pls_cv,
  #   moiraine::so2pls_crossval_o2m(
  #     omicspls_input,
  #     so2pls_cv_adj,
  #     nr_folds = 10,
  #     nr_cores = 6,
  #     seed = 356
  #   )
  # ),
  # tar_target(
  #   so2pls_cv_res,
  #   moiraine::so2pls_get_optim_ncomp(so2pls_cv)
  # ),
  # 
  # ## Plotting standard cross-validation results
  # tar_target(
  #   so2pls_cv_plot,
  #   moiraine::so2pls_plot_cv(so2pls_cv)
  # ),
  # 
  # ## Cross-validation for sparsity parameters
  # tar_target(
  #   so2pls_cv_sparsity,
  #   moiraine::so2pls_crossval_sparsity(
  #     omicspls_input,
  #     n = so2pls_cv_res["n"],
  #     nx = so2pls_cv_res["nx"],
  #     ny = so2pls_cv_res["ny"],
  #     nr_folds = 10,
  #     keepx_seq = c(seq(5, 30, 5), seq(40, 100, 10)),
  #     keepy_seq = c(seq(5, 40, 5)),
  #     seed = -1138855226
  #   )
  # ),
  # tar_target(
  #   so2pls_cv_sparsity_res,
  #   moiraine::so2pls_get_optim_keep(so2pls_cv_sparsity)
  # ),
  # 
  # ## Plotting the results of the cross-validation for the number of features
  # ## to retain from each dataset for the different joint components
  # tar_target(
  #   so2pls_cv_sparsity_plot,
  #   moiraine::so2pls_plot_cv_sparsity(so2pls_cv_sparsity)
  # ),
  # 
  # ## Extracting sparsity results in table format
  # tar_target(
  #   so2pls_cv_sparsity_table,
  #   moiraine::so2pls_print_cv_sparsity(so2pls_cv_sparsity_res)
  # ),
  # 
  # ## Final sO2PLS run
  # tar_target(
  #   so2pls_final_run,
  #   moiraine::so2pls_o2m(
  #     omicspls_input,
  #     so2pls_cv_res,
  #     so2pls_cv_sparsity_res
  #   )
  # ),
  # 
  # tar_target(
  #   so2pls_output,
  #   moiraine::get_output(so2pls_final_run)
  # ),
  # 
  
  ##===============##
  ## MOFA pipeline ----
  ##===============##

  ## Using supervised prefiletered dataset ----

  ## Creating MOFA input
  tar_target(
    mofa_input,
    moiraine::get_input_mofa(
      mo_presel_supervised,
      options_list = list(
        data_options = list(scale_views = TRUE),
        training_options = list(seed = 43)
      ),
      only_common_samples = FALSE
    )
  ),

  ## Training MOFA model
  tar_target(
    mofa_trained,
    MOFA2::run_mofa(
      mofa_input,
      save_data = TRUE,
      use_basilisk = TRUE
    )
  ),

  tar_target(
    mofa_output,
    moiraine::get_output(mofa_trained)
  )


  # ##====================##
  # ## Results comparison ----
  # ##====================##
  # 
  # ## List of formatted output
  # tar_target(
  #   output_list,
  #   list(spls_output, so2pls_output, mofa_output, diablo_output)
  # ),
  # 
  # tar_target(
  #   output_list_mofa,
  #   list(
  #     "MOFA (supervised pref.)" = mofa_output,
  #     "MOFA (unsupervised pref.)" = mofa_unsupervised_output
  #   )
  # )
)