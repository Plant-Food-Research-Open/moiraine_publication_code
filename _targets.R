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

source(here("R/pipeline_helper_functions.R"))

#unlink("autometric.log")
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
      system.file("extdata/genomics_dataset.csv", package = "moiraine"),
      system.file("extdata/transcriptomics_dataset.csv", package = "moiraine"),
      system.file("extdata/metabolomics_dataset.csv", package = "moiraine")
    ),
    col_ids = c("marker", "gene_id", "sample_id"),
    features_as_rowss = c(TRUE, TRUE, FALSE),
    target_name_suffixes = c("geno", "transcripto", "metabo")
  ),

  ## Features metadata import
  tar_target(
    fmetadata_file_geno,
    system.file("extdata/genomics_features_info.csv", package = "moiraine"),
    format = "file"
  ),

  tar_target(
    fmetadata_geno,
    moiraine::import_fmetadata_csv(
      fmetadata_file_geno,
      col_id = "marker",
      col_types = c("chromosome" = "c")
    )
  ),

  moiraine::import_fmetadata_csv_factory(
    files = c(
      system.file("extdata/metabolomics_features_info.csv", package = "moiraine")
    ),
    col_ids = c("feature_id"),
    target_name_suffixes = c("metabo")
  ),

  moiraine::import_fmetadata_gff_factory(
    files = system.file("extdata/bos_taurus_gene_model.gff3", package = "moiraine"),
    feature_types = "genes",
    add_fieldss = c("Name", "description"),
    target_name_suffixes = "transcripto"
  ),

  ## Samples metadata import
  moiraine::import_smetadata_csv_factory(
    files = system.file("extdata/samples_info.csv", package = "moiraine"),
    col_ids = "animal_id",
    target_name_suffixes = "all"
  ),

  ## Creating omics sets for each dataset
  moiraine::create_omics_set_factory(
    datasets = c(data_geno, data_transcripto, data_metabo),
    omics_types = c("genomics", "transcriptomics", "metabolomics"),
    features_metadatas = c(fmetadata_geno, fmetadata_transcripto, fmetadata_metabo),
    samples_metadatas = c(smetadata_all, smetadata_all, smetadata_all)
  ),

  ## Creating the MultiDataSet object
  tar_target(
    mo_set,
    moiraine::create_multiomics_set(
      list(set_geno,
           set_transcripto,
           set_metabo)
    )
  ),

  ##===================================##
  ## Adding information about features ----
  ##===================================##

  ## RNAseq differential expression results file
  tar_target(
    rnaseq_de_res_file,
    system.file(
      "extdata/transcriptomics_de_results.csv",
      package = "moiraine"
    ),
    format = "file"
  ),

  ## Reading the RNAseq differential expression results
  tar_target(
    rnaseq_de_res_df,
    readr::read_csv(rnaseq_de_res_file) |>
      dplyr::rename(feature_id = gene_id) |>
      dplyr::mutate(dataset = "rnaseq")
  ),

  ## Adding the differential expression results to the MultiDataSet object
  tar_target(
    mo_set_de,
    moiraine::add_features_metadata(mo_set, rnaseq_de_res_df)
  ),

  ##=========================##
  ## Datasets transformation ----
  ##=========================##

  # Applying transformations to the datasets
  moiraine::transformation_datasets_factory(
    mo_set_de,
    c("rnaseq" = "vst-deseq2",
      "metabolome" = "logx"),
    log_bases = 2,
    pre_log_functions = zero_to_half_min,
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

  ## Unsupervised feature selection based on MAD score
  moiraine::feature_preselection_mad_factory(
    mo_set_complete,
    to_keep_ns = c("snps" = 1000, "rnaseq" = 1000),
    with_ties = TRUE,
    filtered_set_target_name = "mo_presel_unsupervised"
  ),

  tar_target(
    individual_mad_plot,
    moiraine::plot_feature_preselection_mad(individual_mad_values)
  ),

  ## Supervised feature selection based on disease status
  moiraine::feature_preselection_splsda_factory(
    mo_set_complete,
    group = "status",
    to_keep_ns = c("snps" = 1000, "rnaseq" = 1000),
    filtered_set_target_name = "mo_presel_supervised",
    seed_perf = c("snps" = -1591100874, "rnaseq" = 1791752001)
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
      group = "status"
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
      ncomp = 7
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
      validation = "Mfold",
      folds = 10,
      nrepeat = 5,
      dist = "centroids.dist",
      cpus = 3,
      seed = 1659021768
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

  ##===============##
  ## sPLS pipeline ----
  ##===============##

  ## Creating sPLS input
  tar_target(
    spls_input,
    moiraine::get_input_spls(
      mo_presel_supervised,
      mode = "canonical",
      datasets = c("rnaseq", "metabolome")
    )
  ),

  ## Initial PLS run with no feature selection and large number of components
  tar_target(
    spls_novarsel,
    moiraine::spls_run(
      spls_input,
      ncomp = 4
    )
  ),

  ## Cross-validation for number of components
  tar_target(
    spls_perf_res,
    mixOmics::perf(
      spls_novarsel,
      validation = "Mfold",
      folds = 10,
      nrepeat = 10,
      cpus = 3
    )
  ),

  ## Plotting cross-validation results (for number of components)
  ## Can try criterion = 'Q2.total', 'cor.tpred', 'cor.upred', 'RSS.tpred',
  ## 'RSS.upred' (but avoid 'RSS' and 'PRESS')
  tar_target(
    spls_perf_plot,
    plot(spls_perf_res, criterion = "Q2.total")
  ),

  ## Selected value for ncomp
  tar_target(
    spls_optim_ncomp,
    moiraine::spls_get_optim_ncomp(spls_perf_res, min_ncomp = 2)
  ),

  ## Cross-validation for number of features to retain
  tar_target(
    spls_tune_res,
    moiraine::spls_tune(
      spls_input,
      ncomp = spls_optim_ncomp,
      keepX = seq(10, 100, 10),
      keepY = seq(10, 100, 10),
      validation = "Mfold",
      folds = 10,
      nrepeat = 5,
      measure = "cor",
      cpus = 3,
      seed = -584594170
    )
  ),

  ## Plotting cross-validation results (for number of features)
  tar_target(
    spls_tune_plot,
    moiraine::spls_plot_tune(spls_tune_res)
  ),

  tar_target(
    spls_tune_table,
    spls_table_optim_keepX(spls_tune_res)
  ),

  ## Final sPLS run
  tar_target(
    spls_final_run,
    moiraine::spls_run(
      spls_input,
      ncomp = spls_optim_ncomp,
      keepX = spls_tune_res$choice.keepX,
      keepY = spls_tune_res$choice.keepY
    )
  ),

  tar_target(
    spls_output,
    moiraine::get_output(spls_final_run)
  ),


  ##=================##
  ## sO2PLS pipeline ----
  ##=================##

  ## Creating sO2PLS input
  tar_target(
    omicspls_input,
    moiraine::get_input_omicspls(
      mo_presel_supervised,
      datasets = c("rnaseq", "metabolome")
    )
  ),

  ## Adjusted cross-validation for number of components
  tar_target(
    so2pls_cv_adj,
    moiraine::so2pls_crossval_o2m_adjR2(
      omicspls_input,
      a = 1:5,
      ax = seq(0, 10, by = 2),
      ay = seq(0, 10, by = 2),
      nr_folds = 10,
      nr_cores = 6,
      seed = 127
    )
  ),
  tar_target(
    so2pls_cv_adj_res,
    moiraine::so2pls_get_optim_ncomp_adj(so2pls_cv_adj)
  ),

  ## Plotting adjusted cross-validation results
  tar_target(
    so2pls_cv_adj_plot,
    moiraine::so2pls_plot_cv_adj(so2pls_cv_adj)
  ),

  ## Standard cross-validation for number of components
  tar_target(
    so2pls_cv,
    moiraine::so2pls_crossval_o2m(
      omicspls_input,
      so2pls_cv_adj,
      nr_folds = 10,
      nr_cores = 6,
      seed = 356
    )
  ),
  tar_target(
    so2pls_cv_res,
    moiraine::so2pls_get_optim_ncomp(so2pls_cv)
  ),

  ## Plotting standard cross-validation results
  tar_target(
    so2pls_cv_plot,
    moiraine::so2pls_plot_cv(so2pls_cv)
  ),

  ## Cross-validation for sparsity parameters
  tar_target(
    so2pls_cv_sparsity,
    moiraine::so2pls_crossval_sparsity(
      omicspls_input,
      n = so2pls_cv_res["n"],
      nx = so2pls_cv_res["nx"],
      ny = so2pls_cv_res["ny"],
      nr_folds = 10,
      keepx_seq = c(seq(5, 30, 5), seq(40, 100, 10)),
      keepy_seq = c(seq(5, 40, 5)),
      seed = -1138855226
    )
  ),
  tar_target(
    so2pls_cv_sparsity_res,
    moiraine::so2pls_get_optim_keep(so2pls_cv_sparsity)
  ),

  ## Plotting the results of the cross-validation for the number of features
  ## to retain from each dataset for the different joint components
  tar_target(
    so2pls_cv_sparsity_plot,
    moiraine::so2pls_plot_cv_sparsity(so2pls_cv_sparsity)
  ),

  ## Extracting sparsity results in table format
  tar_target(
    so2pls_cv_sparsity_table,
    moiraine::so2pls_print_cv_sparsity(so2pls_cv_sparsity_res)
  ),

  ## Final sO2PLS run
  tar_target(
    so2pls_final_run,
    moiraine::so2pls_o2m(
      omicspls_input,
      so2pls_cv_res,
      so2pls_cv_sparsity_res
    )
  ),

  tar_target(
    so2pls_output,
    moiraine::get_output(so2pls_final_run)
  ),

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
        model_options = list( likelihoods = c(
          "snps" = "poisson",
          "rnaseq" = "gaussian",
          "metabolome" = "gaussian"
        )
        ),
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
  ),

  ## Using unsupervised prefiltered dataset ----

  ## Creating the input object for the MOFA pipeline
  ## using the unsupervised preselection results
  tar_target(
    mofa_unsupervised_input,
    moiraine::get_input_mofa(
      mo_presel_unsupervised,
      options_list = list(
        data_options = list(scale_views = TRUE),
        model_options = list(
          likelihoods = c(
            "snps" = "poisson",
            "rnaseq" = "gaussian",
            "metabolome" = "gaussian")
        ),
        training_options = list(seed = 72)
      ),
      only_common_samples = FALSE
    )
  ),

  ## Training the model with the MOFA algorithm
  tar_target(
    mofa_unsupervised_trained,
    MOFA2::run_mofa(
      mofa_unsupervised_input,
      save_data = TRUE,
      use_basilisk = TRUE
    )
  ),

  ## Formatting MOFA output
  tar_target(
    mofa_unsupervised_output,
    moiraine::get_output(mofa_unsupervised_trained)
  ),

  ##====================##
  ## Results evaluation ----
  ##====================##

  ## Reading GO annotation file
  tar_target(
    rnaseq_go_terms_file,
    system.file(
      "extdata/transcriptomics_go_annotation.csv",
      package = "moiraine"
    ),
    format = "file"
  ),

  tar_target(
    rnaseq_go_df,
    readr::read_csv(rnaseq_go_terms_file) |>
      filter(go_domain == "Biological process")
  ),

  ## Making GO terms sets
  tar_target(
    go_sets,
    moiraine::make_feature_sets_from_df(
      rnaseq_go_df,
      col_id = "gene_id",
      col_set = "go_id"
    )
  ),

  ## Filtering GO term sets against measured features
  tar_target(
    go_sets_filtered,
    moiraine::reduce_feature_sets_data(go_sets, mo_set_complete)
  ),

  ## Checking genes GO term sets against datasets
  tar_target(
    go_sets_check,
    moiraine::check_feature_sets(
      go_sets_filtered,
      mo_set_complete,
      datasets = "rnaseq"
    )
  ),

  ## Table of information about GO terms
  tar_target(
    go_sets_info,
    rnaseq_go_df |>
      dplyr::select(go_id, go_name) |>
      dplyr::distinct()
  ),

  ## MOFA latent components enrichment analysis
  tar_target(
    mofa_enrichment_results,
    moiraine::evaluate_method_enrichment(
      mofa_output,
      go_sets_filtered,
      datasets = "rnaseq",
      latent_dimensions = "Factor 1",
      use_abs = TRUE,
      min_set_size = 10,
      add_missing_features = TRUE,
      mo_data = mo_set_complete,
      sets_info_df = go_sets_info,
      col_set = "go_id"
    )
  ),


  ##====================##
  ## Results comparison ----
  ##====================##

  ## List of formatted output
  tar_target(
    output_list,
    list(spls_output, so2pls_output, mofa_output, diablo_output)
  ),

  tar_target(
    output_list_mofa,
    list(
      "MOFA (supervised pref.)" = mofa_output,
      "MOFA (unsupervised pref.)" = mofa_unsupervised_output
    )
  ),

  ##===============================##
  ## Manuscript figures and tables ----
  ##===============================##
  
  tar_target(
    suppl_figure_consensus_metrics,
    moiraine::show_consensus_metrics()
  ),

  tar_target(
    suppl_figure_upset,
    moiraine::plot_samples_upset(mo_set)
  ),

  tar_target(
    figure_pca_rnaseq,
    make_composite_pca_plot(pca_runs_list, mo_set_de, "rnaseq", "rnaseq_batch")
  ),

  tar_target(
    suppl_figure_pca_snps,
    make_composite_pca_plot(pca_runs_list, mo_set_de, "snps", "geno_comp_cluster")
  ),

  tar_target(
    suppl_figure_pca_metabo,
    make_composite_pca_plot(pca_runs_list, mo_set_de, "metabolome", "gender")
  ),

  tar_target(
    figure_samples_score_integration,
    make_integration_samples_score_plot(output_list, mo_set_de)
  ),

  tar_target(
    figure_mofa_features,
    make_mofa_features_plot(mofa_output, mo_set_de)
  ),

  tar_target(
    table_mofa_enrichment,
    make_mofa_enrichment_table(mofa_enrichment_results)
  ),

  tar_target(
    figure_methods_comparison,
    make_methods_comparison_plot(output_list, mo_set_de)
  ),

  tar_target(
    suppl_figure_prefiltering_comparison,
    make_prefiltering_comparison_plot(output_list_mofa)
  ),


  ##=================================##
  ## Manuscript supplementary tables ----
  ##=================================##

  tar_target(
    suppl_table_prefiltering,
    make_prefiltering_suppl_table(
      mo_set,
      mo_presel_supervised,
      mo_presel_unsupervised
    )
  ),

  tar_target(
    suppl_tables_integration,
    make_integration_suppl_tables(
      list(
        "sPLS" = spls_output,
        "sO2PLS" = so2pls_output,
        "MOFA" = mofa_output,
        "DIABLO" = diablo_output,
        "MOFA (unsupervised pref.)" = mofa_unsupervised_output
      )
    )
  )
)