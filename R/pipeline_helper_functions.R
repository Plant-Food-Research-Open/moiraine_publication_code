# Manuscript figures and tables -------------------------------------------

status_colors <- c("Control" = "#66C2A5", "BRD" = "#6A3D9A")

#' Formats selected number of features to retain for sPLS into a table
#'
#' @param spls_tune_res List, output of [moiraine::spls_tune()].
#' @returns A tibble with one row per dataset.
spls_table_optim_keepX <- function(spls_tune_res) {
  comps <- paste0("Component ", seq_along(spls_tune_res$choice.keepX))
  tibble::tibble(
    component = rep(comps, 2),
    Dataset = rep(attr(spls_tune_res, "datasets_name"), each = 2),
    value = c(spls_tune_res$choice.keepX, spls_tune_res$choice.keepY)
  ) |>
    tidyr::pivot_wider(
      names_from = component,
      values_from = value
    ) |>
    dplyr::mutate(Total = rowSums(dplyr::across(tidyselect::where(is.numeric))))
}

#' PCA results for an omics dataset
#'
#' Plots the PCA results (screeplot and first 3 PCs) for the one of the omics
#' datasets.
#'
#' @param pca_res List of PCA results for each omics dataset.
#' @param mo_set MultiDataSet object, the multi-omics dataset.
#' @param dataset Character, name of the dataset for which the PCA results
#'   should be shown.
#' @param covariate Character, name of the samples covariate to use when
#'   colouring the samples in the PC space.
#' @returns A patchwork of plots.
make_composite_pca_plot <- function(pca_res, mo_set, dataset, covariate) {

  ## Plot of percentage of variance explained for first 10 PCs
  p_screeplot <- moiraine::plot_screeplot_pca(
    pca_res,
    datasets = dataset
  ) +
    ggplot2::scale_fill_manual(values = "cornflowerblue") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.08)))

  ## Plot of samples in the space spanned by the first 3 PCs
  p_pcs <- moiraine::plot_samples_coordinates_pca(
    pca_res,
    datasets = dataset,
    pcs = 1:3,
    mo_data = mo_set,
    colour_upper = "status",
    scale_colour_upper = ggplot2::scale_colour_manual(values = status_colors),
    colour_lower = covariate,
    scale_colour_lower = ggplot2::scale_colour_brewer(palette = "Paired", na.value = "grey50")
  )

  ## Combining the two plots
  list(p_screeplot, GGally::ggmatrix_gtable(p_pcs)) |>
    purrr::map(patchwork::wrap_elements) |>
    patchwork::wrap_plots(
      ncols = 1,
      heights = c(1, 2)
    ) +
    patchwork::plot_annotation(tag_levels = "A", tag_suffix = ")") &
    ggplot2::theme(plot.tag = ggplot2::element_text(size = 16))
}

#' Samples score plot vs disease status
#'
#' Plots for each integration method the samples score for the first two latent
#' dimensions generated, colouring the samples according to their disease
#' status.
#'
#' @param list_results List of integration results from the different methods,
#'   each as an `output_dimension_reduction` object.
#' @param mo_set MultiDataSet, the multi-omics dataset.
#' @returns A patchwork of plots.
make_integration_samples_score_plot <- function(list_results, mo_set) {
  ## Extract name of the method from each result
  names(list_results) <- purrr::map(list_results, \(x) attr(x, "method"))

  ## For all methods except sO2PLS, extract the name of the first two latent
  ## dimensions from their results
  methods_except_so2pls <- setdiff(names(list_results), "sO2PLS")

  first_two_dimensions <- list_results[methods_except_so2pls] |>
    purrr::map(\(x) moiraine::get_latent_dimensions(x)[1:2])

  ## Creating the samples scatterplot for all methods with 2 latent dimensions
  plots_list <- methods_except_so2pls |>
    rlang::set_names() |>
    purrr::map(\(x) {
      moiraine::plot_samples_score_pair(
        list_results[[x]],
        latent_dimensions = first_two_dimensions[[x]],
        mo_data = mo_set,
        colour_by = "status"
      ) +
        ggplot2::scale_colour_manual(values = status_colors)
    })

  # For sO2PLS, create a boxplot instead (because it only has 1 latent dimension)
  plots_list[["sO2PLS"]] <- moiraine::plot_samples_score_covariate(
    list_results[["sO2PLS"]],
    mo_data = mo_set,
    covariate = "status",
    latent_dimensions = "joint component 1",
    colour_by = "status"
  ) +
    ggplot2::scale_colour_manual(values = status_colors, guide = "none") +
    ggplot2::scale_fill_manual(values = status_colors, guide = "none")

  ## Combining the plots
  patchwork::wrap_plots(
    plots_list[names(list_results)],
    nrow = 2,
    byrow = TRUE,
    guides = "collect"
  ) +
    patchwork::plot_annotation(tag_levels = "A", tag_suffix = ")") &
    ggplot2::theme(
      legend.position = "bottom",
      plot.tag = ggplot2::element_text(size = 16)
    )
}

#' Features importance score vs single-omics results
#'
#' Displays the top features from MOFA factor 1 as well as the distribution of
#' features importance score for factor 1 against the results from single-omics
#' analyses.
#'
#' @param mofa_output MOFA output in standardised format.
#' @param mo_set MultiDataSet object, the multi-omics dataset.
#' @returns A patchwork of plots.
make_mofa_features_plot <- function(mofa_output, mo_set) {

  p1 <- moiraine::plot_top_features(
    mofa_output,
    latent_dimensions = "Factor 1",
    n_features = 10,
    mo_data = mo_set,
    label_cols = list("rnaseq" = "Name", "metabolome" = "name"),
    truncate = 40
  )[[1]] +
    ggplot2::ggtitle("Top 10 features - MOFA Factor 1")

  p2 <- moiraine::plot_features_weight_covariate(
    mofa_output,
    mo_data = mo_set,
    covariate = list(
      "snps" = "qtl_type",
      "rnaseq" = "de_status",
      "metabolome" = "de_status"
    ),
    latent_dimensions = "Factor 1"
  ) +
    ggplot2::labs(
      x = "QTL analysis results (snps), DE results (rnaseq, metabolome)"
    )

  patchwork::wrap_plots(
    list(p1, patchwork::free(p2)),
    nrow = 2,
    byrow = FALSE,
    heights = c(1, 1.3)
  ) +
    patchwork::plot_annotation(tag_levels = "A", tag_suffix = ")") &
    ggplot2::theme(
      plot.tag = ggplot2::element_text(size = 16)
    )
}


#' Format GO terms enrichment results
#'
#' Formats the table of GO terms enrichment results.
#'
#' @param df Tibble, table of GO terms enrichment results.
#' @returns A tibble.
make_mofa_enrichment_table <- function(df) {
  df |>
    dplyr::filter(pvalue < 0.05) |>
    dplyr::select(
      `GO term` = set_id,
      `GO term name` = go_name,
      `Statistics` = stat_mean,
      `p-value` = pvalue,
      `Adjusted p-value` = adj_pvalue,
      `Number of genes in GO term` = set_size
    )
}



#' Comparison of integration methods
#'
#' Displays the matrix of samples score and features weight correlation between
#' the first three latent dimensions of each method, as well as a comparison of
#' the features importance between the first latent dimension of MOFA and
#' DIABLO.
#'
#' @param list_results List of integration results from the different methods,
#'   each as an `output_dimension_reduction` object.
#' @param mo_set MultiDataSet, the multi-omics dataset.
#' @returns A patchwork of plots.
make_methods_comparison_plot <- function(output_list, mo_set) {

  names(output_list) <- purrr::map(output_list, \(x) attr(x, "method"))

  ## Plotting comparison heatmap
  p1 <- moiraine::comparison_heatmap_corr(
    output_list,
    latent_dimensions = list(
      "sO2PLS" = c(
        "joint component 1",
        "metabolome specific component 1",
        "rnaseq specific component 1"
      ),
      "MOFA" = paste("Factor", 1:3),
      "DIABLO" = paste("Component", 1:3)
    )
  ) |>
    ComplexHeatmap::draw() |>
    grid::grid.grabExpr()

  # moiraine::plot_samples_score_pair(
  #   output_list[c("sPLS", "MOFA")],
  #   latent_dimensions = c("sPLS" = "Component 2", "MOFA" = "Factor 4"),
  #   mo_data = mo_set,
  #   colour_by = "gender"
  # )

  p2 <- moiraine::plot_features_weight_pair(
    output_list[c("MOFA", "DIABLO")],
    latent_dimensions = c("MOFA" = "Factor 1", "DIABLO" = "Component 1"),
    mo_data = mo_set,
    features_metric = "importance",
    label_cols = list(
      "rnaseq" = "Name",
      "metabolome" = "name"
    )
  )

  design <- c(
    "#AAAAAAAAA#
     #AAAAAAAAA#
     #AAAAAAAAA#
     BBBBBBBBBBB
     BBBBBBBBBBB"
  )

  patchwork::wrap_plots(
    list(p1, patchwork::free(p2))
    #nrow = 2,
    #byrow = FALSE,
    #heights = c(1.5, 1)
  ) +
    patchwork::plot_layout(design = design) +
    patchwork::plot_annotation(tag_levels = "A", tag_suffix = ")") &
    ggplot2::theme(
      plot.tag = ggplot2::element_text(size = 16)
    )
}

#' Comparison of MOFA results with different features pre-filtering approaches
#'
#' @param list_results List of integration results from the different MOFA runs,
#'   each as an `output_dimension_reduction` object.
#' @returns A plot from ComplexHeatmap.
make_prefiltering_comparison_plot <- function(output_list) {
  moiraine::comparison_heatmap_corr(
    output_list,
    latent_dimensions = list(
      "MOFA (supervised pref.)" = paste0("Factor ", 1:5),
      "MOFA (unsupervised pref.)" = paste0("Factor ", 1:5)
    ),
    legend_ncol = 1
  )
}

#' Generates a table containing the outcome of the supervised and unsupervised
#' features pre-filtering.
#'
#' @param mo_set MultiDataSet object, the non-filtered multi-omics dataset.
#' @param mo_set_sup_filt MultiDataSet object, the multi-omics dataset after
#'   supervised pre-filtering.
#' @param mo_set_unsup_filt MultiDataSet object, the multi-omics dataset after
#'   unsupervised pre-filtering.
#' @returns A tibble with one row per feature.
make_prefiltering_suppl_table <- function(mo_set,
                                          mo_set_sup_filt,
                                          mo_set_unsup_filt) {
  datasets <- c("snps", "rnaseq")

  sup_filt_features <- moiraine::get_features(mo_set_sup_filt)[datasets]
  unsup_filt_features <- moiraine::get_features(mo_set_unsup_filt)[datasets]

  moiraine::get_features(mo_set)[datasets] |>
    purrr::imap(\(x, ds) {
      tibble::tibble(
        feature_id = x,
        supervised_prefiltering = feature_id %in% sup_filt_features[[ds]],
        unsupervised_prefiltering = feature_id %in% unsup_filt_features[[ds]]
      ) |>
        dplyr::mutate(
          dplyr::across(
            .cols = -feature_id,
            .fns = ~ dplyr::case_when(
              .x ~ "retained",
              !.x ~ "discarded"
            )
          )
        )
    }) |>
    purrr::list_rbind(names_to = "dataset")
}

#' Generate a list of tables with features weight and samples score for the
#' different integration methods
#'
#' @param res_list Named list of dimension reduction output objects, the
#'   integration results.
#' @returns A named list of tibbles.
make_integration_suppl_tables <- function(res_list) {
  features_weight_list <- purrr::map(res_list, \(x) x$features_weight)
  names(features_weight_list) <- paste(names(res_list), "features weight")

  samples_score_list <- purrr::map(res_list, \(x) x$samples_score)
  names(samples_score_list) <- paste(names(res_list), "samples score")

  c(features_weight_list, samples_score_list)
}
