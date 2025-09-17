#' Write supplementary tables as sheets in Excel spreadsheet
#'
#' @param df_list Named list of tibbles, names are considered as title of the
#'   suppl. table.
#' @param output_file Character, path and name of file to create (including
#'   extension '.xlsx').
write_suppl_tables <- function(df_list, output_file) {
  wb <- openxlsx::createWorkbook()
  
  meta_df <- tibble::tibble(
    name = names(df_list)
  ) |>
    dplyr::mutate(
      suppl_table = paste0("Table S", 1:dplyr::n()),
      suppl_table = paste0(suppl_table, ": ", name),
      description = ""
    ) |>
    dplyr::select(
      `Suppl. Table` = suppl_table,
      Description = description
    )
  
  #openxlsx::renameWorksheet(wb, "Sheet1", "Information")
  openxlsx::addWorksheet(wb, "Information")
  openxlsx::writeData(wb, "Information", meta_df)
  
  names(df_list) <- meta_df$`Suppl. Table`
  
  worksheet_names <- meta_df$`Suppl. Table` |>
    stringr::str_remove("^Table ") |>
    stringr::str_replace(": ", "-") |>
    stringr::str_replace("unsupervised pref.", "u.p.") |>
    stringr::str_trunc(width = 31) |>
    rlang::set_names(meta_df$`Suppl. Table`)
  
  purrr::iwalk(
    df_list,
    \(x, y) {
      ws_name <- worksheet_names[[y]]
      
      openxlsx::addWorksheet(wb, ws_name)
      openxlsx::writeData(wb, ws_name, y) ## writing table title
      openxlsx::writeData(wb, ws_name, x, startRow = 2) ## writing data
    }
  )
  
  openxlsx::saveWorkbook(wb, output_file, overwrite = TRUE)
}

#' Plot memory usage of targets pipeline.
#'
#' @param logfile Character, path to the log file written by the `autometric`
#'   package.
#' @param tarmetafile Character, path to the csv file containing the
#'   `targets::tar_meta()` output.
#' @returns A ggplot.
make_autoplot_figure <- function(logfile, tarmetafile) {
  df_autometric <- autometric::log_read(logfile, units_memory = "gigabytes") |> 
    dplyr::mutate(time = time / 60)

  df_targets <- readr::read_csv(tarmetafile, show_col_types = FALSE) |> 
    dplyr::filter(!(type %in% c("function", "object"))) |> 
    dplyr::select(name, type, time, bytes, seconds) |> 
    dplyr::mutate(
      type = dplyr::case_when(
        stringr::str_detect(name, "individual_splsda_perf_") ~ "stem",
        TRUE ~ type
      ),
      name = dplyr::case_match(
        name,
        "individual_splsda_perf_c34ab0eca3a30184" ~ "individual_splsda_perf_snps",
        "individual_splsda_perf_782547ad94da4194" ~ "individual_splsda_perf_rnaseq", 
        .default = name
      )
    ) |> 
    dplyr::filter(!(name %in% c("individual_splsda_perf")))
  
  # To figure out which branch corresponds to which dataset
  # df_targets |> 
  #   dplyr::filter(type == "branch") |>
  #   dplyr::pull(name) |> 
  #   rlang::set_names() |> 
  #   purrr::map(\(x) attr(tar_read_raw(x), "dataset_name"))
  
  df_branch <- df_targets |> 
    dplyr::filter(type == "branch") |> 
    dplyr::mutate(name = stringr::str_remove(name, "_[^_]+$")) |> 
    dplyr::group_by(name) |> 
    dplyr::slice_min(time, n = 1, with_ties = FALSE) |> 
    dplyr::ungroup() |> 
    dplyr::select(name, pattern_time = time)
  
  df_targets <- df_targets |> 
    dplyr::left_join(df_branch, by = "name") |> 
    dplyr::mutate(time = dplyr::coalesce(time, pattern_time)) |> 
    dplyr::select(-pattern_time) |> 
    dplyr::filter(type %in% c("stem", "pattern")) |> 
    dplyr::arrange(time) |> 
    dplyr::mutate(
      target_order = rank(time),
      cum_time = cumsum(dplyr::lag(seconds, default = 0)),
      cum_time = cum_time / 60
    )
  
  df_labels <- df_targets |> 
    dplyr::filter(seconds > 300) |> 
    dplyr::mutate(
      n = 1:dplyr::n(),
      resident = (max(df_autometric$resident) / 2) + 0.3 * n
    )
  
  df_autometric |> 
    ggplot2::ggplot(ggplot2::aes(x = time, y = resident)) +
    ggplot2::geom_vline(
      ggplot2::aes(xintercept = cum_time), 
      data = df_targets,
      colour = "purple3",
      alpha = 0.6
    ) +
    ggplot2::geom_line(colour = "grey90") +
    ggplot2::geom_point(alpha = 0.8) +
    ggplot2::geom_label(
      ggplot2::aes(x = cum_time, y = resident, label = name),
      colour = "purple3",
      hjust = 0, 
      na.rm = TRUE,
      data = df_labels
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq(0, max(df_autometric$time), by = 30)
    ) +
    ggplot2::labs(
      x = "Elapsed time (min)",
      y = "RSS (Gb)"
    ) +
    ggplot2::theme_bw()
}

#' Plot running time for reduced vs full genomics dataset.
#'
#' @param logfile Character, path to the log file written by the `autometric`
#'   package.
#' @param tarmetafile Character, path to the csv file containing the
#'   `targets::tar_meta()` output.
#' @returns A ggplot.
make_figure_runtime_genomics <- function(logfile, tarmetafile) {
  readr::read_csv(tarmetafile, show_col_types = FALSE) |> 
    dplyr::filter(stringr::str_detect(name, "((pca_runs_list)|(individual_splsda_perf))$")) |> 
    dplyr::mutate(
      task = stringr::str_extract(name, "(pca)|(splsda)"),
      task = factor(task, levels = c("pca", "splsda"), labels = c("PCA", "sPLSDA")),
      dataset = stringr::str_extract(name, "(full)|(reduced)"),
      dataset = factor(
        dataset, 
        levels = c("reduced", "full"),
        labels = c("23,036 (reduced)", "85,630 (full)")
      ),
      time = seconds / 60,
      label = paste0(round(time, 1), " min"),
      label = dplyr::case_when(
        time == max(time) ~ paste0("\n", label),
        TRUE ~ paste0(label, "\n")
      ),
      label_col = dplyr::case_when(
        time == max(time) ~ "w",
        TRUE ~ "b"
      )
    ) |> 
    ggplot2::ggplot(ggplot2::aes(x = task, y = time, fill = dataset)) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::geom_text(
      ggplot2::aes(label = label, colour = label_col),
      position = ggplot2::position_dodge(width = 0.9),
      show.legend = FALSE
    ) +
    ggplot2::scale_colour_manual(values = c("b" = "black", "w" = "white")) +
    ggplot2::scale_fill_brewer(palette = "Paired") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.05))) +
    ggplot2::labs(
      x = NULL,
      y = "Running time (min)",
      fill = "Number of SNPs"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor.y = ggplot2::element_blank()
    )
}
