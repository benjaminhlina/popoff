#' Make functions
#'
#' These functions make different objects that can be exported like plots and summary tables.
#'
#' @param raw_log a `data.frame` containing raw psat logs for an individual.
#' @param out_dir a `character` string containing the path to save plot or table.
#'
#' @name make_function
#' @export

make_summary_plot <- function(raw_log, out_dir = NULL) {
  pdt_long <- raw_log |>
    dplyr::select(date_time, depth, temperature) |>
    tidyr::pivot_longer(
      cols = c(depth, temperature),
      names_to = "variable",
      values_to = "value"
    ) |>
    # Source - https://stackoverflow.com/a/32758968
    # Posted by irJvV
    # Retrieved 2026-08-29, License - CC BY-SA 3.0
    dplyr::mutate(
      variable = paste0(
        toupper(substr(variable, 1, 1)),
        substr(variable, 2, nchar(variable))
      )
    )

  fish_id <- unique(raw_log$fish_id)

  depth_profile <- ggplot2::ggplot(
    pdt_long,
    ggplot2::aes(x = date_time, y = value, colour = variable)
  ) +
    ggplot2::facet_wrap(~variable, scales = "free", ncol = 1) +
    ggplot2::geom_point(size = 1) +
    ggplot2::scale_x_datetime(
      labels = scales::date_format("%m/%y"),
      date_breaks = "1 month"
    ) +
    ggplot2::theme_bw(base_size = 15) +
    ggplot2::theme(
      axis.title.y = ggplot2::element_blank(),
      legend.position = "none",
      strip.background = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5)
    ) +
    ggplot2::labs(x = "Date (mm/yy)", title = fish_id)

  if (!is.null(out_dir)) {
    path <- file.path(
      out_dir,
      paste0("depth_profile_", fish_id, ".png")
    )
    ggplot2::ggsave(
      filename = path,
      plot = depth_profile,
      width = 11,
      height = 8.5,
    )

    msg_saved_plot(path = path)
  }
  invisible(depth_profile)
}

#' @param raw_log a `data.frame` containing raw psat logs for an individual.
#' @param time_res a `numeric` that is the median time delay between depth logs, in seconds
#' @param bin a `numeric` containing the number of seconds within an aggregation (e.g., 3600 seconds is 1 hour)
#'
#' @name make_function
#' @export

make_summary_table <- function(raw_log, time_res, bin) {
  t0 <- msg_start("summary")

  time_res_chr <- paste(time_res, "sec")

  window_time <- round(bin / time_res)
  msg_window(window_time)

  if (window_time < 1) {
    window_time <- 1
  }

  # ---- first we need to create time resoluation and sequence
  # for eavery day end....
  raw_log_res <- raw_log |>
    dplyr::group_by(date) |>
    dplyr::mutate(
      day_end = as.POSIXct(
        paste(
          date,
          "23:59:59"
        ),
        tz = "UTC"
      )
    ) |>
    dplyr::arrange(date_time)

  # ----- next we need to complete (expand.grid) on a sequency that takes the min date an time
  # for a given day
  raw_log_fill <- raw_log_res |>
    tidyr::complete(
      date_time = seq(
        from = min(date_time, na.rm = TRUE),
        to = day_end[1],
        by = time_res
      )
    )

  raw_log_depth <- raw_log_fill |>
    # ---- compute rolling metrics across regular time grid ----
    dplyr::mutate(
      depth_run = slider::slide_dbl(
        depth,
        .safe_max,
        .before = window_time[1] - 1,
        .after = 0,
        .complete = FALSE
      ),
      depth_run = dplyr::if_else(
        !is.na(depth),
        true = depth_run,
        false = NA
      )
    )

  summary_log <- raw_log_depth |>
    # ---- Summarize depth & temperature per date ------
    dplyr::summarize(
      depth_range = list(.compute_depth_range(depth, depth_run)),
      depth_min = depth_range[[1]][1],
      depth_max = depth_range[[1]][2],

      valid_temps = list(temperature[!is.na(temperature)]),
      temp_subset = list(.compute_temp_subset(
        depth,
        temperature,
        depth_min,
        depth_max
      )),

      temp_min = .compute_temp_bound(temp_subset[[1]], valid_temps[[1]], "min"),
      temp_max = .compute_temp_bound(temp_subset[[1]], valid_temps[[1]], "max"),
      .groups = "drop"
    ) |>
    dplyr::select(-valid_temps, -temp_subset, -depth_range)

  msg_end(chr = "Daily summary", t0 = t0)
  return(summary_log)
}
