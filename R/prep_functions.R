#' Prepare Bathymetic Data
#'
#' This function prepared bahtymetic data to be used in the likelihood estimates.
#' It calculates the standard deviation of each depth value by using the values
#' around a given depth. To control the area used for standard deviation use
#' `focal_dim`.
#'
#' @param bathy a `SpatRaster` generated using `load_raster(())`
#' @param focal_dim a `numerical` value that is odd. See details below.
#' @param ctd a `SpatRaster`that contains conductivity, temperature and depth
#' (ctd) values for a given area. It will be used as informed standard deviation
#' of the bathymetry instead of calculating the stadard deviation.
#'
#' @details the window used (i.e., `focal_dim`) must have odd dimensions.
#' If you need even sides, you can use a matrix and add a column or row of `NA`'s
#' to mask out values. Window values are typically `1` or `NA` to indicate
#' whether a value is used or ignored in computations, respectively.
#' `NA`` values in w can be useful for creating non-rectangular (e.g. circular)
#' windows. A weights matrix of numeric values can also be supplied to `focal_dim`.
#' In the case of a weights matrix, cells with `NA`` weights will be ignored,
#' and the rest of the values in the focal window will be multiplied by the
#' corresponding weight prior to `sd` being applied.
#' Note, `na.rm` does not need to be TRUE if `focal_dim` contains `NA`` values as these cells are ignored in computations.
#'
#' @return A `list()` containing a `SpatRaster` of both bathymetry and the
#' standard deviation (`sdi`), the index values (`idx`) for observed depths and
#' `NA` values (`idx_na`), and the focal dimensions used to calculate the
#' standard deviation.
#' @seealso [terra::focal()] [load_raster()]
#'
#' @export

prep_bathy_likdepth <- function(
  bathy,
  focal_dim = 9,
  ctd = NULL
) {
  rlang::check_required(bathy)

  error_spatrast(bathy)
  error_focal_dim(focal_dim)
  t0 <- msg_start("bathymetric standard deviation")

  idx <- attr(bathy, 'idx')
  idx_na <- attr(bathy, 'idx_na')

  # ---- local bathymetric roughness (sd within a moving window) ----
  if (!is.null(ctd)) {
    error_spatrast(ctd)
    sdi <- ctd
    focal_dim <- NULL
  } else {
    sd_r <- terra::focal(
      bathy,
      w = focal_dim,
      fun = "sd",
      na.rm = TRUE
    )
    sdi <- terra::as.matrix(sd_r, wide = TRUE)

    ssd <- as.vector(sdi)
    fill_sd0 <- stats::quantile(ssd[ssd != 0], probs = 0.01, na.rm = TRUE)

    if (is.na(fill_sd0)) {
      fill_sd0 <- 1e-3
    }
    sdi[sdi == 0] <- fill_sd0
    sdi[is.na(sdi)] <- fill_sd0
    sdi <- terra::rast(sdi)
  }

  # set crs and georeferecne
  terra::crs(sdi) <- terra::crs(bathy)
  msg_end(chr = "Bathymetric standard deviations", t0 = t0)

  out <- list(
    bathy = bathy,
    sdi = sdi,
    idx = idx,
    idx_na = idx_na,
    focal_dim = focal_dim
  )
  return(out)
}

# ----- prep_summmary_likdepth -----

#' Prepare Summary Data
#'
#' This function prepared summary data to be used in the likelihood estimates.
#'
#' @param raw_log a `data.frame` containing the raw pop-off satalitte tag log
#' @param iniloc a `data.frame` containing the iniloc information about when the
#' tag was deployed, mrpat, and pop dates and location.
#' @param bin a `numeric` containing the number of seconds within an
#' aggregation (e.g., 3600 seconds is 1 hour).
#' @param make_plot a `logical` value that controls whether the to make the
#' summary plot. Default is `TRUE`.
#' @param out_dir Defaults to `NULL` and will not export plot.
#' Please provide a `character` string containing the path desired
#' to save the plot.
#'
#' @return A `list()` containing a the `raw_log` as `data.frame`, the summarized
#' log as a `data.frame`, a `ggplot2` object (`depth_profile`), the time
#' resolution, and the time bin used as single value `vectors`.
#'
#' @export

prep_summary_likdepth <- function(
  raw_log,
  iniloc,
  bin = NULL,
  out_dir = NULL,
  make_plot = TRUE
) {
  # ----- first check if arguments needed are supplied ----
  rlang::check_required(raw_log)
  rlang::check_required(iniloc)
  rlang::check_required(bin)

  error_input_data(input_data = raw_log, type = "raw_log")
  error_input_data(input_data = iniloc, type = "iniloc")
  error_out_dir(out_dir)
  error_logical(make_plot)

  error_numeric(bin)

  iniloc_f <- iniloc |>
    dplyr::filter(fish_id %in% raw_log$fish_id)

  tag_date <- iniloc_f$date[iniloc_f$event == "tag"]

  pop_date <- iniloc_f$date[iniloc_f$event == "pop"]

  raw_log_filter <- raw_log |>
    dplyr::filter(date >= as.Date(tag_date), date <= as.Date(pop_date)) |>
    dplyr::mutate(
      depth = depth * -1
    )

  if (isTRUE(make_plot)) {
    depth_profile <- make_summary_plot(
      raw_log = raw_log_filter,
      out_dir = out_dir
    )
  } else {
    depth_profile <- NULL
  }

  # ---- make median lag -----
  median_min_lag <- raw_log_filter |>
    dplyr::arrange(date_time) |>
    dplyr::mutate(
      min_lag = as.numeric(difftime(
        date_time,
        dplyr::lag(date_time),
        units = "secs"
      ))
    ) |>
    dplyr::pull(min_lag) |>
    stats::median(na.rm = TRUE)

  # make summary_table

  summary_table <- make_summary_table(
    raw_log = raw_log_filter,
    time_res = median_min_lag,
    bin = bin
  )

  summary_list <- list(
    raw_log_filter = raw_log_filter,
    summary_table = summary_table,
    depth_profile = depth_profile,
    time_res = median_min_lag,
    bin = bin
  )
  return(summary_list)
}
