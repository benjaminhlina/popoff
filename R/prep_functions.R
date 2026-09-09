#' Prepare Bathymetic Data
#'
#' This function prepared bahtymetic data to be used in the likelihood estimates.
#'
#' @param bathy a `spatrater` generated using `load_bathy_raster(())`
#' @param weighted a `logical` value that determines whther or not sd is weighted. Default is `FALSE`
#' @param focal_dim a `numerical` value that is odd. See details below.
#' @param ctd_dir only required when `weighted = TRUE`. The path to ctd.
#'
#' @details  window used (i.e., `focal_dim`) must have odd dimensions. If you need even sides,
#' you can use a matrix and add a column or row of `NA`'s to mask out values.
#' Window values are typically `1` or `NA` to indicate whether a value is used or ignored in computations, respectively.
#' `NA`` values in w can be useful for creating non-rectangular (e.g. circular) windows.
#' A weights matrix of numeric values can also be supplied to `focal_dim`.
#' In the case of a weights matrix, cells with `NA`` weights will be ignored, and the rest of the values in the
#' focal window will be multiplied by the corresponding weight prior to `sd` being applied.
#'  Note, `na.rm` does not need to be TRUE if `focal_dim` contains `NA`` values as these cells are ignored in computations.
#'
#' @seealso [terra::focal()] [load_bathy_raster()]
#'
#' @export

prep_bathy_likdepth <- function(
  bathy, # SpatRaster, RasterLayer, matrix, or path to .RData/.tif
  weighted = FALSE,
  focal_dim = 9,
  ctd_dir = NULL # only required when weighted = TRUE
) {
  rlang::check_required(bathy)

  error_spatrast(bathy)
  error_logical(weighted)
  error_focal_dim(focal_dim)
  t0 <- msg_start("bathymetric standard deviation")

  idx <- attr(bathy, 'idx')
  idx_na <- attr(bathy, 'idx_na')

  # ---- local bathymetric roughness (sd within a moving window) ----
  if (weighted) {
    if (is.null(ctd.dir)) {
      stop("ctd.dir must be supplied when weighted = TRUE")
    }
    if (!focal_dim %in% c(3, 5, 9)) {
      stop("pre-computed sd matrices only exist for focal_dim = 3, 5, or 9")
    }
    sd_path <- file.path(
      ctd.dir,
      paste0("sdcell.n", focal_dim, ".weighted.RData")
    )
    sdi <- local({
      e <- new.env()
      load(sd_path, envir = e)
      e$sdcell
    })
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

    # set crs and georeferecne
    terra::crs(sdi) <- terra::crs(bathy)
  }
  msg_end(chr = "Bathymetric standard deviations", t0 = t0)

  list(
    bathy = bathy,
    sdi = sdi, # local bathymetric sd matrix
    idx = idx,
    idx_na = idx_na,
    focal_dim = focal_dim,
    weighted = weighted
  )
}

# ----- prep_summmary_likdepth -----

#' Prepare Summary Data
#'
#' This function prepared summary data to be used in the likelihood estimates.
#'
#' @param raw_log a `data.frame` containing the raw pop-off satalitte tag log
#' @param iniloc a `data.frame` containing the iniloc information about when the tag was deployed, mrpat, and pop dates and location.
#' @param bin a `numeric` containing the number of seconds within an aggregation (e.g., 3600 seconds is 1 hour).
#' @param make_plot a `logical` value that controls whether the to make the summary plot. Default is `TRUE`.
#' @param output_dir Defaults to `NULL` and will not export plot. Please provide a `character` string containing the path deseired
#' to save the plot.
#'
#' @export

prep_summary_likdepth <- function(
  raw_log,
  iniloc,
  bin,
  output_dir = NULL,
  make_plot = TRUE
) {
  # ----- first check if arguments needed are supplied ----
  rlang::check_required(raw_log)
  rlang::check_required(iniloc)
  rlang::check_required(bin)

  error_input_data(input_data = raw_log, type = "raw_log")
  error_input_data(input_data = iniloc, type = "iniloc")
  error_out_dir(output_dir)
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
      output_dir = output_dir
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
