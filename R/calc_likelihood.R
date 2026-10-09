#' Calculate Likelihood Functions
#'
#' These are functions that will calculated the likelihood
#' of movement based on additonal parameters such as depth
#'
#' @param summary_log_i `data.frame` of a single date's rows (as produced by splitting `summary_log` on `date` upstream); may be 0-row or have NA/Inf `depth_max`
#' @param zbi a `vector` containing bathymetric values for an area `idx`, `idx_na`, `bathy`, `sdi`
#' @param sdi_dsd a `vector` containing standard deviaation values for an area + the daily_sd_depth
#' @param idx a `vector` containing the index of valid depth values
#' @param idx_na a `vector` containing the index of non-valid depth value/where depth is `NA`
#' @param dims integer vector length 2, c(nrow, ncol) of the original bathymetry grid
#' @param time character/Date, used only for error messages
#' @keywords internal
calc_daily_likdepth <- function(
  summary_log_i,
  zbi,
  sdi_dsd,
  idx,
  idx_na,
  dims,
  time
) {
  # no usable depth range that day -> uninformative likelihood
  if (
    nrow(summary_log_i) == 0 ||
      is.na(summary_log_i$depth_max[1]) ||
      is.infinite(summary_log_i$depth_max[1])
  ) {
    lik <- zbi
    lik[] <- 0
    lik[idx] <- 1
    return(matrix(lik, nrow = dims[1], ncol = dims[2]))
  }

  pulled_range <- summary_log_i[1, c("depth_min", "depth_max")] |>
    as.numeric() |>
    sort()

  # inflate observed depth range by tag error: +-1 m and +-1%
  low <- pulled_range[1] - (1 + 0.01 * abs(pulled_range[1]))
  high <- pulled_range[2] + (1 + 0.01 * abs(pulled_range[2]))

  lik <- tryCatch(
    likint3_matrix(
      w = cbind(
        zbi,
        sdi_dsd,
        low,
        high
      ),
      widx = idx
    ),
    error = function(e) {
      cli::cli_abort(
        c(
          "x" = "Depth-transmission likelihood failed for {.val {as.character(time)}}.",
          "i" = "Depth_min/Depth_max: {.val {pulled_range[1]}} / {.val {pulled_range[2]}}"
        ),
        parent = e
      )
    }
  )

  lik_mat <- matrix(lik, nrow = dims[1], ncol = dims[2])
  lik_mat[is.na(lik_mat)] <- 0
  lik_mat[idx_na] <- NA

  max_val <- max(lik_mat, na.rm = TRUE)
  if (is.finite(max_val) && max_val > 0) {
    lik_mat <- lik_mat / max_val
  }
  return(lik_mat)
}

#' Calculate Likelihood Depth Transmitted
#'
#' These are functions that will calculated the likelihood of depth occupancy based on
#' depth values transmitted and bathymetry
#'
#' @param prepped_log a `list` containing summary objects crated by `prep_summary_likdepth()`.
#' @param prepped_bathy a `list` containing objects created by `prep_bathy_likdepth()`.
#' @param event_log a `data.frame` containing the event log data such as
#' tagging, mrpat, and and pop-off date and time.
#' popoff data.
#' @param daily_sd_depth a numeric value
#'
#' @export

calc_likdepth_transmit <- function(
  prepped_log,
  prepped_bathy,
  event_log,
  daily_sd_depth = 0.57
) {
  # ctd.dir="C:/Users/pgatti/Documents/halibut/R/Data_PGalbraith/toBeLoaded/";
  # ncores=4;weighted=FALSE; focalDim = 3;daily.sd.depth=.57;
  # pdt=sm$sm;depth.list=depth.mix

  # rlang::arg_match(prepped_log)
  # rlang::arg_match(prepped_bathy)

  error_list_prep(prepped_log, type = "log")
  error_list_prep(prepped_bathy, type = "bathy")
  error_numeric(daily_sd_depth)
  error_input_data(input_data = event_log, type = "event_log")

  # ---- msg_start
  t0 <- msg_start(chr = "depth (m) likelihood")
  # grab summary table
  summary_table <- prepped_log$summary_table

  # make sure its ordered by date and then split it
  sum_split <- summary_table |>
    dplyr::arrange(date) |>
    (\(.) split(., .$date))()

  # get the totla number o fdays
  n_days <- length(sum_split)

  msg_dates(df = summary_table, n_days = n_days)

  # grab indx and index na
  idx <- prepped_bathy$idx
  idx_na <- prepped_bathy$idx_na

  # grab bathymetry and convert it to a vector

  zbi <- .convert_spr_to_vec(x = prepped_bathy$bathy)
  zbi_dims <- dim(prepped_bathy$bathy)

  # ----- create index_vector ----
  # this might be same as idx and idx_na
  idx_vec <- .convert_spr_to_vec(
    x = terra::init(prepped_bathy$bathy, fun = "cell")
  )
  # force into matrix on the dims of zbi
  idx_mat <- array(idx_vec, dim = zbi_dims[1:2])

  # grab sdi and convert to vector
  sdi <- .convert_spr_to_vec(x = prepped_bathy$sdi)

  # add daily_sd_depth to sdi
  sdi_dsd <- sdi + daily_sd_depth

  # ----- modify the event_log table -----
  event_log <- event_log |>
    dplyr::mutate(
      date = lubridate::floor_date(date, unit = '1 day') |>
        as.character(),
      raster_cell = terra::cellFromXY(prepped_bathy$bathy, cbind(lon, lat)),
      linear_index = match(raster_cell, idx_mat)
    )

  #---------------------------------------------------------------------
  # start of heavy daily loop
  # here is the place to implement parallelisation
  #--------------------------
  # ---- get date names -----

  dates <- names(sum_split)

  error_linear_index(
    event_date = event_log$date,
    linear_index = event_log$linear_index
  )

  error_missing_dates(event_log$date, dates)

  lik_list <- purrr::imap(
    sum_split,
    function(summary_log_i, date_i) {
      calc_daily_likdepth(
        summary_log_i = summary_log_i,
        zbi = zbi,
        sdi_dsd = sdi_dsd,
        idx = idx,
        idx_na = idx_na,
        dims = zbi_dims,
        time = date_i
      )
    },
    .progress = "Calculating depth likelihood"
  )

  likelihood_array <- simplify2array(lik_list)
  dimnames(likelihood_array)[[3]] <- dates
  valid_mask <- .valid_mask(bathy_vec = zbi, dims = zbi_dims)

  purrr::pwalk(
    event_log,
    function(date, linear_index, ...) {
      # slice <- likelihood_array[,, date]

      # msg_slice(slice, date, linear_index)
      msg_slice(valid_mask, date, linear_index)

      slice <- array(NA, dim = dim(valid_mask))
      slice[valid_mask] <- 0
      slice[linear_index] <- 1
      likelihood_array[,, date] <<- slice
    }
  )

  msg_end(chr = "Depth-transmission likelihood calculation", t0 = t0)

  return(likelihood_array)
}
