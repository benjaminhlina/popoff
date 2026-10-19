#' Message functions
#'
#' These functions produce information to the user.
#'
#' @param bathy a bathymetric feature.
#'
#' @name msg_function
#' @keywords internal

msg_bathy <- function(bathy) {
  if (inherits(bathy, "SpatRaster")) {
    cli::cli_inform(
      "{.arg bathy} detected as {.cls SpatRaster}. Using as-is."
    )
    return(invisible(bathy))
  }

  if (inherits(bathy, "RasterLayer")) {
    cli::cli_inform(
      "{.arg bathy} detected as {.cls RasterLayer}. Converting to {.cls SpatRaster}."
    )
    return(invisible(bathy))
  }

  if (is.matrix(bathy)) {
    cli::cli_inform(
      "{.arg bathy} detected as a {.cls matrix} ({nrow(bathy)}x{ncol(bathy)}). Converting to {.cls SpatRaster}."
    )
    return(invisible(bathy))
  }

  if (is.character(bathy)) {
    if (
      grepl("\\.RData$", bathy, ignore.case = TRUE) ||
        grepl("\\.rds$", bathy, ignore.case = TRUE) ||
        grepl("\\.tif$", bathy, ignore.case = TRUE)
    ) {
      cli::cli_inform(
        "{.arg bathy} detected as a file path: {.file {bathy}}"
      )
      return(invisible(bathy))
    }
  }
}


#' @param df a `data.frame`
#' @param n_days a `numerical` value
#'
#' @name msg_function
#' @keywords internal
msg_dates <- function(df, n_days) {
  d_min <- min(df$date)
  d_max <- max(df$date)
  d_span <- difftime(d_max, d_min)

  cli::cli_alert_info(paste0(
    "Generating depth profile likelihood for {.val {d_min}} through ",
    "{.val {d_max}}, a total of {.val {d_span}} days, while only ",
    "{.val {n_days}} days had observed depth values."
  ))
}
#' @param workers a numerical value
#' @name msg_function
#' @keywords internal
msg_daemon <- function(workers) {
  cli::cli_inform(
    "Daemons already running, ignoring {.arg workers} and using a total of
        {.val {mirai::status()$connections}} workers"
  )
}

#' @param t1 end time
#' @name msg_function
#' @keywords internal

msg_end <- function(chr, t0, t1 = round(Sys.time())) {
  cli::cli_alert_success(
    "{chr} took {round(as.numeric(difftime(t1, t0, units = 'mins')), 2)} minutes"
  )
  invisible(t1)
}

#' @param n_empty number of empty
#' @param out the out object.
#' @name msg_function
#' @keywords internal
msg_extract_no_val <- function(n_empty, out) {
  dates <- out$date[is.na(out$lon)][1]
  n_rows <- nrow(out)
  cli::cli_warn(
    "{n_empty} of {n_rows} date{?s} had no finite values; first at {dates}."
  )
}

#' @param n_fallback a numerical value
#' @name msg_function
#' @keywords internal
msg_fallback <- function(n_fallback) {
  if (n_fallback > 0L) {
    cli::cli_warn(
      "{n_fallback} step{?s} collapsed during smoothing and used the filtered posterior instead. Proceed interpretation with caution."
    )
  }
}
#' @param path a file path.
#' @name msg_function
#' @keywords internal

msg_load_raster <- function(path) {
  cli::cli_inform(
    "{.arg path} detected as a file path: {.file {path}}. Attempting to read with {.fn terra::rast}."
  )
}

#' @param obj an object
#' @param obj_name the name of an object
#' @name msg_function
#' @keywords internal
msg_load_rdata <- function(obj, obj_name) {
  cli::cli_inform(
    "Loaded object {.val {obj_name[1]}} of class {.cls {class(obj)}} from file."
  )
}
#' @param n_bad_days an integer of the bad esimtates
#' @name msg_function
#' @keywords internal
msg_no_mass <- function(n_bad_days) {
  if (n_bad_days > 0) {
    cli::cli_alert_warning(
      "The negative log likelihood has the follwoing nubmer of {.val {n_bad_days}}
      days that produced no estimate. Proceed with using this estimate and 
      value of sigma with caution."
    )
  }
}

#' @param progress a `logical` that controls whether the msg occurs.
#' @param n_time a `numeric` value that is the number of time intervals
#' @param name a `character` vector that contains the name of the progress bar
#' @param type a `character` vector that contains whether filtering or smoothing
#' @param .envir the envirnoment
#' @name msg_function
#' @keywords internal
msg_progress <- function(
  progress,
  n_time,
  name = c("HMM forward filtering", "HMM backwards smoothing"),
  type = c("Filtering", "Smoothing"),
  .envir = parent.frame()
) {
  if (isTRUE(progress)) {
    name <- match.arg(name)
    type <- match.arg(type)

    cli::cli_progress_bar(
      name = name,
      format = paste0(
        "{cli::pb_spin} ",
        type,
        " time step {cli::pb_current}/{cli::pb_total} ",
        "{cli::pb_bar} {cli::pb_percent} | ETA: {cli::pb_eta}"
      ),
      total = n_time - 1,
      clear = FALSE,
      .envir = .envir
    )
  }
}

#' @name msg_function
#' @keywords internal
msg_progress_update <- function(progress) {
  if (isTRUE(progress)) {
    cli::cli_progress_update()
  }
}

#' @name msg_function
#' @keywords internal
msg_saved_plot <- function(path) {
  cli::cli_alert_success(
    "Saved plot at the following location {.path {path}}"
  )
}
#' @param obj an object
#' @param date A date value
#' @param linear_index Position of the known location within a daily slice
#' @name msg_function
#' @keywords internal
msg_slice <- function(obj, date, linear_index) {
  if (is.na(obj[linear_index])) {
    cli::cli_warn(
      "Known location on {.val {date}} falls on an NA (masked/land) cell."
    )
  }
}


#' @param chr a `character` string that tells which likelihood such as e.g., `"Depth"`
#' @param t0 current time
#' @name msg_function
#' @keywords internal

msg_start <- function(chr, t0 = round(Sys.time())) {
  cli::cli_alert_info("Starting {chr} calculations at {.val {t0}}")
  invisible(t0)
}

#' @param window a numerical value
#' @name msg_function
#' @keywords internal

msg_window <- function(window) {
  if (window < 1) {
    cli::cli_alert_info(
      "Time window, which is the `bin` (e.g., 3600 sec = 1 h) divided by the PSAT tranmission time (e.g., 360) is < {.val {1}} and thus has been adjusted to {.val {1}}."
    )
  }
  if (window >= 1) {
    cli::cli_alert_info(
      "Time window, which is the `bin` (e.g., 3600 sec = 1 h) divided by the PSAT tranmission time (e.g., 360), is {.val {window}}."
    )
  }
}
