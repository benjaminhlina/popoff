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
  cli::cli_alert_info(
    "Generating depth profile likelihood for {.val {min(df$date)}} through {.val {max(df$date)}} which is a total of {n_days} days."
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
