#' Download Functions
#'
#' These functions download different data (e.g., bathymetry) from different
#' sources.
#'
#' @param xlim a `vector` that contains two `numerical` values
#' which are the longitudinal extent.
#' @param ylim a `vector` that contains two `numerical` values
#' which are the latitudinal extent.
#' @param land_value a `numerical` value at which to exclude depth values.
#' Defaults to `0`.
#' @param ... additional arguments to pass to `getNOAA.bathy()`
#'
#' @seealso [marmap::getNOAA.bathy()]
#' @return a `matrix` containing a bathymetric values for a given extent.
#'
#' @export

download_noaa_bathy <- function(
  xlim,
  ylim,
  land_value = NULL,
  ...
) {
  suppressMessages(loadNamespace("marmap"))

  error_limits(xlim)
  error_limits(ylim)
  limits <- list(xlim = xlim, ylim = ylim)

  if (is.null(land_value)) {
    land_value <- 0
  }

  error_numeric(land_value)

  bath_matrix <- marmap::getNOAA.bathy(
    lon1 = min(xlim),
    lon2 = max(xlim),
    lat1 = min(ylim),
    lat2 = max(ylim),
    ...
  ) |>
    .cli_messages() |>
    unclass() |>
    t()

  lon <- colnames(bath_matrix) |>
    as.numeric() |>
    order()
  lat <- rownames(bath_matrix) |>
    as.numeric() |>
    order(decreasing = TRUE)

  half_x <- abs(lon[2] - lon[1]) / 2
  half_y <- abs(lat[2] - lat[1]) / 2

  bath_matrix <- bath_matrix[lat, lon, drop = FALSE]

  bath_matrix[bath_matrix > land_value] <- NA

  attr(bath_matrix, "xlim") <- c(min(xlim) - half_x, max(xlim) + half_x)
  attr(bath_matrix, "ylim") <- c(min(ylim) - half_y, max(ylim) + half_y)

  return(bath_matrix)
}
