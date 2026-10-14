#' Extract estimated locations from a normalised posterior array
#'
#' For each time step (day) that contains data, extracts a single estimated
#' location from the posterior surface using the requested statistic.
#'
#' @param obj A `list` containing an element named `normalised_posterior`, a
#' three-dimensional `array` with dimensions `[lat, lon, date]`. The
#' `dimnames` must be set: `lat` and `lon` as numeric-coercible coordinates
#' (e.g., projected metres) and `date` as `Date`-coercible strings.
#' @param type A `character` string giving the statistic used to select the
#' estimated cell for each day. One of `"max"` (default), `"min"`, `"mean"`,
#' or `"median"`.
#'
#' @details
#' Days where every cell is `NA` are dropped from the output.
#'
#' `"max"` and `"min"` return the cell with the highest and lowest posterior
#' value, respectively. A surface-wide mean or median has no location of its
#' own, so `"mean"` and `"median"` return the cell whose value is closest to
#' that statistic. Ties are resolved by taking the first matching cell.
#'
#' @return A `data.frame` with one row per day that has data and the columns:
#' \describe{
#'   \item{date}{`date` of the time step.}
#'   \item{lon}{Longitude (x) of the selected cell, in the units of the array
#'   `dimnames`.}
#'   \item{lat}{Latitude (y) of the selected cell, in the units of the array
#'   `dimnames`.}
#'   \item{value}{Posterior value of the selected cell.}
#' }
#'
#'
#' @examples
#' \dontrun{
#' extract_locations(final, type = "max")
#' extract_locations(final, type = "median")
#' }
#'
#' @export
extract_locations <- function(obj, type = c("max", "min", "mean", "median")) {
  type <- match.arg(type)

  lats <- as.numeric(dimnames(obj)$lat)
  lons <- as.numeric(dimnames(obj)$lon)
  dates <- as.Date(dimnames(obj)$date)

  n_row <- dim(obj)[1]

  has_data <- which(!apply(obj, 3, function(x) all(is.na(x))))

  result <- lapply(has_data, function(i) {
    val <- obj[,, i]

    # index of the cell that matches the requested statistic
    idx <- switch(
      type,
      max = which.max(val),
      min = which.min(val),
      mean = which.min(abs(val - mean(val, na.rm = TRUE))),
      median = which.min(abs(val - stats::median(val, na.rm = TRUE)))
    )

    row <- (idx - 1) %% n_row + 1
    col <- (idx - 1) %/% n_row + 1

    data.frame(
      date = dates[i],
      lon = lons[col],
      lat = lats[row],
      value = val[idx]
    )
  })

  out <- do.call(rbind, result)
  rownames(out) <- NULL

  return(out)
}
