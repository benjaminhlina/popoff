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
extract_locations <- function(
  obj,
  type = c("max", "min", "mean", "median"),
  drop_empty = FALSE
) {
  type <- match.arg(type)

  dn <- dimnames(obj)
  lats <- as.numeric(dn$lat)
  lons <- as.numeric(dn$lon)
  dates <- as.Date(dn$date)
  n_row <- dim(obj)[1]

  result <- lapply(seq_along(dates), function(i) {
    val <- obj[,, i]
    val[!is.finite(val)] <- NA

    if (all(is.na(val))) {
      return(data.frame(
        date = dates[i],
        lon = NA,
        lat = NA,
        value = NA
      ))
    }

    valid <- which(!is.na(val))
    v <- val[valid]

    idx <- switch(
      type,
      max = valid[which.max(v)],
      min = valid[which.min(v)],
      mean = valid[which.min(abs(v - mean(v)))],
      median = valid[which.min(abs(v - stats::median(v)))]
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

  n_empty <- sum(is.na(out$lon))
  if (n_empty > 0) {
    msg_extract_no_val(n_empty = n_empty, out = out)
  }
  if (drop_empty) {
    out <- out[!is.na(out$lon), ]
  }
  out
}
