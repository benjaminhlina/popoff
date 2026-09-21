#' Mask predicted/likelihood surfaces to a threshold-defined bounding region
#'
#'
#' @param predicted a `matrix``. Predicted spatial distribution at time t.
#' @param likelihood a `matrix``. Observation likelihood at time t.
#' @param lon a `numeric`` vector. Grid longitude values (one per column).
#' @param lat a `numeric`` vector. Grid latitude values (one per row).
#' @param par0 a `numeric`. Kernel diameter (odd), used to derive a movement
#'   radius so the bounding box isn't clipped tighter than the diffusion
#'   kernel would allow in one time step.
#' @param bound_threshold a `numeric`.
#' Fraction of likelihood's max used to decide which cells count
#' as a plausible  location (e.g. 0.1 = 10% of peak).
#' @param bound_min an `integer``. Minimum bounding-box side length (in cells)
#'   in each dimension; box is widened symmetrically if smaller.
#'
#' @return a `matrix``, same dims as predicted/likelihood: predicted * likelihood
#' outside the bounding box zeroed out.
#'
#' @export
mask_likelihood_array <- function(
  predicted,
  likelihood,
  lon,
  lat,
  par0,
  bound_threshold = 0.1,
  bound_min = 10
) {
  n_row <- length(lat)
  n_col <- length(lon)

  max_likelhood <- max(likelihood, na.rm = TRUE)

  above_thresh <- which(
    likelihood >= bound_threshold * max_likelhood,
    arr.ind = TRUE
  )

  if (nrow(above_thresh) == 0) {
    # ASSUMPTION: fallback to unmasked product rather than returning all-NA
    # or all-zero, matching the "in case of data update led to an empty
    # matrix" comment at the call site — but the actual fallback behavior
    # is unverified.
    return(predicted * likelihood)
  }

  row_range <- range(above_thresh[, "row"])
  col_range <- range(above_thresh[, "col"])

  # ASSUMPTION: par0 is the kernel's odd diameter (as built in
  # gaussian_kernel()), so (par0 - 1) / 2 is the radius in cells — expand
  # the box by this much so a step's worth of diffusion isn't clipped.
  kernel_radius <- (par0 - 1) / 2

  row_min <- row_range[1] - kernel_radius
  row_max <- row_range[2] + kernel_radius
  col_min <- col_range[1] - kernel_radius
  col_max <- col_range[2] + kernel_radius

  # ASSUMPTION: widen symmetrically, one cell at a time each side, until
  # bound_min is met in both dimensions. Untested edge case: this can
  # widen past the grid extent before clamping below.
  while ((row_max - row_min + 1) < bound_min) {
    row_min <- row_min - 1
    row_max <- row_max + 1
  }
  while ((col_max - col_min + 1) < bound_min) {
    col_min <- col_min - 1
    col_max <- col_max + 1
  }

  # Clamp to grid extent
  row_min <- max(1, floor(row_min))
  row_max <- min(n_row, ceiling(row_max))
  col_min <- max(1, floor(col_min))
  col_max <- min(n_col, ceiling(col_max))

  keep <- matrix(FALSE, nrow = n_row, ncol = n_col)
  keep[row_min:row_max, col_min:col_max] <- TRUE

  pred_masked <- predicted
  likelihood_masked <- likelihood
  pred_masked[!keep] <- 0
  likelihood_masked[!keep] <- 0

  out <- pred_masked * likelihood_masked
  return(out)
}
