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
#' @name mask_functions
#' @keywords internal
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

  max_likelihood <- max(likelihood, na.rm = TRUE)

  # Avoid which(arr.ind = TRUE) on a 3D array: it materialises a huge
  # (n_hits x 3) matrix (~38M rows here). Instead, collapse the date
  # dimension and find the occupied rows/cols directly.
  above <- likelihood >= bound_threshold * max_likelihood
  above[is.na(above)] <- FALSE

  if (length(dim(above)) == 3) {
    above_2d <- apply(above, c(1, 2), any)
  } else {
    above_2d <- above
  }

  if (!any(above_2d)) {
    # Fallback: unmasked product (unverified against call-site intent)
    return(predicted * likelihood)
  }

  rows_hit <- which(rowSums(above_2d) > 0)
  cols_hit <- which(colSums(above_2d) > 0)

  row_range <- range(rows_hit)
  col_range <- range(cols_hit)

  # par0 is the kernel's odd diameter, so (par0 - 1) / 2 is the radius in cells
  kernel_radius <- (par0 - 1) / 2

  row_min <- row_range[1] - kernel_radius
  row_max <- row_range[2] + kernel_radius
  col_min <- col_range[1] - kernel_radius
  col_max <- col_range[2] + kernel_radius

  # Widen symmetrically until bound_min is met (can overshoot grid; clamped below)
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

  # 2D keep mask, recycled across the date dimension if present
  keep <- matrix(FALSE, nrow = n_row, ncol = n_col)
  keep[row_min:row_max, col_min:col_max] <- TRUE

  pred_masked <- predicted
  likelihood_masked <- likelihood

  if (length(dim(likelihood)) == 3) {
    for (k in seq_len(dim(likelihood)[3])) {
      pred_k <- pred_masked[,, k]
      lik_k <- likelihood_masked[,, k]
      pred_k[!keep] <- 0
      lik_k[!keep] <- 0
      pred_masked[,, k] <- pred_k
      likelihood_masked[,, k] <- lik_k
    }
  } else {
    pred_masked[!keep] <- 0
    likelihood_masked[!keep] <- 0
  }

  out <- pred_masked * likelihood_masked
  return(out)
}
