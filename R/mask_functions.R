#' Mask predicted/likelihood surfaces to a threshold-defined bounding region
#'
#'
#' @param predicted a `matrix``. Predicted spatial distribution at time t.
#' @param likelihood a `matrix``. Observation likelihood at time t.
#' @param par0 a `numeric`. Kernel diameter (odd), used to derive a movement
#' radius so the bounding box isn't clipped tighter than the diffusion
#' kernel would allow in one time step.
#' @param bound_threshold a `numeric`.
#' Fraction of likelihood's max used to decide which cells count
#' as a plausible  location (e.g. 0.1 = 10% of peak).
#' @param bound_min an `integer``. Minimum bounding-box side length (in cells)
#' in each dimension; box is widened symmetrically if smaller.
#'
#' @return a `matrix``, same dims as predicted/likelihood: predicted * likelihood
#' outside the bounding box zeroed out.
#'
#' @name mask_functions
#' @keywords internal
mask_likelihood_array <- function(
  predicted,
  likelihood,
  par0,
  bound_threshold = 0.1,
  bound_min = 10
) {
  if (!(identical(dim(predicted), dim(likelihood)))) {
    cli::cli_abort(
      "The dimensions of predicted are {.val {dim(predicted}} and do not match {.val {dim(likelihood)}}"
    )
  }

  # ------ gab number of rows and columns
  n_row <- nrow(likelihood)
  n_col <- ncol(likelihood)

  # ----- max likelihood this should be 1 ------
  max_likelihood <- max(likelihood, na.rm = TRUE)

  # ----- create threashoold if likelihoohds are above 0.1
  above <- likelihood >= bound_threshold * max_likelihood
  above[is.na(above)] <- FALSE

  if (!any(above)) {
    # Fallback: unmasked product (unverified against call-site intent)
    return(predicted * likelihood)
  }

  # ------ row and col range for cells above ------
  row_range <- range(which(rowSums(above) > 0))
  col_range <- range(which(colSums(above) > 0))

  # par0 is the kernel's odd diameter, so (par0 - 1) / 2 is the radius in cells
  kernel_radius <- (par0 - 1) / 2

  # ------ we then get the row min and max based on the kernel_readisu
  row_min <- row_range[1] - kernel_radius
  row_max <- row_range[2] + kernel_radius
  col_min <- col_range[1] - kernel_radius
  col_max <- col_range[2] + kernel_radius

  # Widen symmetrically until bound_min is met (clamped below)
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

  rows <- row_min:row_max
  cols <- col_min:col_max

  # Zero outside the box, product inside it
  out <- matrix(0, nrow = n_row, ncol = n_col)
  # ----- return out masked predicted * likelihood-----
  out[rows, cols] <- predicted[rows, cols] * likelihood[rows, cols]
  return(out)
}
