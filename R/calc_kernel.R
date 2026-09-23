#' Gaussian Kernel
#'
#' Determines the Gaussian Kernel for an area
#'
#' @param kernal_size a `numerical` value that is the size of the areas
#' to calculate the Gaussian Kernel.
#' @param sigma a `numerical` value that is the standard deviation
#' @param muadv a `numerical` value that adds to the mean (i.e, mu)
#' @param normalize a `logical` that checks whether to normalize the Kernel estimate.
#' Default is `FALSE`.
#' @param mask a `logical` that creates a spatial mask based on the elucidan
#' distance. Default is `TRUE`
#'
#' @return a `matrix` that is Gaussian Kernel over an area.
#'
#' @export

gaussian_kernel <- function(
  kernal_size,
  sigma,
  muadv = 0,
  normalize = FALSE,
  mask = TRUE
) {
  error_numeric(kernal_size)
  error_numeric(sigma)
  error_numeric(muadv)
  error_logical(normalize)
  error_logical(mask)
  error_kernel(mask, sigma)
  t0 <- msg_start("Gaussian kernel")

  if (round(kernal_size) < 1) {
    kernal_size <- 1
  }
  x <- 1:round(kernal_size)
  mu <- c(mean(x), mean(x)) + muadv

  # cal elducidan distance from mu, x and sigma
  fx <- exp(-0.5 * ((x - mu[1]) / sigma)^2) / sqrt((2 * pi) * sigma**2)
  fy <- exp(-0.5 * ((x - mu[2]) / sigma)^2) / sqrt((2 * pi) * sigma**2)
  fx[!is.finite(fx)] <- 0
  fy[!is.finite(fy)] <- 0
  kern <- (fx %*% t(fy))

  if (isTRUE(normalize)) {
    kern <- kern / (sum(kern, na.rm = T))
  } else {
    kern
  }

  kern[is.nan(kern)] <- 0

  if (isTRUE(mask)) {
    # get dimensions of the input kernel matrix
    n_cols <- ncol(kern)
    n_rows <- nrow(kern)

    # calculate the center index radius (e.g., for a 5x5 matrix, center offset is 2)
    center_offset <- trunc(n_cols / 2) + 1

    # create x-coordinate offsets relative to the center (column-wise grid)
    coord_kx <- matrix(
      data = rep(1:n_cols, each = n_rows),
      ncol = n_cols,
      nrow = n_rows
    ) -
      center_offset

    # create y-coordinate offsets relative to the center (row-wise grid)
    coord_ky <- matrix(
      data = rep(1:n_rows, times = n_cols),
      ncol = n_cols,
      nrow = n_rows
    ) -
      center_offset

    # calculate Euclidean distance from (0,0) center for every matrix cell (vectorized)
    coord <- sqrt(coord_kx^2 + coord_ky^2)

    # mask cells where the rounded distance exceeds the maximum allowed radius
    coord[round(coord) > trunc(n_cols / 2)] <- NA
    mask_kern <- which(is.na(coord))

    # zero out the masked positions in the original kernel
    kern[mask_kern] <- 0

    if (isTRUE(normalize)) {
      # normalize kernel values so the remaining elements sum to 1
      kern <- kern / sum(kern, na.rm = TRUE)
    } else {
      kern
    }
  }
  msg_end("Calculated Gaussian kernel", t0 = t0)
  return(kern)
}
