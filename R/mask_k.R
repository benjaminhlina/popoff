#' Mask k
#'
#' This fucntion takes coordates and calcualtes eludian distance between
#' locations
#'
#' @param k is `matrix` continaing coordinates
#'
#' @return k as a `matrix` that are elucdian distances between locations with maks
#'
#' @export
mask_k <- function(k) {
  # Get dimensions of the input kernel matrix
  n_cols <- ncol(k)
  n_rows <- nrow(k)

  # Calculate the center index radius (e.g., for a 5x5 matrix, center offset is 2)
  center_offset <- trunc(n_cols / 2) + 1

  # Create X-coordinate offsets relative to the center (column-wise grid)
  coord_kx <- matrix(
    data = rep(1:n_cols, each = n_rows),
    ncol = n_cols,
    nrow = n_rows
  ) -
    center_offset

  # Create Y-coordinate offsets relative to the center (row-wise grid)
  # FIXED BUG: Using times instead of each to correctly build the vertical axis
  coord_ky <- matrix(
    data = rep(1:n_rows, times = n_cols),
    ncol = n_cols,
    nrow = n_rows
  ) -
    center_offset

  coord <- cbind(as.vector(coord_kx), as.vector(coord_ky))
  coord <- apply(coord, 1, function(x) {
    ed <- sqrt(x[1]**2 + x[2]**2)
    return(ed)
  }) # distance euclidienne avec le point de coordonnes 0 0

  coord <- matrix(data = coord, ncol = ncol(k), nrow = nrow(k))
  coord[round(coord) > (trunc(ncol(k) / 2))] <- NA
  mask <- which(is.na(coord))
  k[mask] <- 0
  k <- k / sum(k, na.rm = T)
  return(k)
}


#' Gaussian Kernel
#'
#' Determines the Gaussian Kernel for an area
#'
#' @param size a `numerical` value that is the size of the areas to calculate the
#' Kernel
#' @param sigma a `numerical` value that is the standard deviation
#' @param muadv a `numerical` value that adds to the mean (i.e, mu)
#' @param norm a `logical` that checks whether to normalize the Kernel estimate
#'
#' @return a `matrix` that is Gaussian Kernel over an area.
#'
#' @export
#'
gaussian_kernel <- function(size, sigma, muadv = 0, norm = FALSE) {
  if (round(size) < 1) {
    size <- 1
  }
  x <- 1:round(size)
  mu <- c(mean(x), mean(x)) + muadv
  # options(digits = 5)
  fx <- exp(-0.5 * ((x - mu[1]) / sigma)^2) / sqrt((2 * pi) * sigma**2)
  fy <- exp(-0.5 * ((x - mu[2]) / sigma)^2) / sqrt((2 * pi) * sigma**2)
  fx[!is.finite(fx)] <- 0
  fy[!is.finite(fy)] <- 0
  kern <- (fx %*% t(fy))
  if (isTRUE(norm)) {
    kern <- kern / (sum(kern, na.rm = T))
  } else {
    kern
  }
  kern[is.nan(kern)] <- 0

  return(kern)
}
