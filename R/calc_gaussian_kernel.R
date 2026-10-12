#' Gaussian Kernel
#'
#' Determines the Gaussian Kernel for a given area.
#'
#' @param sigma a `numerical` value that is the standard deviation
#' @param muadv a `numerical` value that adds to the mean (i.e, mu)
#' @param normalize a `logical` that checks whether to normalize the Kernel estimate.
#' Default is `FALSE`.
#' @param mask a `logical` that creates a spatial mask based on the elucidan
#' distance. Default is `TRUE`
#' @param back_transform a `loctical` value that controls back transformation
#' @param sigma_range a `vector` containing two `numerical` values
#' indicating the range of the `sigma`. Defaults to `NULL` and is only
#' needed if `back_transform` is `TRUE`.
#'
#' @return a `matrix` that is Gaussian Kernel over an area.
#'
#' @export

gaussian_kernel <- function(
  sigma,
  muadv = 0,
  normalize = FALSE,
  mask = TRUE,
  back_transform = FALSE,
  sigma_range = NULL
) {
  error_numeric_min(sigma)
  error_numeric(muadv)
  # error_numeric(truncate)
  error_logical(normalize)
  error_logical(mask)
  error_kernel(mask, sigma)
  error_logical(back_transform)
  error_transform(
    back_transform = back_transform,
    parameter_range = sigma_range
  )

  # if (!is.null(radius)) {
  #   error_numeric_min(radius)
  # }

  # the following blog post walks th
  if (isTRUE(back_transform)) {
    sigma <- .transpose_back(
      parameter_real = sigma,
      parameter_range = sigma_range
    )
  }

  # ----- old method -----
  kernel_size <- ceiling(sigma * 4)
  if (kernel_size %% 2 == 0) {
    kernel_size <- kernel_size + 1
  }
  x <- 1:round(kernel_size)

  mu <- c(mean(x), mean(x)) + muadv

  # ----- new method -----
  # if (is.null(radius)) {
  #   radius <- ceiling(truncate * sigma)
  # }

  # kernel_size <- 2 * radius + 1

  # x <- seq_len(kernel_size)
  # mu <- (kernel_size + 1) / 2 + muadv

  # sigma * 4 is the kernel’s full width, so the radius is only about 2 sigma.
  # In 2D, with your circular mask, that keeps only about 86% of the mass
  #  (1 - exp(-2)), not the 90-95% in your comment.
  # The ceiling and the odd-size bump then change the
  # effective cutoff as sigma moves:
  # sigma	size	radius	radius in sigma units
  # 1.01	5	2	1.98σ
  # 1.25	5	2	1.6σ
  # 1.30	7	3	2.3σ
  #
  # Sigma 1.25 and sigma 1.3 get quite different truncations.
  # After normalization the kernel shape changes in a sawtooth pattern, which is
  # the source of the neg_likelihood jumps I flagged earlier.
  # The truncation is cutting off enough mass that the jumps show up in the NLL.

  # si resultat pair on rajoute 1 pour assurer une cellule centrale (point de depart)
  # this ensure to gather around 90-95% of diffusion
  # The round(kernel_size) < 1 check is redundant since ceiling()
  # already returns an integer, and for sigma > 0 it’s at least 1.
  # It only matters if sigma is 0 or negative, and in that case you’d rather
  # error out than silently return a 1x1 kernel.
  # ajout 22/03
  # avoid that the fish go further in diagnal that in straight line

  # muadv is added to both x and y, so it shifts the kernel diagonally only.
  # If you want independent x/y drift, it needs to be a length-2 vector.
  # mu <- c(mean(x), mean(x)) + muadv

  # cal elducidan distance from mu, x and sigma
  fx <- exp(-0.5 * ((x - mu[1]) / sigma)^2) / sqrt((2 * pi) * sigma**2)
  fy <- exp(-0.5 * ((x - mu[2]) / sigma)^2) / sqrt((2 * pi) * sigma**2)
  fx[!is.finite(fx)] <- 0
  fy[!is.finite(fy)] <- 0
  kern <- (fx %*% t(fy))

  if (isTRUE(normalize)) {
    kern <- kern / (sum(kern, na.rm = T))
  }
  # ern[is.nan(kern)] <- 0 comes after normalization. If sum(kern) is 0,
  # normalizing produces NaN, so the cleanup is in the right place, but
  # the same risk exists in the second normalization after masking.

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
  }

  if (isTRUE(normalize)) {
    kern <- kern / (sum(kern, na.rm = TRUE))
  }
  # ern[is.nan(kern)] <- 0 comes after normalization. If sum(kern) is 0,
  # normalizing produces NaN, so the cleanup is in the right place, but
  # the same risk exists in the second normalization after masking.

  kern[is.nan(kern)] <- 0

  return(kern)
}
