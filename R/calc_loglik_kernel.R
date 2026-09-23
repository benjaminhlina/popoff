#' Calculates the negative log-likelihood for a Gaussian Kernel
#'
#' This function takes a parameter, a likelihood and a grid and determines
#' the negative log-lkelihood of a Guassian Kernel
#'
#' @param parameter a `numerical` parameter
#' @param grid a `vector` containing latitude and longitude?
#' @param likelihood a `matrix` containing the liklihood
#' @param normalize a `logical` value that controls whether
#' Guassian Kernel values are normalized. Default is `FALSE`.
#' @param mask a `logical` value that controls whether Guassian Kernel values
#'  are masked to a specific geolocation. Default value is `TRUE`.
#' @param back_transform a `loctical` value that controls back transformation
#' @param parameter_range a `vector` containing two `numerical` values
#' indicating the range of the `parameter`. Defaults to `NULL` and is only
#' needed if `back_transform` is `TRUE`.
#'
#' @export

neg_loglik_kernel <- function(
  parameter,
  grid,
  likelihood,
  normalize = FALSE,
  mask = TRUE,
  back_transform = FALSE,
  parameter_range = NULL,
  mask_likelihood = FALSE
) {
  error_matrix(likelihood)
  error_logical(normalize)
  error_logical(mask)
  error_logical(mask_likelihood)
  error_logical(back_transform)
  error_transform(
    back_transform = back_transform,
    parameter_range = parameter_range
  )
  ## Calculating the likelihood of the parameters
  #size=ceiling(par0[1])*2*2+1 #!!!!!
  # gausskern param is gausskern diameter & SD (sigma)
  # kernel radius max of 2 SD (instead of one as Braun, because
  # our D value is way smaller)
  # thus diameter= 2*(2*sd)
  if (isTRUE(back_transform)) {
    parameter <- .transpose_back(
      parameter_real = parameter,
      parameter_range = parameter_range
    )
  } else {
    parameter
  }

  kernal_size <- ceiling(parameter * 4)

  if (kernal_size %% 2 == 0) {
    kernal_size <- kernal_size + 1
  } # si resultat pair on rajoute 1 pour assurer une cellule centrale (point de depart)
  # this ensure to gather around 90-95% of diffusion

  # ajout 22/03
  # avoid that the fish go further in diagnal that in straight line
  gaus_kern <- gaussian_kernel(
    kernal_size = kernal_size,
    sigma = parameter,
    muadv = 0,
    normalize = normalize,
    mask = mask
  )

  hmm_filter <- hmm_filter_1B_deal_NA(
    grid = grid,
    likelihood_array = likelihood,
    kernal = gaus_kern,
    mask_likelihood = mask_likelihood
  )
  # ----- Negative log-likelihood -----
  neg_log_lik <- -sum(log(hmm_filter$psi))

  return(neg_log_lik)
}
