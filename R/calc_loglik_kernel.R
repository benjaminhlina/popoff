#' Calculates the negative log-likelihood for a Gaussian Kernel
#'
#' This function takes a parameter, a likelihood and a location and determines
#' the negative log-lkelihood of a Guassian Kernel
#'
#' @param parameter a `numerical` parameter
#' @param g a `vector` containing latitude and longitude?
#' @param L a `matrix` containing the liklihood
#' @param normalize a `logical` value
#' @param mask a `logical` value
#' @param bt a `loctical` value that controls back transformation
#'
#' @export

neg_loglik_kernel <- function(
  parameter,
  g,
  L,
  normalize = FALSE,
  mask = TRUE,
  bt = FALSE
) {
  ## Calculating the likelihood of the parameters
  #size=ceiling(par0[1])*2*2+1 #!!!!!
  # gausskern param is gausskern diameter & SD (sigma)
  # kernel radius max of 2 SD (instead of one as Braun, because our D value is way smaller)
  # thus diameter= 2*(2*sd)
  if (isTRUE(bt)) {
    parameter <- .transpose_back(
      parameter_real = parameter,
      parameter_range = parameter_range
    )
  } else {
    parameter
  }

  kernal_size <- ceiling(parameter * 4)
  #sigma=par
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

  hmm_filter <- hmm_filter_1B_deal_NA(g, L, gaus_kern, maskL = FALSE)
  # ----- Negative log-likelihood -----
  neg_log_lik <- -sum(log(hmm_filter$psi))

  return(neg_log_lik)
}
