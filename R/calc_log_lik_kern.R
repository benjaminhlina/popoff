neg_log_lik_kern <- function(par, g, L, norm = FALSE, mask = TRUE, bt = FALSE) {
  ## Calculating the likelihood of the parameters
  #size=ceiling(par0[1])*2*2+1 #!!!!!
  # gausskern param is gausskern diameter & SD (sigma)
  # kernel radius max of 2 SD (instead of one as Braun, because our D value is way smaller)
  # thus diameter= 2*(2*sd)
  if (isTRUE(bt)) {
    par <- .transpose_back(
      parameter_real = parameter,
      parameter_range = parameter_range
    )
  } else {
    par
  }

  kernal_size <- ceiling(par * 4)
  #sigma=par
  if (kernal_size %% 2 == 0) {
    kernal_size <- kernal_size + 1
  } # si resultat pair on rajoute 1 pour assurer une cellule centrale (point de depart)
  # this ensure to gather around 90-95% of diffusion

  # ajout 22/03
  # avoid that the fish go further in diagnal that in straight line
  gaus_kern <- gaussian_kernel(
    kernal_size = kernal_size,
    sigma = par,
    muadv = 0,
    norm = norm,
    mask = mask
  )

  hmm_filter <- hmm.filter.1B.dealNA(g, L, gaus_kern, maskL = F)
  # ----- Negative log-likelihood -----
  neg_log_lik <- -sum(log(hmm_filter$psi))
  return(neg_log_lik)
}
