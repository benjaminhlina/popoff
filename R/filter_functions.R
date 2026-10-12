#' Filter step functions
#'
#' These functions provide the backbone of the hidden Markov model as they
#' use the kernel size and evaluate the likelihoods that fall within a forward
#' direction from the previous time step to determine a likely path of movement.
#'
#' @param likellihood this is the `likelihood` array
#' @param kernel the kernel given by [gaussian_kernel()].
#' @param state a given state to eveulate the function in. This is used during
#' the itterative processes.
#' @param time_step a given index value of time step to ittearate over
#'
#' @return a list containing the `normalized posterior`, the `predicated` depth
#' location, the likelihoods as `psi` and whether `use_likelihood` was used.
#'
#' @name filter_functions
#' @keywords internal
filter_step <- function(likelihood, kernel, tate, time_step) {
  pred_t <- imagine::convolution2D(state$normalised_posterior, kernel)
  lik_t <- likelihood[,, time_step]

  # Only update with the likelihood if it has mass and does not zero out the
  # prediction (na.rm because the likelihood can still contain NAs).
  overlap <- sum(pred_t * lik_t, na.rm = TRUE)
  use_likelihood <- sum(lik_t, na.rm = TRUE) > 1e-06 &&
    is.finite(overlap) &&
    overlap != 0

  post <- if (!use_likelihood) {
    pred_t
  } else if (isTRUE(mask_likelihood)) {
    mask_likelihood_array(
      predicted = pred_t,
      likelihood = lik_t,
      lon = grid$lon,
      lat = grid$lat,
      par0 = dim(kernel)[1],
      bound_threshold = bound_threshold,
      bound_min = bound_min
    )
  } else {
    pred_t * lik_t
  }

  psi_t <- sum(post, na.rm = TRUE)
  if (isTRUE(progress)) {
    cli::cli_progress_update(id = pb)
  }
  out <- list(
    normalised_posterior = post / psi_t,
    predicted = pred_t,
    psi = psi_t,
    updated = use_likelihood
  )
  return(out)
}
