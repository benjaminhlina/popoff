#' HMM forward filter with handling for empty or missing likelihood updates
#'
#' Runs the forward filter over time. Each step convolves the previous
#' posterior with the movement kernel (prediction), then multiplies by the
#' likelihood surface (update) of that step.
#' If the likelihood has no mass, or the update would wipe out the prediction,
#' the prediction is carried forward unchanged. This currently is for one
#' behaviour. Futher itterations of this function may have an argument
#' `behaviour` that allows the user to switch this to be used for two behaviours.
#'
#' @param likelihood_array Array of likelihoods, calculated using
#' [calc_likdepth_transmit()] with the following dimensions
#' `[n_col, n_row, n_time]`.
#' @param kernel Movement kernel (matrix) created by [gaussian_kernel()]
#' @param grid a `data.frame` containing `lon` and `lat` that are evenly spaced
#' across the extent of the study area.
#' @param bound_threshold,bound_min Passed to [mask_likelihood_array()] when
#' `mask_likelihood = TRUE`. `bound_threshold` is fraction of
#' likelihood's max used to decide which cells count as a plausible location
#' (e.g. 0.1 = 10% of peak).
#' @param mask_likelihood Logical. Apply [mask_likelihood_array()] to the update
#' step instead of a plain product.
#' @param progress a `logical` value that turns on and off the progress bar.
#' Default is `TRUE`.
#'
#' @return List with `normalised_posterior` (`[n_col, n_row, n_time]`),
#' `pc_predicted` (post-convolution prediction, same dims) and
#' `ncps`(normalising constant per step, length `n_time - 1`).
#'
#' @export
hmm_forward_filter <- function(
  likelihood_array,
  kernel,
  grid,
  bound_threshold = 0.1,
  bound_min = 10,
  mask_likelihood = FALSE,
  progress = TRUE
) {
  t0 <- msg_start("HMM forward filter")
  n_time <- dim(likelihood_array)[3]

  if (isTRUE(progress)) {
    pb <- cli::cli_progress_bar(
      name = "HMM forward filter",
      format = paste0(
        "{cli::pb_spin} Filtering time step {cli::pb_current}/{cli::pb_total} ",
        "{cli::pb_bar} {cli::pb_percent} | ETA: {cli::pb_eta}"
      ),
      total = n_time - 1,
      clear = FALSE
    )
  }
  # One filter step: previous state to the next state. The state carries the
  # normalised posterior (`normalised_posterior`),
  # the prediction and the normaliser (`normalising constant per step`).
  filter_step <- function(state, t) {
    pred_t <- imagine::convolution2D(state$normalised_posterior, kernel)
    lik_t <- likelihood_array[,, t]

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
    list(
      normalised_posterior = post / psi_t,
      predicted = pred_t,
      psi = psi_t,
      updated = use_likelihood
    )
  }

  # initial state is the first likelihood surface; psi - normalizing constant per step
  #  is undefined at t = 1

  init <- list(
    normalised_posterior = likelihood_array[,, 1],
    predicted = likelihood_array[,, 1],
    psi = NA
  )

  states <- purrr::accumulate(
    .x = seq_len(n_time)[-1],
    .f = filter_step,
    .init = init
  )

  # Look at your use_likelihood guard. When the predicted surface doesn’t
  # overlap the likelihood (overlap == 0), you fall back to post <- pred_t.
  # Then psi_t = sum(pred_t), which is roughly 1, so log(psi) ≈ 0.
  # That step adds nothing to the NLL, as if it were a free pass.
  #
  # This matters for sigma. A small sigma produces a narrow kernel that can’t
  # reach where the likelihood moved, so more steps get skipped, and those
  # steps are never penalized. A too-small sigma can therefore
  # look better than it should. Track this as well:
  #
  # Reject a sigma if n_bad > 0, or if
  #  n_updated is much lower than at larger sigmas.
  # Or compare only sigmas with the same n_updated, or use NLL per
  # updated step (neg_log_lik / n_updated).
  #
  #
  n_updated <- purrr::map_lgl(states[-1], "updated") |>
    sum()

  dn <- dimnames(likelihood_array)

  # ----- get noramlized posteriors ----
  normalised_posterior_calc <- purrr::map(states, "normalised_posterior")

  # ---- get predicted post convuteles -----
  predicted_calc <- purrr::map(
    states,
    "predicted",
  )

  # ---- get normalizing constant per setp -----
  ncps_calc <- purrr::map_dbl(states[-1], "psi") |>
    stats::setNames(dn[[3]][-1])

  # ---- calculate the log likelihood and then calcualte the
  #  negative log likelihood ----

  bad <- !is.finite(ncps_calc) | ncps_calc <= 0
  n_bad <- sum(bad)
  msg_no_mass(n_bad)

  # floor instead of drop: each collapsed step costs ~690 nats
  safe_ncps <- pmax(ncps_calc, 1e-300)
  safe_ncps[!is.finite(safe_ncps)] <- 1e-300

  log_lik <- safe_ncps |>
    log() |>
    sum()

  neg_log_lik <- log_lik * -1

  # ---- results -----
  msg_end(chr = "HMM forward filter", t0 = t0)
  results <- list(
    normalised_posterior = array(
      unlist(
        normalised_posterior_calc
      ),
      dim = dim(likelihood_array),
      dimnames = dn
    ),
    pc_predicted = array(
      unlist(
        predicted_calc
      ),
      dim = dim(likelihood_array),
      dimnames = dn
    ),
    ncps = ncps_calc,
    log_lik = log_lik,
    neg_log_lik = neg_log_lik,
    n_bad = n_bad,
    n_updated = n_updated
  )

  return(results)
}
