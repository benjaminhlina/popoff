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
#' @param likelihood_array Arrxy of likelihoods, calculated using
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
#' @param return_maps a `logical` value that returns maps of posteriors and
#' predictions which are large. Defaults to `FALSE`.
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
  progress = TRUE,
  return_maps = FALSE
) {
  t0 <- msg_start("HMM forward filter")
  n_time <- dim(likelihood_array)[3]

  # Preallocate vectors for performance and minimal memory overhead
  ncps_vec <- numeric(n_time - 1)
  updated_vec <- logical(n_time - 1)

  # Optional storage lists only if maps are explicitly requested
  if (isTRUE(return_maps)) {
    post_list <- vector("list", n_time)
    pred_list <- vector("list", n_time)
    post_list[[1]] <- likelihood_array[,, 1]
    pred_list[[1]] <- likelihood_array[,, 1]
  }

  current_posterior <- likelihood_array[,, 1]
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
  # Iterative forward loop (Memory efficient)
  for (t in 2:n_time) {
    pred_t <- imagine::convolution2D(current_posterior, kernel)
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
    current_posterior <- post / psi_t

    # Record metrics
    ncps_vec[t - 1] <- psi_t
    updated_vec[t - 1] <- use_likelihood

    if (isTRUE(return_maps)) {
      post_list[[t]] <- current_posterior
      pred_list[[t]] <- pred_t
    }

    if (isTRUE(progress)) {
      cli::cli_progress_update(id = pb)
    }
  }

  # initial state is the first likelihood surface; psi - normalizing constant per step
  #  is undefined at t = 1

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
  dn <- dimnames(likelihood_array)
  if (!is.null(dn[[3]])) {
    names(ncps_vec) <- dn[[3]][-1]
  }

  n_updated <- sum(updated_vec)
  bad <- !is.finite(ncps_vec) | ncps_vec <= 0
  n_bad <- sum(bad)

  # Floor instead of drop: each collapsed step costs ~690 nats
  safe_ncps <- pmax(ncps_vec, 1e-300)
  safe_ncps[!is.finite(safe_ncps)] <- 1e-300

  log_lik <- safe_ncps |>
    log() |>
    sum()

  neg_log_lik <- log_lik * -1

  # Base results payload (lightweight)
  results <- list(
    ncps = ncps_vec,
    log_lik = log_lik,
    neg_log_lik = neg_log_lik,
    n_bad = n_bad,
    n_updated = n_updated
  )

  # Conditionally attach heavy matrices only if requested
  if (isTRUE(return_maps)) {
    results$normalised_posterior <- array(
      unlist(post_list),
      dim = dim(likelihood_array),
      dimnames = dn
    )
    results$pc_predicted <- array(
      unlist(pred_list),
      dim = dim(likelihood_array),
      dimnames = dn
    )
  }

  msg_end(chr = "HMM forward filter", t0 = t0)
  return(results)
}
