hmm_forward_filter <- function(
  likelihood_array,
  kernel,
  bound_threshold = 0.1,
  bound_min = 10,
  mask_likelihood = TRUE,
  progress = TRUE,
  return_maps = FALSE
) {
  t0 <- msg_start("HMM forward filter")
  dims <- dim(likelihood_array)
  n_time <- dims[3]

  # Preallocate vectors for performance and minimal memory overhead
  ncps_vec <- numeric(n_time - 1)
  updated_vec <- logical(n_time - 1)

  dn <- dimnames(likelihood_array)
  # Optional storage lists only if maps are explicitly requested
  if (isTRUE(return_maps)) {
    post_list <- vector("list", n_time)
    pred_list <- vector("list", n_time)
    post_list[[1]] <- likelihood_array[,, 1]
    pred_list[[1]] <- likelihood_array[,, 1]
  }

  # ----- grab land values ------
  land <- is.na(likelihood_array[,, 1])

  land_arr <- array(land, dim = dims, dimnames = dn)

  current_posterior <- likelihood_array[,, 1]
  current_posterior[land] <- 0

  if (isTRUE(progress)) {
    msg_progress(
      n_time = n_time,
      name = "HMM forward filtering",
      type = "Filtering"
    )
  }
  # One filter step: previous state to the next state. The state carries the
  # normalised posterior (`normalised_posterior`),
  # the prediction and the normaliser (`normalising constant per step`).
  # Iterative forward loop (Memory efficient)
  for (t in 2:n_time) {
    pred_t <- imagine::convolution2D(current_posterior, kernel)
    pred_t[land] <- 0

    lik_t <- likelihood_array[,, t]
    lik_t[land] <- 0

    # pred_mass <- sum(pred_t)
    # if (!is.finite(pred_mass) || pred_mass <= 0) {
    #   cli::cli_abort(
    #     "Prediction has no mass at step {t}: the state is dead before the likelihood is applied."
    #   )
    # }
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
        par0 = dim(kernel)[1],
        bound_threshold = bound_threshold,
        bound_min = bound_min
      )
    } else {
      pred_t * lik_t
    }

    psi_t <- sum(post, na.rm = TRUE)
    if (!is.finite(psi_t) || psi_t <= 0) {
      # Likelihood step collapsed: carry the prediction forward so the
      # filtered state stays a valid surface. psi_t is still recorded
      # (floored later), so the step is penalised in the NLL.
      use_likelihood <- FALSE
      post <- pred_t
      current_posterior <- post / sum(post, na.rm = TRUE)
    } else {
      current_posterior <- post / psi_t
    }

    # Record metrics
    ncps_vec[t - 1] <- psi_t
    updated_vec[t - 1] <- use_likelihood

    if (isTRUE(return_maps)) {
      post_list[[t]] <- current_posterior
      pred_list[[t]] <- pred_t
    }

    if (isTRUE(progress)) {
      cli::cli_progress_update()
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

  if (!is.null(dn[[3]])) {
    names(ncps_vec) <- dn[[3]][-1]
  }

  n_updated <- sum(updated_vec)
  bad <- !is.finite(ncps_vec) | ncps_vec <= 0
  n_bad_days <- sum(bad)

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
    n_bad_days = n_bad_days,
    n_updated = n_updated
  )

  # Conditionally attach heavy matrices only if requested
  if (isTRUE(return_maps)) {
    results$normalised_posterior <- array(
      unlist(post_list),
      dim = dims,
      dimnames = dn
    )
    results$pc_predicted <- array(
      unlist(pred_list),
      dim = dims,
      dimnames = dn
    )

    results$normalised_posterior[land_arr] <- NA
    results$pc_predicted[land_arr] <- NA
  }

  msg_end(chr = "HMM forward filter", t0 = t0)
  return(results)
}
