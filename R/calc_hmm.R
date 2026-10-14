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
  dn <- dimnames(likelihood_array)
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

#' HMM backward filter with handling for empty or missing likelihood updates
#'
#' Runs the backward recursion of a one-state, grid-based hidden Markov model
#' (HMM) for geolocation. It combines the forward-filtered probability surfaces
#' with the movement kernel (i.e., gaussian kernel) and
#' the observation likelihood to give smoothed location probabilities at every
#' time step. Strictly speaking this is the backward *smoothing* step that follows
#' the forward filter.
#'
#'
#' @param forward_filter_obj A `list` containing at least two elements,
#'  both numeric arrays with dimensions `[x, y, t_step]`: The first
#' object being `normalised_posterior` which is the forward-filtered posterior
#' probability surfaces derived from [hmm_forward_filter()].
#'
#' The second is `pc_predicted` whcih is one-step-ahead predicted surfaces
#' the filtered surface after applying the movement kernel0.
#' @param kernel A numeric `matrix`` giving the movement (diffusion) kernel used
#' for the 2D convolution created by [gaussian_kernel()]. This is the same kernel
#' used in `hmm_forward_filter()`.
#' @param likelihood_array A `numeric` array with dimensions `[x, y, t_step]` holding the
#' observation likelihood surfaces. It must have the same dimensions as both
#' `normalised_posterior` and `pc_predicted`.
#'
#' @details
#' All arrays use the layout `[x, y, t_step]`, with time as the last dimension.
#' `n_t_steps` is the number of time steps.
#'
#' The final time step is seeded from the observation likelihood
#' `likelihood_array` rather than from the forward-filtered surface.
#' For the last step, the smoothed surface is `likelihood_array[, , n_t_steps]`
#' normalised to sum to 1. The step before it is
#' `likelihood_array[, , n_t_steps] * normalised_posterior[, , n_t_steps_t - 1]`,
#'  also normalised. For each earlier time_step (i.e. `t_step`)
#' from `n_t_steps - 1` down to 2, the function:
#' \enumerate{
#'   \item divides the smoothed surface at `t_step` by the predicted surface at `t_step`,
#'   \item convolves that ratio with the movement kernel `kernel`
#'     ([imagine::convolution2D()]), and
#'   \item multiplies the result by the forward-filtered surface at `t_step - 1`
#'     and normalises it to sum to 1.
#' }
#'
#' Normalisation ignores `NA` cells (`na.rm = TRUE`), so `NA` cells, for
#' example land or cells outside the model domain, do not affect the sums.
#' Cells where `pc_predicted` is exactly 0 produce `NaN` or `Inf` in the ratio.
#' Make sure `pc_predicted` is strictly positive wherever `smooth` is non-zero.
#' @return A numeric array with the same dimensions as `normalised_posterior`,
#' `[x, y, t_steps]`, containing the smoothed probability surfaces.
#' Each time slice sums to 1 (ignoring `NA` cells), except slice 1, which is
#' left at 0 because the loop stops at `t_step = 2`.
#'
#' @references
#' Pedersen, M. W., Righton, D., Thygesen, U. H., Andersen, K. H., and
#' Madsen, H. (2008). Geolocation of North Sea cod (*Gadus morhua*) using
#' hidden Markov models and behavioural switching. *Canadian Journal of
#' Fisheries and Aquatic Sciences*, 65(11), 2367-2377.
#'
#' Braun, C. D., Galuardi, B., and Thorrold, S. R. (2018). HMMoce: An R
#' package for improved geolocation of archival-tagged fishes using a hidden
#' Markov method. *Methods in Ecology and Evolution*, 9(5), 1212-1220.
#'
#
#' @export

hmm_backward_smooth <- function(
  forward_filter_obj,
  kernel,
  likelihood_array,
  progress = TRUE
) {
  t0 <- msg_start("HMM backwards smoothing")

  error_forward_filter_obj(forward_filter_obj)

  normalised_posterior <- forward_filter_obj$normalised_posterior
  pc_predicted <- forward_filter_obj$pc_predicted
  dims <- dim(forward_filter_obj$normalised_posterior)
  dn <- dimnames(normalised_posterior)
  time_step <- dims[3]

  smooth <- array(0, dim = dims, dimnames = dn)

  smooth[,, time_step] <- normalised_posterior[,, time_step]
  # premier pas de temps du smooth est le dernier du backward !!
  time_seq <- time_step:2

  if (isTRUE(progress)) {
    msg_progress(
      n_time = time_step,
      name = "HMM backwards smoothing",
      type = "Smoothing"
    )
  }

  for (t in time_seq) {
    rat <- smooth[,, t] / pc_predicted[,, t]
    rat[!is.finite(rat)] <- 0

    post <- imagine::convolution2D(rat, kernel)

    if (t == time_step) {
      post <- likelihood_array[,, t]
      smooth[,, t] <- post / sum(post, na.rm = TRUE)
      post <- post * normalised_posterior[,, t - 1]
    } else {
      post <- post * normalised_posterior[,, t - 1]
    }

    smooth[,, t - 1] <- post / sum(post, na.rm = TRUE)

    if (isTRUE(progress)) cli::cli_progress_update()
  }

  if (isTRUE(progress)) {
    cli::cli_progress_done()
  }
  msg_end("HMM backwards smoothing", t0)
  return(smooth)
}


hmm_backward_smooth_old <- function(
  forward_filter_obj,
  kernel,
  likelihood_array,
  progress = TRUE
) {
  # ------ messages and error-----
  t0 <- msg_start("HMM backwards smoothing")
  error_forward_filter_obj(forward_filter_obj)

  # ----- get objects we need -----
  normalised_posterior <- forward_filter_obj$normalised_posterior
  pc_predicted <- forward_filter_obj$pc_predicted
  dims <- dim(normalised_posterior)
  dn <- dimnames(normalised_posterior)
  time_step <- dims[3]
  # row <- dims[1]
  # col <- dims[2]
  # ----- create smooth array ----
  smooth <- array(0, dim = dims, dimnames = dn)

  smooth[,, time_step] <- normalised_posterior[,, time_step]

  time_seq <- time_step:2
  if (isTRUE(progress)) {
    msg_progress(
      n_time = time_step,
      name = "HMM backwards smoothing",
      type = "Smoothing"
    )
  }
  for (t in time_seq) {
    rat <- smooth[,, t] / pc_predicted[,, t]
    rat[!is.finite(rat)] <- 0 # 0/0 and x/0 -> 0 so NaN can't spread

    post <- imagine::convolution2D(rat, kernel)
    post <- post * normalised_posterior[,, t - 1]

    post[!is.finite(post)] <- 0

    fac <- sum(post, na.rm = TRUE)
    if (!is.finite(fac) || fac <= 0) {
      cli::cli_abort(
        "Smoothing collapsed at time step {t - 1}: posterior sums to {fac}."
      )
    }

    smooth[,, t - 1] <- post / fac

    if (isTRUE(progress)) cli::cli_progress_update()
  }

  if (isTRUE(progress)) {
    cli::cli_progress_done()
  }
  msg_end("HMM backwards smoothing", t0)

  return(smooth)
}
