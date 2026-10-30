#' Convert a known-events data frame to time-step indices
#'
#' Internal helper. Matches the dates in `known_events` to the third
#' dimname (dates) of a `[n_col, n_row, n_time]` array.
#'
#' @param known_events A `data.frame` with a date column (`Date`, `POSIXct` or
#' ISO `character`), or `NULL`.
#' @param time_names Character vector, `dimnames(array)[[3]]`.
#' @param date_col Name of the date column in `known_events`.
#'
#' @return Integer vector of time-step indices (possibly empty).
#' @noRd
known_event_steps <- function(known_events, time_names, date_col = "date") {
  if (is.null(known_events)) {
    return(integer(0))
  }
  if (!is.data.frame(known_events) || !date_col %in% names(known_events)) {
    cli::cli_abort(
      "{.arg known_events} must be a data frame with a {.field {date_col}} column."
    )
  }
  if (is.null(time_names)) {
    cli::cli_abort(
      "{.arg known_events} needs dates in {.code dimnames(likelihood_array)[[3]]}."
    )
  }

  ev_dates <- as.character(as.Date(known_events[[date_col]]))
  step_dates <- as.character(as.Date(time_names))

  idx <- match(ev_dates, step_dates)

  if (anyNA(idx)) {
    cli::cli_warn(
      "{sum(is.na(idx))} known event date{?s} not found in the array's time dimension and ignored: {ev_dates[is.na(idx)]}."
    )
  }

  sort(unique(idx[!is.na(idx)]))
}

#' Chebyshev distance (in grid cells) from every cell to a known location
#'
#' Internal helper. The known cell(s) are the max of the likelihood slice
#' (the "1 at the nearest cell"). Chebyshev distance matches a square kernel
#' footprint: one step moves at most `r` cells in x and y.
#'
#' @param lik_slice `[n_col, n_row]` likelihood matrix for a known date.
#' @return Matrix of distances, same dims as `lik_slice`.
#' @noRd
known_reach_distance <- function(lik_slice) {
  hit <- which(lik_slice == max(lik_slice, na.rm = TRUE), arr.ind = TRUE)
  rr <- row(lik_slice)
  cc <- col(lik_slice)
  d <- matrix(Inf, nrow(lik_slice), ncol(lik_slice))
  for (i in seq_len(nrow(hit))) {
    d <- pmin(d, pmax(abs(rr - hit[i, 1]), abs(cc - hit[i, 2])))
  }
  d
}

#' HMM forward filter with handling for empty or missing likelihood updates
#'
#' Runs the forward filter over time. Each step convolves the previous
#' posterior with the movement kernel (prediction), then multiplies by the
#' likelihood surface (update) of that step.
#' If the likelihood has no mass, or the update would wipe out the prediction,
#' the prediction is carried forward unchanged, **except** on dates supplied in
#' `known_events`. On those dates the likelihood is treated as a known location
#' (e.g. a MRPAT or acoustic detection, or the tagging / pop-off location) and
#' is always applied: no mass guard and no masking. If the prediction cannot
#' reach the known cell, the filter snaps to the known location and the step is
#' counted as bad (so its normalising constant is floored and penalised in
#' `neg_log_lik`). This penalises sigmas whose kernel is too narrow to reach
#' where the animal is known to have been.
#'
#' @param likelihood_array Array of likelihoods, calculated using
#' [calc_likdepth_transmit()] with the following dimensions
#' `[n_col, n_row, n_time]`. The third dimname must hold the dates.
#' @param kernel Movement kernel (matrix) created by [gaussian_kernel()]
#' @param bound_threshold,bound_min Passed to [mask_likelihood_array()] when
#' `mask_likelihood = TRUE`. `bound_threshold` is fraction of
#' likelihood's max used to decide which cells count as a plausible location
#' (e.g. 0.1 = 10% of peak).
#' @param mask_likelihood Logical. Apply [mask_likelihood_array()] to the update
#' step instead of a plain product. Default is `TRUE`. Ignored on known-event
#' dates.
#' @param progress a `logical` value that turns on and off the progress bar.
#' Default is `TRUE`.
#' @param return_maps a `logical` value that returns maps of posteriors and
#' predictions which are large. Defaults to `TRUE`.
#' @param known_events Optional `data.frame` with one row per known-location
#' date. Dates are matched to `dimnames(likelihood_array)[[3]]`. Default `NULL`.
#' @param date_col Name of the date column in `known_events`. Default `"date"`.
#' @param reach_constraint Logical. When `known_events` is supplied, each
#' day's posterior is restricted to cells from which the next known location
#' is still reachable: Chebyshev distance to that location must be at most
#' `(steps remaining) * r`, where `r` is the kernel radius in cells. This
#' makes the forward filter informed by future known locations (e.g. MRPAT
#' and pop-off), not only by past data. The constraint is applied after the
#' normalising constant is recorded, so `ncps` still reflects the
#' unconstrained filter. If the constraint would remove all mass on a day it
#' is not applied that day and the day is counted in `n_unreachable`.
#' Default `TRUE`.
#'
#' @return List with `normalised_posterior` (`[n_col, n_row, n_time]`),
#' `pc_predicted` (post-convolution prediction, same dims),
#' `ncps` (normalising constant per step, length `n_time - 1`),
#' `n_unreachable` (days where no cell satisfied the reachability constraint;
#' any value above 0 means this sigma / kernel cannot connect the known
#' locations and should be rejected),
#' `infeasible` (`TRUE` if `n_unreachable > 0` or `n_known_missed > 0`),
#' `nll_per_update` (`neg_log_lik / n_updated`),
#' `n_known` (number of known-event steps applied after step 1),
#' `n_known_missed` (known steps the prediction could not reach) and
#' `known_idx` (time-step indices of known events, used by
#' [hmm_backward_smooth()]).
#'
#' @export
hmm_forward_filter_new_2 <- function(
  likelihood_array,
  kernel,
  bound_threshold = 0.1,
  bound_min = 10,
  mask_likelihood = TRUE,
  progress = TRUE,
  return_maps = TRUE,
  known_events = NULL,
  date_col = "date",
  reach_constraint = TRUE
) {
  error_likelihood_array(likelihood_array)
  error_logical(reach_constraint)
  error_matrix(kernel)
  error_numeric(bound_threshold)
  error_numeric(bound_min)
  error_logical(mask_likelihood)
  error_logical(progress)
  error_logical(return_maps)

  t0 <- msg_start("HMM forward filter")
  dims <- dim(likelihood_array)
  n_time <- dims[3]

  # Preallocate vectors for performance and minimal memory overhead
  ncps_vec <- numeric(n_time - 1)
  updated_vec <- logical(n_time - 1)
  known_missed_vec <- logical(n_time - 1)
  unreachable_vec <- logical(n_time - 1)

  dn <- dimnames(likelihood_array)

  # ----- known-location dates -> time-step indices ------
  known_idx <- known_event_steps(known_events, dn[[3]], date_col)

  # ----- reachability bridge: distance to each known location ------
  use_reach <- isTRUE(reach_constraint) && length(known_idx) > 0
  if (use_reach) {
    reach_r <- (dim(kernel)[1] - 1) / 2 # max cells moved per step
    reach_dist <- lapply(known_idx, function(k) {
      s <- likelihood_array[,, k]
      s[is.na(s)] <- 0
      known_reach_distance(s)
    })
    names(reach_dist) <- as.character(known_idx)
    # next known step strictly after each t (NA if none)
    next_known <- vapply(
      seq_len(n_time),
      function(t) {
        k <- known_idx[known_idx > t]
        if (length(k)) k[1] else NA_integer_
      },
      integer(1)
    )
  }

  # ----- grab land values ------
  land <- is.na(likelihood_array[,, 1])

  land_arr <- array(land, dim = dims, dimnames = dn)

  current_posterior <- likelihood_array[,, 1]
  current_posterior[land] <- 0

  current_posterior <- current_posterior / sum(current_posterior, na.rm = TRUE)
  if (isTRUE(return_maps)) {
    post_list <- vector("list", n_time)
    pred_list <- vector("list", n_time)
    post_list[[1]] <- current_posterior # Uses the normalized version
    pred_list[[1]] <- current_posterior # Uses the normalized version
  }

  msg_progress(
    progress = progress,
    n_time = n_time,
    name = "HMM forward filtering",
    type = "Filtering"
  )

  for (t in 2:n_time) {
    pred_t <- imagine::convolution2D(current_posterior, kernel)
    pred_t[land] <- 0

    lik_t <- likelihood_array[,, t]
    lik_t[land] <- 0

    if (t %in% known_idx) {
      # ----- known location: always apply the likelihood -----
      lik_sum <- sum(lik_t, na.rm = TRUE)
      if (!is.finite(lik_sum) || lik_sum <= 0) {
        cli::cli_abort(
          "Known event at time step {t} ({dn[[3]][t]}) has an empty likelihood surface."
        )
      }

      post <- pred_t * lik_t
      psi_t <- sum(post, na.rm = TRUE)
      use_likelihood <- TRUE

      if (!is.finite(psi_t) || psi_t <= 0) {
        # Kernel could not reach the known cell: snap to the known location
        # and let psi_t = 0 be floored/penalised downstream.
        known_missed_vec[t - 1] <- TRUE
        psi_t <- 0
        current_posterior <- lik_t / lik_sum
      } else {
        current_posterior <- post / psi_t
      }
    } else {
      # Only update with the likelihood if it has mass and does not zero out
      # the prediction (na.rm because the likelihood can still contain NAs).
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

      # Normalising constant from the unconstrained update (this is the
      # quantity that goes into the NLL).
      psi_t <- sum(post, na.rm = TRUE)

      # Reachability bridge: drop cells from which the next known location
      # cannot be reached in the remaining steps.
      if (use_reach && !is.na(next_known[t])) {
        k <- next_known[t]
        allowed <- reach_dist[[as.character(k)]] <= (k - t) * reach_r
        post_r <- post
        post_r[!allowed] <- 0
        if (isTRUE(sum(post_r, na.rm = TRUE) > 0)) {
          post <- post_r
        } else {
          unreachable_vec[t - 1] <- TRUE
        }
      }

      current_posterior <- post / sum(post, na.rm = TRUE)
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

  # initial state is the first likelihood surface; psi (normalising constant
  # per step) is undefined at t = 1.
  #
  # Steps where the prediction does not overlap the likelihood fall back to
  # post <- pred_t, so psi_t ~ 1 and log(psi) ~ 0: a free pass. A too-small
  # sigma can therefore look better than it should. Known-event dates are
  # exempt from this: they are always applied and penalised if unreachable.
  # Also compare n_updated across sigmas, or use neg_log_lik / n_updated.

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
    n_updated = n_updated,
    n_known = sum(known_idx > 1),
    n_known_missed = sum(known_missed_vec),
    n_unreachable = sum(unreachable_vec),
    infeasible = sum(unreachable_vec) > 0 || sum(known_missed_vec) > 0,
    nll_per_update = if (n_updated > 0) neg_log_lik / n_updated else NA_real_,
    known_idx = known_idx
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

    results$kernel <- kernel
    results$likelihood_array <- likelihood_array

    results$normalised_posterior[land_arr] <- NA
    results$pc_predicted[land_arr] <- NA
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
#' @param forward_filter_obj A `list` from [hmm_forward_filter()] containing
#' `normalised_posterior`, `pc_predicted`, `likelihood_array` and `kernel`
#' (so `return_maps = TRUE`). If it also contains `known_idx`, those known
#' dates are honoured by the smoother.
#' @param progress a `logical` that controls progress bar.
#' @param known_events Optional `data.frame` of known-location dates. If
#' supplied it overrides `forward_filter_obj$known_idx`. Use the same data
#' frame given to [hmm_forward_filter()]. Default `NULL`.
#' @param date_col Name of the date column in `known_events`. Default `"date"`.
#' @details
#' All arrays use the layout `[x, y, t_step]`, with time as the last dimension.
#'
#' The final time step is seeded from the observation likelihood
#' `likelihood_array` rather than from the forward-filtered surface. For each
#' earlier time step from `n_t_steps - 1` down to 1, the function:
#' \enumerate{
#'   \item divides the smoothed surface at `t_step` by the predicted surface at `t_step`,
#'   \item convolves that ratio with the movement kernel `kernel`
#'     ([imagine::convolution2D()]), and
#'   \item multiplies the result by the forward-filtered surface at `t_step - 1`
#'     and normalises it to sum to 1.
#' }
#'
#' On known-event dates the smoothed surface is set to the normalised
#' likelihood at that date, overriding the recursion. This guarantees the
#' smoother is pinned to the known location even when the recursion collapses
#' (e.g. the kernel could not reach the known cell), and everything between
#' two known dates is smoothed conditional on both ends.
#'
#' Normalisation ignores `NA` cells (`na.rm = TRUE`), so land does not affect
#' the sums. Cells where `pc_predicted` is exactly 0 produce `NaN` or `Inf` in
#' the ratio; these are set to 0.
#'
#' @return A numeric array with the same dimensions as `normalised_posterior`,
#' `[x, y, t_steps]`, containing the smoothed probability surfaces.
#' Each time slice sums to 1 (ignoring `NA` cells).
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
#' @export
hmm_backward_smooth_new_2 <- function(
  forward_filter_obj,
  progress = TRUE,
  known_events = NULL,
  date_col = "date"
) {
  error_forward_filter_obj(forward_filter_obj)
  error_logical(progress)
  t0 <- msg_start("HMM backwards smoothing")

  # Filtered posteriors (alpha_t) and one-step predictions from the forward pass
  alpha <- forward_filter_obj$normalised_posterior
  pred <- forward_filter_obj$pc_predicted
  likelihood_array <- forward_filter_obj$likelihood_array
  kernel <- forward_filter_obj$kernel

  dims <- dim(alpha)
  dn <- dimnames(alpha)
  n_time <- dims[3]

  # Known-location dates: explicit argument wins, else what the forward
  # filter used.
  known_idx <- if (!is.null(known_events)) {
    known_event_steps(known_events, dn[[3]], date_col)
  } else if (!is.null(forward_filter_obj$known_idx)) {
    forward_filter_obj$known_idx
  } else {
    integer(0)
  }

  # Land mask: NA in the first slice. Land is static, so one slice is enough.
  land <- is.na(alpha[,, 1])

  smooth <- array(0, dim = dims, dimnames = dn)

  # Helper: normalised likelihood slice at step k (land = 0)
  known_slice <- function(k) {
    s <- likelihood_array[,, k]
    s[land] <- 0
    tot <- sum(s, na.rm = TRUE)
    if (is.finite(tot) && tot > 0) s / tot else NULL
  }

  # ----- Derive pop-off boundary condition from the last slice -----
  pop_off_slice <- known_slice(n_time)

  if (!is.null(pop_off_slice)) {
    smooth[,, n_time] <- pop_off_slice
  } else {
    # Fallback to the forward filter's last state if the slice has no mass
    smooth[,, n_time] <- alpha[,, n_time]
  }
  n_fallback <- 0L
  collapsed_dates <- character(0)
  time_labels <- dn[[3]]

  msg_progress(
    progress = progress,
    n_time = n_time,
    name = "HMM backwards smoothing",
    type = "Smoothing"
  )

  for (t in n_time:2) {
    # Ratio of smoothed to predicted. Land (NA), 0/0 and x/0 become 0.
    rat <- smooth[,, t] / pred[,, t]
    rat[!is.finite(rat)] <- 0

    # Propagate backwards (kernel assumed symmetric), then weight by the
    # filtered posterior at t - 1
    post <- imagine::convolution2D(rat, kernel) * alpha[,, t - 1]
    post[!land & !is.finite(post)] <- 0
    post[land] <- NA_real_

    fac <- sum(post, na.rm = TRUE)

    if (!is.finite(fac) || fac <= 0) {
      # Collapsed step: fall back to the filtered posterior for this step
      post <- alpha[,, t - 1]
      post[!land & !is.finite(post)] <- 0
      fac <- sum(post, na.rm = TRUE)
      n_fallback <- n_fallback + 1L
      collapsed_dates <- c(collapsed_dates, time_labels[t - 1])
      # error_smoothing_collapsed(facs = fac, time_step = t)
    }

    smooth[,, t - 1] <- post / fac

    # Known location at t - 1: pin the smoothed surface to the likelihood
    if ((t - 1L) %in% known_idx) {
      k_slice <- known_slice(t - 1L)
      if (!is.null(k_slice)) {
        smooth[,, t - 1] <- k_slice
      }
    }

    if (isTRUE(progress)) {
      cli::cli_progress_update()
    }
  }

  if (isTRUE(progress)) {
    cli::cli_progress_done()
  }

  # Land stays NA in every slice
  smooth[array(land, dim = dims)] <- NA

  msg_fallback(n_fallback = n_fallback, dates = rev(collapsed_dates))
  return(smooth)
}
