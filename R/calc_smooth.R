smooth_hmm <- function(
  sigma,
  filter_obj,
  normalize = FALSE,
  mask = TRUE,
  likelihood_array,
  ...
) {
  dots <- list(...)

  # route extra arguments by name ----
  if (length(dots) > 0) {
    dot_names <- names(dots)

    if (is.null(dot_names) || any(dot_names == "")) {
      cli::cli_abort(
        "All arguments passed through `...` must be named.",
      )
    }

    kern_formals <- setdiff(
      names(formals(gaussian_kernel)),
      c("sigma", "normalize", "mask")
    )
    hmm_formals <- setdiff(
      names(formals(hmm_backwards_smooth)),
      c("likelihood_array", "kernel", )
    )

    unknown <- setdiff(dot_names, c(kern_formals, hmm_formals))

    if (length(unknown) > 0) {
      cli::cli_abort(
        "Unused argument(s) in `...`: {.arg unknown}. See `gaussian_kernel()` and `hmm_backwards_smooth()` for valid arguments.",
      )
    }

    kern_args <- dots[dot_names %in% kern_formals]
    hmm_args <- dots[dot_names %in% hmm_formals]
  } else {
    kern_args <- list()
    hmm_args <- list()
  }

  # ------ run kernel ------
  kern <- do.call(
    gaussian_kernel,
    c(
      list(
        sigma = sigma,
        normalize = normalize,
        mask = mask
      ),
      kern_args
    )
  )

  # ----- run hmmm filter -----
  hmm_smooth <- do.call(
    hmm_backwards_smooth,
    c(
      list(
        likelihood_array = likelihood_array,
        kernel = kern,
        return_maps = return_maps
      ),
      hmm_args
    )
  )

  if (isFALSE(return_maps)) {
    out <- data.frame(
      sigma = sigma,
      log_lik = hmm_forward$log_lik,
      neg_log_lik = hmm_forward$neg_log_lik,
      n_bad_days = hmm_forward$n_bad_days,
      n_updated = hmm_forward$n_updated,
      n_steps = length(hmm_forward$ncps)
    )

    rm(hmm_forward, kern)
    gc(verbose = FALSE)
  } else {
    out <- hmm_forward
    rm(kern)
    gc(verbose = FALSE)
  }
  return(out)
}
