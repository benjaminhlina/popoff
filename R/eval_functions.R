#' Evaluate sigma
#'
#' This function provides the ablity to easily give a value of sigma which
#' will first calcualte the Gaussian Kernel then use that kernel and
#' the likelihood values to use a hidden Markov model foward filter that
#' produces a negative log likelihood value that we can evaluate.
#'
#' @param sigma a `numerical` value to use in the evaluation.
#' @param normalize a `logical` value that changes whether the kernel is
#' normalized. See [gaussian_kernel()]. Defults to `FALSE`.
#' @param mask a `logical` that changes whether a maks is applied to the
#' kernel. See [gaussian_kernel()]. Defults to `TRUE`.
#' @param likelihood_arra the `likelihood_array` produced by
#' [calc_likdepth_transmit()].
#' @param grid a grid to evaluate and filter the foward filtering process
#'
#' @param ... additional named arguments passed to either
#' [gaussian_kernel()] (`muadv`, `back_transform`, `sigma_range`) or
#' [hmm_forward_filter()] (`bound_threshold`, `bound_min`,
#' `mask_likelihood`, `progress`). Arguments are routed by name. Unnamed or
#' unrecognized arguments cause an error.
#'
#'
#' @return a `list` containing the `sigma` value, the `neg_log_lik`,
#' `n_bad` which is the number of days that did not have a valid estimate,
#' and `n_updated` the number of times the the updated `likelihood` was used.
#'
#' @seealso [gaussian_kernel()], [hmm_forward_filter()]
#' @export
eval_sigma <- function(
  sigma,
  normalize = FALSE,
  mask = TRUE,
  likelihood_array,
  grid,
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
      names(formals(hmm_forward_filter)),
      c("likelihood_array", "kernel", "grid")
    )

    unknown <- setdiff(dot_names, c(kern_formals, hmm_formals))

    if (length(unknown) > 0) {
      cli::cli_abort(
        "Unused argument(s) in `...`: {.arg unknown}. See `gaussian_kernel()` and `hmm_forward_filter()` for valid arguments.",
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
  hmm_forward <- do.call(
    hmm_forward_filter,
    c(
      list(
        likelihood_array = likelihood_array,
        kernel = kern,
        grid = grid,
        return_maps = FALSE
      ),
      hmm_args
    )
  )

  out <- data.frame(
    sigma = sigma,
    log_lik = hmm_forward$log_lik,
    neg_log_lik = hmm_forward$neg_log_lik,
    n_bad = hmm_forward$n_bad,
    n_updated = hmm_forward$n_updated,
    n_steps = length(hmm_forward$ncps)
  )
  rm(hmm_forward, kern)
  gc(verbose = FALSE)
  return(out)
}

#' Evaluate the forward-filter negative log-likelihood over a vector of sigmas
#'
#' Evaluates `eval_sigma()` for each value in `sigma_vec` (in parallel if
#' mirai daemons are set) so the best-fitting sigma can be identified.
#'
#' @param sigma_vec A numeric vector of sigma values.
#' @param likelihood_array Likelihood array.
#' @param grid A `data.frame` containing the grid.
#' @param ... Further arguments passed to `eval_sigma()`.
#' @param progress Logical; show a progress bar over sigma values.
#'
#' @return A tibble with columns `sigma` and `neg_log_lik`.
#' @export
eval_sigma_parallel <- function(
  sigma_vec,
  likelihood_array,
  grid,
  ...,
  progress = TRUE
) {
  t0 <- msg_start("sigma evaluation in parallel")

  if (mirai::daemons_set()) {
    likelihoods_shared <- mori::share(likelihood_array)
    grid_shared <- mori::share(grid)

    res <- mirai::mirai_map(
      sigma_vec,
      function(eval_sigma, sigma_val, likelihoods_shared, grid_shared, ...) {
        eval_sigma(
          sigma_val,
          likelihoods_shared,
          grid_shared,
          ...,
          progress = FALSE
        )
      },
      .args = list(
        likelihoods = likelihoods_shared,
        grid = grid_shared,
        eval_sigma = eval_sigma,
        ...
      )
    )

    if (isTRUE(progress)) {
      out <- res[.progress]
    }
    out <- res[]
  } else {
    out <- purrr::map(
      sigma_vec,
      eval_sigma,
      likelihood_array = likelihood_array,
      grid = grid,
      ...,
      progress = FALSE,
      .progress = if (progress) "Evaluating sigma" else FALSE
    )
  }

  out_df <- out |>
    dplyr::bind_rows()

  msg_end('Sigma evaluation in parallel', t0 = t0)
  return(out_df)
}
