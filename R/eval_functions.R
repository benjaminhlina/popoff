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
      c("likelihood_array", "kernal", "grid")
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
        grid = grid
      ),
      hmm_args
    )
  )

  out <- list(
    sigma = sigma,
    neg_log_lik = hmm_foward$neg_log_lik,
    n_bad = hmm_foward$n_bad,
    n_updated = hmm_foward$n_updated
  )
}
