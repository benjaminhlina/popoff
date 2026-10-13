#' Refine an optimal sigma estimate by iterative grid narrowing
#'
#' Takes the results of an initial coarse sigma search and zooms in on the
#' best value over one or more rounds. Each round builds a finer grid of sigma
#' values around the current best estimate, evaluates the likelihood at each
#' value with [eval_sigma_iterate()], and keeps the result for the next round.
#'
#' @details
#' In each round the search window is the best sigma from the previous grid
#' plus or minus one grid step. The step is the spacing of the previous grid,
#' calculated as `(max(sigma) - min(sigma)) / (nrow(res) - 1)`. The window is
#' clipped to `[0, max_sigma]`, and `n_points` values are spaced evenly across
#' it. Because the window shrinks each round, the effective resolution
#' increases geometrically with `n_rounds`.
#'
#' `n_points` should be odd so that the previous best value falls on the new
#' grid and the refined estimate cannot be worse than the previous one. This
#' holds only when the window is not truncated by `0` or `max_sigma`. If the
#' best value sits at a boundary, the window is asymmetric and the previous
#' best may no longer be included.
#'
#' The returned optimum comes from the final round only. If `n_rounds = 0`, the
#' best row of `sigma_iterate_obj` is returned unchanged.
#'
#' @param sigma_iterate_obj A data frame of initial sigma search results, as
#'   returned by [eval_sigma_iterate()]. It must contain the columns `sigma`
#'   (the sigma values evaluated) and `neg_log_lik` (the negative
#'   log-likelihood at each value), and must have at least two rows so a grid
#'   step can be calculated.
#' @param max_sigma Numeric. Upper bound on sigma. Refined grids are never
#'   extended beyond this value (the lower bound is fixed at 0).
#' @param likelihood_array The likelihood array passed to
#'   [eval_sigma_iterate()] for evaluating each candidate sigma.
#' @param n_points Integer. Number of sigma values to evaluate in each
#'   refinement round. Use an odd number so that the previous best value is
#'   included in the new grid. Default is `17`.
#' @param n_rounds Integer. Number of refinement rounds. Default is `2`.
#' @param normalize Logical. Passed to [eval_sigma_iterate()]. Default is
#'   `FALSE`.
#' @param mask_likelihood Logical. Passed to [eval_sigma_iterate()]. Default
#'   is `FALSE`.
#' @param .progress Logical. Whether to show a progress bar during parallel
#'   evaluation. Passed to [eval_sigma_iterate()]. Default is `FALSE`.
#'
#' @return A one-row data frame with columns:
#' \describe{
#'   \item{`optimal_sigma`}{The sigma value with the lowest negative
#'     log-likelihood in the final refinement round.}
#'   \item{`neg_log_lik`}{The negative log-likelihood at `optimal_sigma`.}
#' }
#'
#' @seealso [eval_sigma_iterate()] for the underlying likelihood evaluation.
#'
#' @examples
#' \dontrun{
#' # Coarse initial search
#' coarse <- eval_sigma_iterate(
#'   sigma_vec = seq(0.5, 10, length.out = 20),
#'   likelihood_array = lik_array
#' )
#'
#' # Refine around the best coarse value
#' refine_sigma(
#'   sigma_iterate_obj = coarse,
#'   max_sigma = 10,
#'   likelihood_array = lik_array,
#'   n_points = 17,
#'   n_rounds = 3
#' )
#' }
#'
#' @export

refine_sigma <- function(
  sigma_iterate_obj,
  max_sigma,
  likelihood_array,
  n_points = 17, # odd, so the previous best is on the new grid
  n_rounds = 2,
  normalize = FALSE,
  mask_likelihood = FALSE,
  workers = NULL,
  parallel = TRUE,
  .progress = FALSE
) {
  res <- sigma_iterate_obj

  for (i in seq_len(n_rounds)) {
    best <- res$sigma[which.min(res$neg_log_lik)]
    step <- (max(res$sigma) - min(res$sigma)) / (nrow(res) - 1)

    res <- eval_sigma_iterate(
      sigma_vec = seq(
        max(0, best - step),
        min(max_sigma, best + step),
        length.out = n_points
      ),
      parallel = parallel,
      workers = workers,
      likelihood_array = likelihood_array,
      normalize = normalize,
      mask_likelihood = mask_likelihood,
      .progress = .progress
    )
  }

  top <- res[which.min(res$neg_log_lik), ]

  data.frame(
    optimal_sigma = top$sigma,
    neg_log_lik = top$neg_log_lik
  )
}
