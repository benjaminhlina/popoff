#' Likelihood integration functions
#'
#' @param w a `vector` containing the mean of an area
#' @param wsd a `vector` containing the standard deviation
#' @param minT a `vector` containing the minmum temperature
#' @param maxT a `vector` containing the maximum temperature
#'
#' @name  likelihood_integration
#' @keywords internal

likint3 <- function(w, wsd, minT, maxT) {
  widx <- !is.na(w)
  wdf <- data.frame(w = as.vector(w[widx]), wsd = as.vector(wsd[widx]))
  wint <- .integrate_dnorm_rows(w[widx], wsd[widx], minT, maxT)
  w <- w * 0
  w[widx] <- wint
  w
}

#' Integrate N(mean, sd) between the minimum and maxium depth for a subset of cells
#'
#' Vectorized replacement for row-wise stats::integrate(dnorm, ...) --
#' exact, not an approximation: integral of dnorm from low to high is
#' pnorm(high) - pnorm(low) by definition of the CDF.
#'
#' @param w a `matrix`/`data.frame` with columns: mean, sd, lower bound, upper bound
#' @param widx a logical vector indicating which rows of `w` to integrate - an index
#'
#' @details `likint3_matrix` use to be called `linkint3.sub`
#' @name  likelihood_integration
#' @keywords internal
#'
likint3_matrix <- function(w, widx) {
  error_likelihood_matrix(matrix = w, idx = widx)
  wdf <- w[widx, , drop = FALSE]

  wint <- .integrate_dnorm_rows(
    mean = wdf[, 1],
    sd = wdf[, 2],
    lower = wdf[, 3],
    upper = wdf[, 4]
  )
  w[, 1] <- 0
  w[widx, 1] <- wint
  w[, 1]
}
