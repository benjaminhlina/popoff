#' Convert functions
#'
#' Convert an object of one class to another (e.g., `SpatRaster` to `matrix`
#'
#' @param x a `SpatRaster` object
#'
#' @name convert_transformers
#' @keywords internal

.convert_spr_to_vec <- function(x) {
  error_spatrast(x)
  # if (inherits(x, "SpatRaster")) {
  convert <- terra::as.matrix(x, wide = TRUE) |>
    as.vector()
  # wide = TRUE returns matrix matching spatial grid [rows, cols]
  # as.vector() on matrix converts it in standard column-major order
  return(convert)
  # }
  # return(as.vector(x))
}

#' Check whether a list element exists and satisfies a type check
#'
#' @param x A list.
#' @param field Name of the element to check.
#' @param check Either a class string (checked via `inherits()`) or a
#'   predicate function (e.g. `is.numeric`) applied to the element.
#'
#' @return `TRUE`/`FALSE`.
#' @keywords internal
.has_valid_field <- function(x, field, check) {
  if (!field %in% names(x)) {
    return(FALSE)
  }

  val <- x[[field]]

  if (is.character(check)) {
    inherits(val, check)
  } else if (is.function(check)) {
    isTRUE(check(val))
  } else {
    rlang::abort(
      "`check` must be a class string or a predicate function.",
      call = rlang::caller_env()
    )
  }
}

#' Integrate dnorm
#'
#' @param mean the mean value
#' @param sd the sd
#' @param lower lower bounds
#' @param upper upper bounds
#'
#' @name intergration_helper
#' @keywords internal

.integrate_dnorm_rows <- function(mean, sd, lower, upper) {
  integrated <- pmax(
    0,
    stats::pnorm(q = upper, mean = mean, sd = sd) -
      stats::pnorm(q = lower, mean = mean, sd = sd)
  )
  return(integrated)
}
#' Safe base functions
#'
#' These functions check if x is `NA` then
#' return a stat (e.g, range, min, max)
#'
#' @param x a `numerical` vector
#'
#' @name safe_functions
#' @keywords internal
.safe_max <- function(x) {
  if (all(is.na(x))) {
    return(NA_real_)
  }
  max(x, na.rm = TRUE)
}

#' @param type a `character` determining whether `min()` or `max()` is used
#'
#' @name safe_functions
#' @keywords internal
.safe_range_val <- function(x, type = c("min", "max")) {
  x <- x[!is.na(x)]
  if (length(x) == 0) {
    return(NA_real_)
  }
  if (type == "min") min(x) else max(x)
}

#' CRS validation
#'
#' Checks if a valid crs is provided
#'
#' @param crs a coordinate reference system
#' @name crs_functions
#' @keywords internal
.valid_crs <- function(crs) {
  if (!inherits(crs, "character") || length(crs) != 1 || is.na(crs)) {
    return(FALSE)
  }

  result <- tryCatch(
    terra::crs(crs),
    error = function(e) NA_character_,
    warning = function(w) NA_character_
  )
  valid <- !is.na(result) && nzchar(result)

  return(valid)
}


#' Index functions
#'
#' Provides an index number.
#'
#' @param x an object to get a index
#' @param na a `TRUE/FALSE` that controls whether index locations that are `NA` are returned. Default is `FALSE`.`
#' @name index_functions
#' @keywords internal

.valid_idx <- function(x, na = FALSE) {
  vals <- .convert_spr_to_vec(x)

  if (isFALSE(na)) {
    idx <- which(!is.na(vals))
  } else {
    idx <- which(is.na(vals))
  }

  return(idx)
}
