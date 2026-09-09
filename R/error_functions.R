#' Error functions
#'
#' These functions produce error message displayed to the user.
#' @param chr
#'

#' @param crs a character string of a crs
#' @param arg_name The name of the argument
#'
#' @name error_functions
#' @keywords internal

error_crs <- function(crs, arg_name = NULL) {
  if (is.null(arg_name)) {
    arg_name <- rlang::as_label(rlang::enexpr(crs))
  }

  if (!.valid_crs(crs)) {
    cli::cli_abort(
      c(
        "{arg_name} is not a valid coordinate reference system (crs).",
        "x" = "Please provide {arg_name} with a proper crs, e.g. {.val EPSG:4326}"
      )
    )
  }
}


#' @param df a `data.frame`

#' @name error_functions
#' @keywords internal

error_df <- function(df, arg_name = NULL) {
  if (is.null(arg_name)) {
    arg_name <- rlang::as_label(rlang::enexpr(df))
  }
  if (!(inherits(df, "data.frame"))) {
    cli::cli_abort(
      c(
        "{arg_name} is not {.class data.frame}",
        "x" = "Please provide {arg_name} with a {.class data.frame}"
      )
    )
  }
}


#' @param path The file path to a file.
#' @name error_functions
#' @keywords internal

error_empty_file <- function(path) {
  cli::cli_abort(
    c(
      "The file did not contain any objects.",
      "x" = "File path: {.file {path}}"
    ),
  )
}

#' @param obj an object either `vector` or `matrix`
#' @name error_functions
#' @keywords internal
error_focal_dim <- function(obj, arg_name = NULL) {
  if (is.null(arg_name)) {
    arg_name <- rlang::as_label(rlang::enexpr(obj))
  }

  if (inherits(obj, "matrix")) {
    if (!is.numeric(obj)) {
      error_numeric(obj)
    }
    error_matrix(obj)
  } else {
    error_numeric(obj)
    if (any(obj %% 2 == 0)) {
      cli::cli_abort(
        c(
          "Supplied value in {arg_name} has to be odd",
          "x" = "You supplied: {.val {obj}}"
        )
      )
    }
  }
}

#' @param e the condition object caught by `tryCatch()`
#' @name error_functions
#' @keywords internal
error_load <- function(path, e) {
  cli::cli_abort(
    c(
      "Failed to read {.arg path} as a raster.",
      "x" = "{.file {path}} could not be loaded by {.fn terra::rast}.",
      "i" = "{.arg path} was type {.cls {class(path)}}.",
      "i" = "Original error: {conditionMessage(e)}"
    ),
    class = "bathy_rast_read_failed"
  )
}

#' @param vec a `vector` to check if numeric
#' @name error_functions
#' @keywords internal
error_logical <- function(vec, arg_name = NULL) {
  if (is.null(arg_name)) {
    arg_name <- rlang::as_label(rlang::enexpr(vec))
  }

  if (!(inherits(vec, "logical"))) {
    cli::cli_abort(
      c(
        "{arg_name} needs to be a logical value",
        "x" = "Please supply {arg_name} as a {.cls logic}"
      )
    )
  }
}
#' @param matrix a `matrix` to check if valid
#' @name error_functions
#' @keywords internal

error_matrix <- function(matrix, arg_name = NULL) {
  if (is.null(arg_name)) {
    arg_name <- rlang::as_label(rlang::enexpr(matrix))
  }
  dims <- dim(matrix)

  if (any(dims %% 2 == 0)) {
    cli::cli_abort(
      c(
        "Supplied matrix in {arg_name} must have odd dimensions",
        "x" = "You supplied a matrix with dimensions: {.val {dims[1]}} x {.val {dims[2]}}"
      )
    )
  }
}

#' @name error_functions
#' @keywords internal

error_numeric <- function(vec, arg_name = NULL) {
  if (is.null(arg_name)) {
    arg_name <- rlang::as_label(rlang::enexpr(vec))
  }

  if (!(inherits(vec, "numeric"))) {
    cli::cli_abort(
      c(
        "{arg_name} needs to be a numeric",
        "x" = "Please supply {arg_name} as a {.cls numeric}"
      )
    )
  }
}

#' @name error_functions
#' @keywords internal

error_path <- function(path) {
  if (!file.exists(path)) {
    cli::cli_abort(
      c(
        "File not found.",
        "x" = "No such file: {.file {path}}",
        "i" = "Please provide the correct path to file you are wanting to use "
      ),
    )
  }
}

#' @param raster a raster object to assess
#' @name error_functions
#' @keywords internal

error_raster <- function(raster) {
  cli::cli_abort(
    c(
      "{.arg raster} must be a {.cls SpatRaster}, {.cls RasterLayer}, matrix, or a file path to one of those.",
      "x" = "You supplied an object of class {.cls {class(raster)}}.",
      "i" = "{.arg raster} had {.field typeof}: {.val {typeof(raster)}}."
    ),
    class = "bathy_invalid_input"
  )
}

#' @param spatrast a raster object to assess
#' @name error_functions
#' @keywords internal

error_spatrast <- function(spatrast, arg_name = NULL) {
  if (is.null(arg_name)) {
    arg_name <- rlang::as_label(rlang::enexpr(spatrast))
  }

  if (!(inherits(spatrast, "SpatRaster"))) {
    cli::cli_abort(
      c(
        "{arg_name} is not {.class SpatRaster}",
        "x" = "Please provide {arg_name} with a {.class SpatRaster}"
      )
    )
  }
}
