#' Error functions
#'
#' These functions produce error message displayed to the user.
#' @param chr a character to check
#' @param col a character containing the column name
#'
#' @name error_functions
#' @keywords internal
error_chr <- function(chr, col = NULL) {
  if (!is.character(chr)) {
    cli::cli_abort(
      c(
        "Invalid column type in {.field {col}}.",
        "x" = "{.field {col}} must be a {.cls character} vector, not {.obj_type_friendly {chr}}."
      )
    )
  }
}

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


#' @param date a date to check
#' @param expected_cls a class to be expected e.g., `Date` or `POSIXct`
#'
#' @name error_functions
#' @keywords internal

error_date <- function(date, col, expected_cls) {
  if (!inherits(date, expected_cls)) {
    cli::cli_abort(
      c(
        "Invalid column type for {.field {col}}.",
        "x" = "{.field {col}} must inherit from {.cls {expected_cls}}, not {.obj_type_friendly {date}}."
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
#' @param required_col a `vector` containing `characters` of the columns to check.
#' @name error_functions
#' @keywords internal

error_missing_col <- function(df, required_col) {
  df_name <- rlang::as_label(rlang::enquo(df))
  missing_cols <- setdiff(required_cols, names(df))

  if (length(missing_cols) > 0) {
    cli::cli_abort(
      c(
        "{.var {df_name}} is missing required column{?s}.",
        "x" = "Missing column{?s}: {.val {missing_cols}}",
        "i" = "{.var {df_name}} must contain: {.val {required_cols}}"
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

#' @param min_val a numerical mimimum value allowed
#' @name error_functions
#' @keywords interna

error_numeric_min <- function(vec, col, min_val = 0) {
  error_numeric(vec)

  if (any(vec <= min_val, na.rm = TRUE)) {
    cli::cli_abort(
      c(
        "Invalid values in {.field {col}} column.",
        "x" = "{.field {col}} must contain values strictly greater than {min_val}.",
        "i" = "Found values <= {min_val}."
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

#' @param output_dir a chracter string containing the out directory
#' @name error_functions
#' @keywords internal

error_out_dir <- function(output_dir) {
  if (is.null(output_dir)) {
    return(invisible(NULL))
  }

  if (!is.character(output_dir) || length(output_dir) != 1) {
    cli::cli_abort(
      c(
        "Invalid `output_dir` input.",
        "x" = "`output_dir` must be a single string, not {.obj_type_friendly {output_dir}}."
      )
    )
  }

  if (!dir.exists(output_dir)) {
    cli::cli_abort(
      c(
        "Output directory not found.",
        "x" = "No such directory: {.file {output_dir}}",
        "i" = "Please provide a valid, existing directory path."
      )
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

#' @param input_data a `data.frame` to assess
#' @param type the type of input data. A `character` usually `"raw_log"` or `"iniloc"`.
#' @name error_functions
#' @keywords internal
error_input_data <- function(input_data, type) {
  error_df(input_data)

  if (type == "raw_log") {
    error_missing_col(
      df = input_data,
      required_col = c("fish_id", "date", "date_time", "depth")
    )
    # --- Type and Value Validations ---

    error_chr(input_data$fish_id, col = "fish_id")
    error_date(input_data$date, col = "date", expected_cls = "Date")
    error_date(input_data$date_time, col = "date_time", expected_cls = "POSIXt")
    error_numeric_min(input_data$depth, col = "depth", min_val = 0)
  }
  if (type == "iniloc") {
    error_missing_col(
      df = input_data,
      required_col = c("fish_id", "date", "event")
    )
    error_chr(input_data$fish_id, col = "fish_id")
    error_chr(input_data$event, col = "event")
    error_date(input_data$date, col = "date", expected_cls = "POSIXt")
  }
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
