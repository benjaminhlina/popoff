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

#' @param dot_names an argument supplied to `...`.`
#' @name error_functions
#' @keywords internal

error_dots <- function(dot_names) {
  if (is.null(dot_names) || any(dot_names == "")) {
    cli::cli_abort(
      "All arguments passed through `...` must be named.",
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
#' @name error_functions
#' @keywords internal
error_extent <- function(vec, arg_name = NULL) {
  if (is.null(arg_name)) {
    arg_name <- rlang::as_label(rlang::enexpr(vec))
  }
  if (!is.numeric(vec) || length(vec) != 4 || anyNA(vec)) {
    cli::cli_abort(c(
      "{.arg {arg_name}} must be a numeric vector of length 4.",
      "x" = "You supplied {.cls {class(vec)}} of length {length(vec)}.",
      "i" = "Use the form {.code {arg_name} = c(xmin, xmax, ymin, ymax)} for the longtidue and latitude extents."
    ))
  }
}

#' @param array an array to check
#' @param name name of array
#' @name error_functions
#' @keywords internal
error_filter_array <- function(array, name) {
  dims <- dim(array)

  if (length(dims) != 3) {
    cli::cli_abort(
      c(
        "{.field {name}} must be a 3-dimensional array with dimensions {.code [x, y, t_steps]}.",
        "x" = "It has {length(dims)} dimension{?s}."
      )
    )
  }

  if (dims[3] < 2) {
    cli::cli_abort(
      c(
        "{.field {name}} must have at least 2 time steps.",
        "x" = "The third dimension (time) has length {dims[3]}."
      )
    )
  }

  invisible(dims)
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
    error_matrix(obj, dim_check = c(1, 2))
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
#' @param forward_filter_obj a `list` produced by [hmm_forward_filter()]
#' @name error_functions
#' @keywords internal
error_forward_filter_obj <- function(forward_filter_obj) {
  error_list(forward_filter_obj)
  required <- c(
    "normalised_posterior",
    "pc_predicted",
    "kernel",
    "likelihood_array"
  )
  missing <- setdiff(required, names(forward_filter_obj))

  if (length(missing) > 0) {
    cli::cli_abort(
      c(
        "{.arg forward_filter_obj} is missing required array{?s}: {.field {missing}}.",
        "i" = "Did you run {.fn hmm_forward_filter} without {.code return_maps = TRUE}?",
        "i" = "Re-run {.fn hmm_forward_filter} with {.code return_maps = TRUE}."
      )
    )
  }
  for (nm in required) {
    error_filter_array(forward_filter_obj[[nm]], nm)
  }
}

#' @param input_data a `data.frame` to assess
#' @param type the type of input data. A `character` usually `"raw_log"` or `"event_log"`.
#' @name error_functions
#' @keywords internal
error_input_data <- function(input_data, type) {
  error_df(input_data)

  if (type == "raw_log") {
    error_missing_names(
      obj = input_data,
      required_val = c("fish_id", "date", "date_time", "depth")
    )
    # --- Type and Value Validations ---

    error_chr(input_data$fish_id, col = "fish_id")
    error_date(input_data$date, col = "date", expected_cls = "Date")
    error_date(input_data$date_time, col = "date_time", expected_cls = "POSIXt")
    error_numeric_min(input_data$depth, col = "depth", min_val = 0)
  }
  if (type == "event_log") {
    error_missing_names(
      obj = input_data,
      required_val = c("fish_id", "date", "event")
    )
    error_chr(input_data$fish_id, col = "fish_id")
    error_chr(input_data$event, col = "event")
    error_date(input_data$date, col = "date", expected_cls = "POSIXt")
  }
}


#' @param list a `list`.
#' @name error_functions
#' @keywords internal

error_list <- function(list) {
  obj_name <- rlang::as_label(rlang::enexpr(list))
  if (!is.list(list)) {
    cli::cli_abort(
      c(
        "{.arg {obj_name}} must be a {.cls list}, not {.cls {class(list)}}."
      ),
      class = "error_list"
    )
  }
}

#' @param idx index value
#' @name error_functions
#' @keywords internal
error_likelihood_matrix <- function(matrix, idx) {
  error_matrix(matrix)

  n <- nrow(matrix)
  bad_idx <- idx < 1 | idx > n

  if (any(bad_idx)) {
    cli::cli_abort(
      c(
        "x" = "{.arg idx} has {sum(bad_idx)} index/indices out of bounds for a matrix with {n} row{?s}.",
        "i" = "Valid range is {.val 1} to {.val {n}}.",
        "i" = "First offending index: {.val {idx[which(bad_idx)[1]]}}"
      )
    )
  }

  bad <- matrix[, 3] > matrix[, 4]

  if (any(bad, na.rm = TRUE)) {
    cli::cli_abort(c(
      "x" = "{sum(bad, na.rm = TRUE)} row{?s} have {.code low > high} in depth integration bounds.",
      "i" = "Check for reversed or malformed {.field Depth_min}/{.field Depth_max} upstream.",
      "i" = "First offending row (of {idx}): {.val {which(bad)[1]}}"
    ))
  }
}
#' @name error_functions
#' @keywords internal
error_limits <- function(vec, arg_name = NULL) {
  if (is.null(arg_name)) {
    arg_name <- rlang::as_label(rlang::enexpr(vec))
  }
  if (!is.numeric(vec) || length(vec) != 2 || anyNA(vec)) {
    cli::cli_abort(c(
      "{.arg {arg_name}} must be a numeric vector of length 2.",
      "x" = "You supplied {.cls {class(vec)}} of length {length(vec)}.",
      "i" = "Use the form {.code {arg_name} = c(min, max)} for the longtidue and latitude extents."
    ))
  }
}

#' @param array a array
#' @name error_functions
#' @keywords internal

error_likelihood_array <- function(array, arg_name = NULL) {
  if (is.null(arg_name)) {
    arg_name <- rlang::as_label(rlang::enexpr(array))
  }

  if (!inherits(array, "array")) {
    cli::cli_abort(
      c(
        "{arg. arg_name} requires a 3-dimensional array.",
        "x" = "Please provide a 3-dimensional array."
      )
    )
  }
  dims <- dim(array)

  if (!length(dims) %in% 3) {
    cli::cli_abort(
      c(
        "{.arg arg_name} does not have the correct dimensions",
        "x" = "Please provide a likelihood array that consists of 3-dmensions that are [lon, lat, time_interval]"
      )
    )
  }

  if (dims[3] < 2) {
    cli::cli_abort("{.arg arg_name} must have at least 2 time slices.")
  }
}


#' @param linear_index the linear index column of event_log `data.frame` to check
#' @param event_date a `vector` containing event dates in event_log
#' @name error_functions
#' @keywords internal
error_linear_index <- function(event_date, linear_index) {
  if (anyNA(linear_index)) {
    cli::cli_abort(
      "Known location{?s} on {.val {event_date[is.na(linear_index)]}} fall outside the bathymetry."
    )
  }
}

#' Validate a prepped list input
#'
#' Checks that `x` is a list with the structure expected from either
#' `prep_summary_likdepth()` (`type = "log"`) or `prep_bathy_raster()`
#' (`type = "bathy"`), and throws an informative error if not.
#'
#' @param list The object to check.
#' @param type One of `"log"` or `"bathy"`, indicating which upstream
#'   function `x` is expected to have come from.
#'
#' @return Invisibly returns `TRUE` if `x` passes validation; otherwise
#'   aborts with an informative error.
#' @keywords internal
#'
#'
error_list_prep <- function(
  list,
  type = c("log", "bathy")
) {
  type <- rlang::arg_match(type)
  obj_name <- rlang::as_label(rlang::enexpr(list))

  # ---- generic "not a list at all" check --------------------------------
  error_list(list = list)

  # ---- type-specific structural checks -----------------------------------
  if (type == "log") {
    error_missing_names(
      obj = list,
      required_val = c(
        "raw_log_filter",
        "summary_table",
        "depth_profile",
        "time_res",
        "bin"
      )
    )

    summary_ok <- .has_valid_field(list, "summary_table", "data.frame") &&
      all(c("date", "depth_min", "depth_max") %in% names(list$summary_table))

    if (!summary_ok) {
      cli::cli_abort(
        c(
          "{.arg {obj_name}} doesn't look like the output of {.fn prep_summary_likdepth}.",
          "x" = "{.field summary_table} must be a data frame with columns {.field date}, {.field depth_min}, and {.field depth_max}.",
          "i" = "Run {.arg {obj_name}} through {.fn prep_summary_likdepth} first."
        ),
        class = "error_list_log"
      )
    }
  } else if (type == "bathy") {
    error_missing_names(
      obj = list,
      required_val = c(
        "bathy",
        "sdi",
        "idx",
        "idx_na",
        "focal_dim"
      )
    )

    rast_ok <- all(vapply(
      c("bathy", "sdi"),
      function(field) .has_valid_field(list, field, "SpatRaster"),
      logical(1)
    ))

    idx_ok <- all(vapply(
      c("idx", "idx_na"),
      function(field) .has_valid_field(list, field, is.numeric),
      logical(1)
    ))

    if (!rast_ok || !idx_ok) {
      cli::cli_abort(
        c(
          "{.arg {obj_name}} doesn't look like the output of {.fn prep_bathy_raster}.",
          "x" = if (!rast_ok) {
            "{.field bathy} and {.field sdi} must both be {.cls SpatRaster} objects."
          },
          "x" = if (!idx_ok) {
            "{.field idx} and {.field idx_na} must both be numeric vectors."
          },
          "i" = "Run {.arg {obj_name}} through {.fn prep_bathy_raster} first."
        ),
        class = "error_list_bathy"
      )
    }
  }

  invisible(TRUE)
}

#' @param obj an object to assess
#' @name error_functions
#' @keywords internal

error_load <- function(obj, arg_name = NULL) {
  if (is.null(arg_name)) {
    arg_name <- rlang::as_label(rlang::enexpr(obj))
  }

  valid_classes <- c("SpatRaster", "RasterLayer", "matrix", "character")

  if (!inherits(obj, valid_classes)) {
    cli::cli_abort(
      c(
        "{.arg obj} must be a {.cls SpatRaster}, {.cls RasterLayer}, matrix, or a file path to one of those.",
        "x" = "You supplied an object of class {.cls {class(obj)}}.",
        "i" = "{.arg obj} had {.field typeof}: {.val {typeof(obj)}}."
      ),
      class = "bathy_invalid_input"
    )
  }
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

error_matrix <- function(matrix, dim_check = NULL, arg_name = NULL) {
  if (is.null(arg_name)) {
    arg_name <- rlang::as_label(rlang::enexpr(matrix))
  }

  if (!is.matrix(matrix)) {
    cli::cli_abort("{.arg {arg_name}} must be a matrix.")
  }

  if (!is.null(dim_check)) {
    dims <- dim(matrix)
    if (any(dim_check > length(dims)) || any(dim_check < 1)) {
      cli::cli_abort(
        "{.arg dim_check} contains index/indices out of range for a matrix with {length(dims)} dimension(s)."
      )
    }
    checked_dims <- dims[dim_check]

    if (any(checked_dims %% 2 == 0)) {
      cli::cli_abort(
        c(
          "The supplied matrix in {.field {arg_name}} must have odd dimensions",
          "x" = "You supplied a matrix with dimensions: {.val {checked_dims}}"
        )
      )
    }
  }

  return(matrix)
}

#' @param mask a `logical`.
#' @param sigma a `numeric`.
#' @name error_functions
#' @keywords internal

error_mask <- function(mask, sigma) {
  if (isTRUE(mask) & sigma == 0) {
    cli::cli_abort(
      c(
        "When mask is {.val {mask}} and sigma is {.val {sigma}} Gaussian Kernel produces {.val {NaN}}",
        "x" = "If using a value of {.val {mask}} for mask then sigma has to be > {.val 0}"
      )
    )
  }
}
#' @param event_log_dates a `vector` containing event dates
#' @param summary_dates a `vector` containing sumamry dates
#' @name error_functions
#' @keywords internal

error_missing_dates <- function(event_log_dates, summary_dates) {
  missing_dates <- setdiff(event_log_dates, summary_dates)
  if (length(missing_dates) > 0) {
    cli::cli_abort(
      "Known date{?s} not in the likelihood array: {.val {missing_dates}}."
    )
  }
}


#' @param required_val a `vector` containing `characters` of the columns to check.
#' @name error_functions
#' @keywords internal

error_missing_names <- function(obj, required_val) {
  obj_name <- rlang::as_label(rlang::enquo(obj))
  missing_names <- setdiff(required_val, names(obj))

  if (is.list(obj)) {
    obj_type <- "list object"
  }
  if (is.data.frame(obj)) {
    obj_type <- "column"
  }
  if (length(missing_names) > 0) {
    cli::cli_abort(
      c(
        "{.var {obj_name}} is missing required {obj_type}{?s}.",
        "x" = "Missing {obj_type}{?s}: {.val {missing_names}}",
        "i" = "{.var {obj_name}} must contain: {.val {required_val}}"
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
        "i" = "Please provide the correct path to file you are wanting to use."
      )
    )
  }
}

#' @param predicted a 3-dimensional array of predicted values
#' @param likelihood a 3-dimensional array of likelihood values
#' @name error_functions
#' @keywords internal
error_pred <- function(predicted, likelihood) {
  pred_dims <- dim(predicted)
  lik_dim <- dim(likelihood)

  if (!(identical(pred_dims, lik_dim))) {
    cli::cli_abort(
      "The dimensions of predicted are {.val {pred_dims}} and do not match {.val {lik_dims}}"
    )
  }
}

#' @param out_dir a chracter string containing the out directory
#' @name error_functions
#' @keywords internal

error_out_dir <- function(out_dir) {
  if (is.null(out_dir)) {
    return(invisible(NULL))
  }

  if (!is.character(out_dir) || length(out_dir) != 1) {
    cli::cli_abort(
      c(
        "Invalid `out_dir` input.",
        "x" = "`out_dir` must be a single string, not {.obj_type_friendly {out_dir}}."
      )
    )
  }

  if (!dir.exists(out_dir)) {
    cli::cli_abort(
      c(
        "Output directory not found.",
        "x" = "No such directory: {.file {out_dir}}",
        "i" = "Please provide a valid, existing directory path."
      )
    )
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
#' @param back_transform a `logical` value
#' @param parameter_range a `vector` containing two `numerical` values.
#' @name error_functions
#' @keywords internal

error_transform <- function(back_transform, parameter_range) {
  if (isTRUE(back_transform) & is.null(parameter_range)) {
    cli::cli_abort(
      c(
        "When {.arg back_transform} is {.val TRUE}, {.arg parameter_range} 
        cannot be {.val NULL}",
        "x" = "Please supply {.arg parameter_range} with a {.cls vector} 
        containing two numerical values (e.g., {.val c(1, 5)})"
      )
    )
  }
}
#' @param workers a `numeric` that contains the number of works
#' @name error_functions
#' @keywords internal
error_workers <- function(workers) {
  if (is.null(workers) || workers < 2) {
    cli::cli_abort(
      "Parallel execution requested, but no mirai daemons are set. Please provide {.arg workers} (> 1) or call {.fn mirai::daemons} first."
    )
  }
}

#' @param unknown a character that is an argument.
#' @name error_functions
#' @keywords internal
error_wrong_arg <- function(unknown) {
  if (length(unknown) > 0) {
    cli::cli_abort(
      "Unused argument(s) in `...`: {.arg unknown}. See `gaussian_kernel()` and `hmm_forward_filter()` for valid arguments.",
    )
  }
}
