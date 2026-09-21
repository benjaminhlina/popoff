test_that("error_list_prep errors on a valid log list [KNOWN BUG: undefined `x`]", {
  local_mocked_bindings(
    error_list = function(list) invisible(TRUE),
    error_missing_names = function(obj, required_val) invisible(TRUE)
  )

  valid_log <- list(
    raw_log_filter = TRUE,
    summary_table = data.frame(date = 1, depth_min = 1, depth_max = 2),
    depth_profile = TRUE,
    time_res = TRUE,
    bin = TRUE
  )

  # KNOWN BUG: `summary_ok` reads `x$summary_table` instead of
  # `list$summary_table`. `x` doesn't exist in this scope, so even a
  # perfectly valid log list errors instead of returning invisible(TRUE).
  expect_true(
    error_list_prep(valid_log, type = "log")
  )
})

test_that("error_list_prep does not error on bathy list [KNOWN BUG: undefined `x`]", {
  local_mocked_bindings(
    error_list = function(list) invisible(TRUE),
    error_missing_names = function(obj, required_val) invisible(TRUE),
    .has_valid_field = function(x, field, check) TRUE
  )

  valid_bathy <- list(
    bathy = terra::rast(),
    sdi = terra::rast(),
    idx = 1:5,
    idx_na = NA_integer_,
    focal_dim = 3,
    weighted = TRUE
  )
  expect_true(
    error_list_prep(valid_bathy, type = "bathy")
  )
})

test_that("error_list_prep propagates error_list()'s failure for non-list input", {
  local_mocked_bindings(
    error_list = function(list) {
      cli::cli_abort("not a list", class = "error_list")
    }
  )

  expect_error(
    error_list_prep("not a list", type = "log"),
    class = "error_list"
  )
})

test_that("error_list_prep propagates error_missing_names()'s failure", {
  local_mocked_bindings(
    error_list = function(list) invisible(TRUE),
    error_missing_names = function(obj, required_val) {
      cli::cli_abort("missing names", class = "error_missing_names")
    }
  )

  expect_error(
    error_list_prep(list(a = 1), type = "log"),
    class = "error_missing_names"
  )
})

test_that("error_list_prep requires `type` to be 'log' or 'bathy'", {
  expect_error(
    error_list_prep(list(a = 1), type = "nonsense")
  )
})

test_that("error_list_prep's `call = call` in the log branch is an undefined-variable bug [KNOWN BUG, currently masked]", {
  # This is unreachable today: to hit cli_abort() in the "log" branch,
  # `summary_ok` must be FALSE, but *evaluating* `summary_ok` already errors
  # due to the undefined `x` bug above. Once that's fixed, this test should
  # be revisited -- cli_abort(call = call) will itself fail with
  # "object 'call' not found" since `call` is never assigned anywhere in
  # error_list_prep()'s body.
  local_mocked_bindings(
    error_list = function(list) invisible(TRUE),
    error_missing_names = function(obj, required_val) invisible(TRUE)
  )

  bad_log <- list(
    raw_log_filter = TRUE,
    summary_table = data.frame(x = 1), # missing required columns
    depth_profile = TRUE,
    time_res = TRUE,
    bin = TRUE
  )

  expect_error(
    error_list_prep(bad_log, type = "log"),
    "doesn't look like the output of",
    fixed = TRUE
  )
})

describe("error_list_prep() - bathy branch structural failures", {
  it("aborts with class error_list_bathy when bathy/sdi are not SpatRaster", {
    testthat::local_mocked_bindings(
      error_list = function(list) invisible(TRUE),
      error_missing_names = function(obj, required_val) invisible(TRUE),
      .has_valid_field = function(obj, field, type) {
        # rast fields fail, idx fields pass
        !(field %in% c("bathy", "sdi"))
      }
    )

    dummy_list <- list(
      bathy = "not_a_raster",
      sdi = "not_a_raster",
      idx = 1:5,
      idx_na = NA_integer_,
      focal_dim = 3,
      weighted = TRUE
    )

    err <- rlang::catch_cnd(
      error_list_prep(list = dummy_list, type = "bathy")
    )

    expect_s3_class(err, "error_list_bathy")
    msg <- conditionMessage(err)

    expect_true(grepl("doesn't look like the output of", msg, fixed = TRUE))
    expect_true(grepl("bathy", msg, fixed = TRUE))
    expect_true(grepl("sdi", msg, fixed = TRUE))
    expect_true(grepl("must both be", msg, fixed = TRUE))
    expect_true(grepl("SpatRaster", msg, fixed = TRUE))
    # idx branch message should NOT appear since idx_ok is TRUE
    expect_false(grepl("numeric vectors", msg, fixed = TRUE))
  })

  it("aborts with class error_list_bathy when idx/idx_na are not numeric", {
    testthat::local_mocked_bindings(
      error_list = function(list) invisible(TRUE),
      error_missing_names = function(obj, required_val) invisible(TRUE),
      .has_valid_field = function(obj, field, type) {
        # idx fields fail, rast fields pass
        !(field %in% c("idx", "idx_na"))
      }
    )

    dummy_list <- list(
      bathy = "raster_stub",
      sdi = "raster_stub",
      idx = "not_numeric",
      idx_na = "not_numeric",
      focal_dim = 3,
      weighted = TRUE
    )

    err <- rlang::catch_cnd(
      error_list_prep(list = dummy_list, type = "bathy")
    )

    expect_s3_class(err, "error_list_bathy")
    msg <- conditionMessage(err)

    expect_true(grepl("doesn't look like the output of", msg, fixed = TRUE))
    expect_true(grepl("idx", msg, fixed = TRUE))
    expect_true(grepl("idx_na", msg, fixed = TRUE))
    expect_true(grepl("numeric vectors", msg, fixed = TRUE))
    # rast branch message should NOT appear since rast_ok is TRUE
    expect_false(grepl("SpatRaster", msg, fixed = TRUE))
  })

  it("includes both failure messages when rast and idx checks both fail", {
    testthat::local_mocked_bindings(
      error_list = function(list) invisible(TRUE),
      error_missing_names = function(obj, required_val) invisible(TRUE),
      .has_valid_field = function(obj, field, type) FALSE
    )

    dummy_list <- list(
      bathy = "bad",
      sdi = "bad",
      idx = "bad",
      idx_na = "bad",
      focal_dim = 3,
      weighted = TRUE
    )

    err <- rlang::catch_cnd(
      error_list_prep(list = dummy_list, type = "bathy")
    )

    expect_s3_class(err, "error_list_bathy")
    msg <- conditionMessage(err)

    expect_true(grepl("SpatRaster", msg, fixed = TRUE))
    expect_true(grepl("numeric vectors", msg, fixed = TRUE))
  })
})
