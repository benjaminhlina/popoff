# ---------------------------------------------------------------------------
# Tests for calc_likdepth_transmit() and its daily helper calc_daily_likdepth()
# ---------------------------------------------------------------------------

# ---- shared fixtures --------------------------------------------------------

make_bathy_fixture <- function() {
  # 2 x 2 grid. One NA cell (bathy) so idx / idx_na are both exercised.
  bathy <- terra::rast(nrows = 2, ncols = 2, vals = c(10, 20, NA, 40))
  sdi <- terra::rast(nrows = 2, ncols = 2, vals = c(2, 2, 2, 2))

  # derive idx / idx_na from the *actual* vectorised order used internally
  # (.convert_spr_to_vec uses as.matrix(wide = TRUE) |> as.vector(), which is
  # not guaranteed to match terra::values() ordering) so the fixture stays
  # correct even if that internal ordering ever changes.
  zbi_vec <- .convert_spr_to_vec(bathy)
  idx <- which(!is.na(zbi_vec))
  idx_na <- which(is.na(zbi_vec))

  list(
    bathy = bathy,
    sdi = sdi,
    idx = idx,
    idx_na = idx_na,
    zbi_vec = zbi_vec,
    dims = dim(bathy)
  )
}

make_summary_table_2 <- function(dates, depth_min, depth_max) {
  data.frame(
    date = dates,
    depth_min = depth_min,
    depth_max = depth_max
  )
}

# =============================================================================
# calc_daily_likdepth() -- unit tests (no arg_match involved, testable direct)
# =============================================================================

test_that("calc_daily_likdepth returns an uninformative (all-cell) likelihood when there is no usable depth range", {
  fx <- make_bathy_fixture()
  sdi_dsd <- .convert_spr_to_vec(fx$sdi) + 0.57

  empty_log <- make_summary_table_2(as.Date("2020-01-01"), NA_real_, NA_real_)[
    0,
  ]

  out <- calc_daily_likdepth(
    summary_log_i = empty_log,
    zbi = fx$zbi_vec,
    sdi_dsd = sdi_dsd,
    idx = fx$idx,
    idx_na = fx$idx_na,
    dims = fx$dims,
    time = as.Date("2020-01-01")
  )

  expect_true(is.matrix(out))
  expect_equal(dim(out), fx$dims[c(1, 2)])
  expect_equal(out[fx$idx], rep(1, length(fx$idx)))

  # KNOWN BUG / INCONSISTENCY: in the "no usable range" early-return branch,
  # NA bathymetry cells are left at 0 (never set to NA), unlike the normal
  # path below, where idx_na cells are explicitly set to NA. This test
  # documents the current (inconsistent) behaviour so a fix is caught here.
  expect_equal(out[fx$idx_na], rep(0, length(fx$idx_na)))
})

test_that("calc_daily_likdepth treats NA depth_max the same as an empty log", {
  fx <- make_bathy_fixture()
  sdi_dsd <- .convert_spr_to_vec(fx$sdi) + 0.57

  na_log <- make_summary_table_2(as.Date("2020-01-01"), 5, NA_real_)

  out <- calc_daily_likdepth(
    summary_log_i = na_log,
    zbi = fx$zbi_vec,
    sdi_dsd = sdi_dsd,
    idx = fx$idx,
    idx_na = fx$idx_na,
    dims = fx$dims,
    time = as.Date("2020-01-01")
  )

  expect_equal(out[fx$idx], rep(1, length(fx$idx)))
})

test_that("calc_daily_likdepth treats Inf depth_max the same as an empty log", {
  fx <- make_bathy_fixture()
  sdi_dsd <- .convert_spr_to_vec(fx$sdi) + 0.57

  inf_log <- make_summary_table_2(as.Date("2020-01-01"), 5, Inf)

  out <- calc_daily_likdepth(
    summary_log_i = inf_log,
    zbi = fx$zbi_vec,
    sdi_dsd = sdi_dsd,
    idx = fx$idx,
    idx_na = fx$idx_na,
    dims = fx$dims,
    time = as.Date("2020-01-01")
  )

  expect_equal(out[fx$idx], rep(1, length(fx$idx)))
})

test_that("calc_daily_likdepth computes normalized normal-integration likelihood, sets NA cells to NA, and normalizes to max 1", {
  fx <- make_bathy_fixture()
  daily_sd_depth <- 0.57
  sdi_dsd <- .convert_spr_to_vec(fx$sdi) + daily_sd_depth

  log_i <- make_summary_table_2(as.Date("2020-01-01"), 15, 25)

  out <- calc_daily_likdepth(
    summary_log_i = log_i,
    zbi = fx$zbi_vec,
    sdi_dsd = sdi_dsd,
    idx = fx$idx,
    idx_na = fx$idx_na,
    dims = fx$dims,
    time = as.Date("2020-01-01")
  )

  # reproduce the expected bounds/integration independently from stats::pnorm
  low <- 15 - (1 + 0.01 * abs(15))
  high <- 25 + (1 + 0.01 * abs(25))
  expected_raw <- pmax(
    0,
    stats::pnorm(high, mean = fx$zbi_vec[fx$idx], sd = sdi_dsd[fx$idx]) -
      stats::pnorm(low, mean = fx$zbi_vec[fx$idx], sd = sdi_dsd[fx$idx])
  )
  expected_normalized <- expected_raw / max(expected_raw)

  expect_true(is.matrix(out))
  expect_equal(dim(out), fx$dims[1:2])
  expect_equal(out[fx$idx], expected_normalized, tolerance = 1e-8)

  # bathymetry NA cell must be NA in the output (normal-path behaviour)
  expect_true(all(is.na(out[fx$idx_na])))

  # peak of the normalized likelihood should be exactly 1
  expect_equal(max(out, na.rm = TRUE), 1)
})

test_that("calc_daily_likdepth does not divide by zero / produce NaN when all integrated values are exactly zero", {
  fx <- make_bathy_fixture()
  # tiny sd + a depth window nowhere near any bathymetry value => pnorm
  # differences underflow to exactly 0 for every idx cell
  sdi_dsd <- rep(0.001, length(fx$zbi_vec))

  far_log <- make_summary_table_2(as.Date("2020-01-01"), 100000, 100001)

  out <- calc_daily_likdepth(
    summary_log_i = far_log,
    zbi = fx$zbi_vec,
    sdi_dsd = sdi_dsd,
    idx = fx$idx,
    idx_na = fx$idx_na,
    dims = fx$dims,
    time = as.Date("2020-01-01")
  )

  expect_false(any(is.nan(out)))
  expect_equal(out[fx$idx], rep(0, length(fx$idx)))
})

test_that("calc_daily_likdepth wraps likint3_matrix failures in an informative 
cli_abort naming the date", {
  fx <- make_bathy_fixture()
  sdi_dsd <- .convert_spr_to_vec(fx$sdi) + 0.57

  # idx out of range for zbi/sdi_dsd length forces an internal subsetting
  # error inside likint3_matrix()
  bad_idx <- fx$idx + 1000L

  log_i <- make_summary_table_2(as.Date("2020-06-15"), 15, 25)

  expect_error(
    calc_daily_likdepth(
      summary_log_i = log_i,
      zbi = fx$zbi_vec,
      sdi_dsd = sdi_dsd,
      idx = bad_idx,
      idx_na = fx$idx_na,
      dims = fx$dims,
      time = as.Date("2020-06-15")
    ),
    regexp = "2020-06-15",
    fixed = TRUE
  )
})


# =============================================================================
# calc_likdepth_transmit() -- integration tests
# =============================================================================

test_that("calc_likdepth_transmit errors immediately because of the rlang::arg_match misuse (KNOWN BUG)", {
  # `preped_log` and `preped_bathy` have no default values in the function
  # signature, so rlang::arg_match() cannot derive `values` and fails before
  # error_list_prep()/error_numeric() -- i.e. any real input -- ever runs.
  # This test exists to catch a regression if that bug is ever silently
  # "fixed" by something even more wrong, and to make the current broken
  # state visible.
  expect_error(
    calc_likdepth_transmit(preped_log = list(), preped_bathy = list())
  )
})

test_that("calc_likdepth_transmit produces a correctly-shaped, correctly-normalized likelihood array end to end", {
  fx <- make_bathy_fixture()

  preped_bathy <- list(
    bathy = fx$bathy,
    sdi = fx$sdi,
    idx = fx$idx,
    idx_na = fx$idx_na
  )

  dates <- as.Date(c("2020-01-01", "2020-01-02", "2020-01-03", "2020-01-04"))
  preped_log <- list(
    summary_table = make_summary_table_2(
      dates = dates,
      depth_min = c(15, NA, 20, 25),
      depth_max = c(25, NA, 30, 35)
    )
  )

  event_log_test <- data.frame(
    fish_id = rep('1234', 2),
    date = c(
      as.POSIXct("2020-01-01 00:00:00"),
      as.POSIXct("2020-01-04 00:00:00")
    ),
    lon = c(-59.10835, -59.04600),
    lat = c(44.33739, 44.23701),
    event = c("tag", "pop")
  )
  # bypass the broken validation lines and stub out the msg_* helpers, whose
  # implementations live elsewhere in the package and aren't under test here
  testthat::local_mocked_bindings(
    arg_match = function(arg, ...) invisible(NULL),
    .package = "rlang"
  )
  testthat::local_mocked_bindings(
    error_list_prep = function(...) invisible(NULL),
    error_numeric = function(...) invisible(NULL),
    msg_start = function(...) Sys.time(),
    msg_dates = function(...) invisible(NULL),
    msg_end = function(...) invisible(NULL)
  )

  out <- calc_likdepth_transmit(
    preped_log,
    preped_bathy,
    event_log = event_log_test,
    daily_sd_depth = 0.57
  )

  # dims: aperm(c(2,1,3)) swaps rows/cols relative to the raw
  # array(dims = c(nrow, ncol, n_days))
  expect_equal(dim(out), c(fx$dims[2], fx$dims[1], length(dates)))

  # day 1 has a usable range -> peak likelihood should be 1 somewhere
  expect_equal(max(out[,, 1], na.rm = TRUE), 1)

  # day 2 has NA depth_max -> uninformative likelihood: every non-NA
  # bathymetry cell should equal 1
  day2 <- out[,, 2]
  expect_true(all(day2 %in% c(0, 1)))
  expect_equal(sum(day2 == 1), length(fx$idx))
  expect_equal(sum(day2 == 0), length(fx$idx_na))
})

test_that("calc_likdepth_transmit orders days by date regardless of input row order", {
  fx <- make_bathy_fixture()
  preped_bathy <- list(
    bathy = fx$bathy,
    sdi = fx$sdi,
    idx = fx$idx,
    idx_na = fx$idx_na
  )

  # deliberately out of chronological order
  dates <- as.Date(c("2020-01-03", "2020-01-01", "2020-01-02", "2020-01-04"))
  preped_log <- list(
    summary_table = make_summary_table_2(
      dates = dates,
      depth_min = c(NA, 15, NA, 20),
      depth_max = c(NA, 25, NA, 30)
    )
  )

  event_log_test <- data.frame(
    fish_id = rep('1234', 2),
    date = c(
      as.POSIXct("2020-01-01 00:00:00"),
      as.POSIXct("2020-01-04 00:00:00")
    ),
    lon = c(-59.10835, -59.04600),
    lat = c(44.33739, 44.23701),
    event = c("tag", "pop")
  )

  testthat::local_mocked_bindings(
    arg_match = function(arg, ...) invisible(NULL),
    .package = "rlang"
  )
  testthat::local_mocked_bindings(
    error_list_prep = function(...) invisible(NULL),
    error_numeric = function(...) invisible(NULL),
    msg_start = function(...) Sys.time(),
    msg_dates = function(...) invisible(NULL),
    msg_end = function(...) invisible(NULL)
  )

  out <- calc_likdepth_transmit(
    preped_log,
    preped_bathy,
    event_log = event_log_test,
    daily_sd_depth = 0.57
  )

  # after arrange(date), 2020-01-01 (the only informative day) is slice 1
  slice1 <- out[,, 1]
  expect_equal(max(slice1, na.rm = TRUE), 1)

  # the two NA-range days (now slices 2 and 3) go through the early-return
  # branch, where idx_na cells come through as 0 rather than NA (documented
  # bug) -- check both values explicitly rather than filtering on is.na()
  slice2 <- out[,, 2]
  expect_true(all(slice2 %in% c(0, 1)))
  expect_equal(sum(slice2 == 1), length(fx$idx))
  expect_equal(sum(slice2 == 0), length(fx$idx_na))

  slice3 <- out[,, 3]
  expect_true(all(slice3 %in% c(0, 1)))
  expect_equal(sum(slice3 == 1), length(fx$idx))
  expect_equal(sum(slice3 == 0), length(fx$idx_na))
})

test_that("calc_daily_likdepth attaches date/depth context to underlying likelihood errors", {
  fx <- make_bathy_fixture()
  sdi_dsd <- .convert_spr_to_vec(fx$sdi) + 0.57
  bad_idx <- fx$idx + 1000L
  log_i <- make_summary_table_2(as.Date("2020-06-15"), 15, 25)

  err <- rlang::catch_cnd(
    calc_daily_likdepth(
      summary_log_i = log_i,
      zbi = fx$zbi_vec,
      sdi_dsd = sdi_dsd,
      idx = bad_idx,
      idx_na = fx$idx_na,
      dims = fx$dims,
      time = as.Date("2020-06-15")
    )
  )

  expect_match(conditionMessage(err), "2020-06-15", fixed = TRUE)
  expect_match(conditionMessage(err), "15", fixed = TRUE)
  expect_match(conditionMessage(err), "25", fixed = TRUE)
  # confirm it's actually chained to the underlying out-of-bounds error,
  # not just coincidentally mentioning the date
  expect_match(conditionMessage(err$parent), "out of bounds", fixed = TRUE)
})
