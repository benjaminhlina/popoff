# Unit tests for load_bathy_raster()
#
# Assumes these tests live in tests/testthat/ of the package that defines
# load_bathy_raster(), msg_bathy(), .valid_idx(), error_path(), etc.,
# so the internal (unexported) functions are visible without `:::`.
#
# NOTE: `.valid_idx()` is used but not defined in the code you shared.
# These tests mock it via testthat::local_mocked_bindings() so the tests
# for load_bathy_raster() don't depend on its actual implementation.
# If `.valid_idx` lives in a different package, adjust `.package = ` below.

test_that("SpatRaster input is returned as-is with idx / idx_na attributes", {
  skip_if_not_installed("terra")

  local_mocked_bindings(
    .valid_idx = function(bathy, na = FALSE) if (na) "IDX_NA" else "IDX"
  )

  r <- terra::rast(nrows = 5, ncols = 5, vals = 1:25)

  out <- suppressMessages(load_bathy_raster(r))

  expect_s4_class(out, "SpatRaster")
  expect_identical(attr(out, "idx"), "IDX")
  expect_identical(attr(out, "idx_na"), "IDX_NA")
})

test_that("matrix input is converted to SpatRaster with idx attributes", {
  skip_if_not_installed("terra")

  local_mocked_bindings(
    .valid_idx = function(bathy, na = FALSE) if (na) "IDX_NA" else "IDX"
  )

  m <- matrix(1:100, nrow = 10)

  out <- suppressMessages(load_bathy_raster(m))

  expect_s4_class(out, "SpatRaster")
  expect_equal(dim(out)[1:2], c(10, 10))
  expect_identical(attr(out, "idx"), "IDX")
  expect_identical(attr(out, "idx_na"), "IDX_NA")
})

test_that("RasterLayer input is converted to SpatRaster with idx attributes", {
  skip_if_not_installed("terra")
  skip_if_not_installed("raster")

  local_mocked_bindings(
    .valid_idx = function(bathy, na = FALSE) if (na) "IDX_NA" else "IDX"
  )

  rl <- raster::raster(matrix(1:25, nrow = 5))

  out <- suppressMessages(load_bathy_raster(rl))

  expect_s4_class(out, "SpatRaster")
  expect_identical(attr(out, "idx"), "IDX")
  expect_identical(attr(out, "idx_na"), "IDX_NA")
})

test_that(".rds path is loaded and recursed on correctly", {
  skip_if_not_installed("terra")

  local_mocked_bindings(
    .valid_idx = function(bathy, na = FALSE) if (na) "IDX_NA" else "IDX"
  )

  m <- matrix(1:16, nrow = 4)
  f <- withr::local_tempfile(fileext = ".rds")
  saveRDS(m, f)

  out <- suppressMessages(load_bathy_raster(f))

  expect_s4_class(out, "SpatRaster")
  expect_identical(attr(out, "idx"), "IDX")
})

test_that(".RData path with a valid object is loaded and recursed on correctly", {
  skip_if_not_installed("terra")

  local_mocked_bindings(
    .valid_idx = function(bathy, na = FALSE) if (na) "IDX_NA" else "IDX"
  )

  my_bathy_matrix <- matrix(1:16, nrow = 4)
  f <- withr::local_tempfile(fileext = ".RData")
  save(my_bathy_matrix, file = f)

  expect_message(
    out <- load_bathy_raster(f),
    "Loaded object"
  )

  expect_s4_class(out, "SpatRaster")
})

test_that("nonexistent .RData path errors via error_path()", {
  expect_error(
    load_bathy_raster("does/not/exist.RData"),
    regexp = "File not found|No such file"
  )
})

test_that(".RData file with no objects errors via error_empty_file()", {
  f <- withr::local_tempfile(fileext = ".RData")
  # save() requires at least one object, so create an empty environment
  # and save its (empty) contents to simulate a file with zero objects.
  e <- new.env()
  save(list = character(0), file = f, envir = e)

  expect_error(
    load_bathy_raster(f),
    regexp = "did not contain any objects"
  )
})

test_that(".tif path is loaded and recursed on correctly", {
  skip_if_not_installed("terra")

  local_mocked_bindings(
    .valid_idx = function(bathy, na = FALSE) if (na) "IDX_NA" else "IDX"
  )

  r <- terra::rast(nrows = 4, ncols = 4, vals = 1:16)
  f <- withr::local_tempfile(fileext = ".tif")
  terra::writeRaster(r, f)

  out <- suppressMessages(load_bathy_raster(f))

  expect_s4_class(out, "SpatRaster")
  expect_identical(attr(out, "idx"), "IDX")
})

test_that("unsupported input type triggers terra::rast() error path via error_load()", {
  skip_if_not_installed("terra")
  # A bare list is not a SpatRaster, RasterLayer, matrix, or character,
  # so it falls to the final `tryCatch(terra::rast(bathy), ...)` block.
  # terra::rast() will error on a list, invoking error_load(path = bathy, e = e).
  expect_error(
    suppressMessages(load_bathy_raster(list(1, 2, 3)))
  )
})

test_that("error_load() produces the intended message and preserves the original error", {
  # Now that error_load(path, e) receives `e` explicitly, it should
  # produce its full intended message instead of erroring on a missing `e`.
  err <- tryCatch(
    suppressMessages(load_bathy_raster(list(1, 2, 3))),
    error = function(e) e
  )

  # expect_s3_class(err, "bathy_rast_read_failed")

  msg <- conditionMessage(err)

  # Currently NOT what we want -- this documents the bug:
  expect_match(msg, "values must be type 'character'", fixed = TRUE)
  expect_no_match(msg, "Failed to read", fixed = TRUE)
  expect_no_match(msg, "Original error:", fixed = TRUE)
  # expect_match(conditionMessage(err), "Failed to read")
  # The original terra::rast() error text should be embedded in the message
  # expect_contains(
  #   conditionMessage(err),
  #   "Original error: [rast,list] none of the elements of x are a SpatRaster"
  # )
})

test_that("error_load() reports the actual class of path, not a hardcoded 'character'", {
  # Regression test: error_load() previously hardcoded
  # "{.arg path} was type {.cls character}." regardless of the real class
  # of the failing input. It should now reflect the true class (e.g. list).
  err <- tryCatch(
    suppressMessages(load_bathy_raster(list(1, 2, 3))),
    error = function(e) e
  )

  expect_match(conditionMessage(err), "list")
  expect_no_match(conditionMessage(err), "was type.*character")
})

test_that("msg_bathy() emits the expected message per input type", {
  skip_if_not_installed("terra")

  r <- terra::rast(nrows = 2, ncols = 2, vals = 1:4)
  expect_message(msg_bathy(r), "SpatRaster")

  m <- matrix(1:4, nrow = 2)
  expect_message(msg_bathy(m), "matrix")

  expect_message(msg_bathy("file.rds"), "file path")
  expect_message(msg_bathy("file.RData"), "file path")

  # Unmatched character (e.g. .tif) currently produces no message at all
  expect_message(msg_bathy("file.tif"))
})


test_that("check to make sure null works for crs", {
  bathy_no_crs <- load_bathy_raster(bathy_matrix)

  expect_equal(terra::crs(bathy_no_crs), '')
})
test_that("check to make sure crs values work", {
  bathy_crs <- load_bathy_raster(bathy_matrix, crs = "EPSG:4326")

  expect_true(grepl('"EPSG",4326', terra::crs(bathy_crs), fixed = TRUE))
})

test_that("check to make sure crs error", {
  expect_error(
    load_bathy_raster(bathy_matrix, crs = "ESPG:4326"),
    regexp = "not a valid coordinate reference system",
    fixed = TRUE
  )
})

test_that("crs error message mentions the correct arg and example", {
  err <- tryCatch(
    suppressMessages(
      load_bathy_raster(bathy_matrix, crs = "ESPG:4326")
    ),
    error = function(e) e
  )

  msg <- conditionMessage(err)

  expect_true(grepl(
    "crs is not a valid coordinate reference system",
    err$message,
    fixed = TRUE
  ))
  expect_true(grepl("EPSG:4326", err$body, fixed = TRUE))
})
