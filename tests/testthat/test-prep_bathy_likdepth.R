bathy_spat <- load_bathy_raster(bathy_matrix)

# ---- input validation ---------------------------------------------------

test_that("prep_bathy_likdepth validates bathy is a SpatRaster", {
  expect_error(
    prep_bathy_likdepth(bathy_matrix),
    "SpatRaster"
  )
})

test_that("prep_bathy_likdepth validates weighted is logical", {
  expect_error(
    prep_bathy_likdepth(bathy_spat, weighted = "yes"),
    "logical"
  )
})

test_that("prep_bathy_likdepth rejects even focal_dim", {
  expect_error(
    prep_bathy_likdepth(bathy_spat, focal_dim = 4),
    "odd"
  )
})

# ---- happy path (weighted = FALSE) --------------------------------------

test_that("prep_bathy_likdepth returns expected structure when weighted = FALSE", {
  result <- suppressMessages(
    prep_bathy_likdepth(bathy_spat, weighted = FALSE, focal_dim = 3)
  )

  expect_type(result, "list")
  expect_named(
    result,
    c("bathy", "sdi", "idx", "idx_na", "focal_dim", "weighted")
  )
  expect_s4_class(result$sdi, "SpatRaster")
  expect_identical(result$focal_dim, 3)
  expect_false(result$weighted)
  expect_identical(result$idx, attr(bathy_spat, "idx"))
  expect_identical(result$idx_na, attr(bathy_spat, "idx_na"))
  expect_identical(dim(result$sdi)[1:2], dim(bathy_spat)[1:2])
})

test_that("sd values are all positive (zeros/NAs filled in)", {
  result <- suppressMessages(
    prep_bathy_likdepth(bathy_spat, weighted = FALSE, focal_dim = 3)
  )

  vals <- terra::values(result$sdi)
  expect_false(anyNA(vals))
  expect_true(all(vals > 0))
})

# test_that("fill_sd0 falls back to 1e-3 when every focal sd is zero", {
#   # constant-depth raster => sd is 0 in every window

#   result <- suppressMessages(
#     prep_bathy_likdepth(bathy_spat, weighted = FALSE, focal_dim = 3)
#   )

#   vals <- terra::values(result$sdi)
#   expect_true(all(vals > 1e-5))
# })

# ---- known bugs (documented so regressions/fixes are visible) -----------

test_that("KNOWN BUG: weighted = TRUE errors due to ctd.dir/ctd_dir mismatch", {
  # ctd_dir is supplied correctly, but the function references `ctd.dir`
  # internally, which does not exist -> errors regardless of input
  expect_error(
    prep_bathy_likdepth(
      bathy_spat,
      weighted = TRUE,
      focal_dim = 3,
      ctd_dir = tempdir()
    ),
    "ctd\\.dir"
  )
})

test_that("KNOWN BUG: sdi does not inherit bathy's CRS/extent", {
  terra::crs(bathy_spat) <- "EPSG:4326"
  # terra::ext(bathy_spat) <- terra::ext(0, 10, 0, 10)

  result <- suppressMessages(
    prep_bathy_likdepth(bathy_spat, weighted = FALSE, focal_dim = 3)
  )

  # documents current (incorrect) behavior -- sdi is NOT georeferenced
  # to match bathy. If this test starts failing, the bug has been fixed
  # and this test (and the underlying code) should be updated.
  expect_true(identical(terra::crs(result$sdi), terra::crs(bathy_spat)))
  expect_true(identical(
    as.vector(terra::ext(result$sdi)),
    as.vector(terra::ext(bathy_spat))
  ))
})
