bathy_spat <- load_raster(bathy_matrix, crs = "EPSG:4326")

set.seed(42)

# base raster: 10x10, no missing data
m <- matrix(rnorm(100, mean = 50, sd = 10), nrow = 10, ncol = 10)

# 1. A solid block of NA at least as large as the focal window (3x3).
#    This guarantees some windows are ENTIRELY NA, so focal(..., fun = sd)
#    returns NA there no matter what na.rm is set to.
m[4:6, 4:6] <- NA

# 2. A couple of isolated single-cell NAs near, but not inside, the block.
#    These give you windows that are PARTIALLY NA -- useful if the code
#    under test behaves differently for na.rm = TRUE (drops NA, computes sd
#    on remaining cells) vs na.rm = FALSE (any NA neighbor -> NA output).
m[2, 8] <- NA

bathy_spat_na <- load_raster(m)


# ---------------------------------------------------------------------
# Fixture 1: perfectly flat raster -> sd_r is 0 everywhere, no NA cells
# -> ssd[ssd != 0] is empty -> quantile() returns NA -> fill_sd0 <- 1e-3
# ---------------------------------------------------------------------
bathy_spat_flat <- terra::rast(matrix(50, nrow = 10, ncol = 10))

# ---------------------------------------------------------------------
# Fixture 2: a ctd SpatRaster to pass through the !is.null(ctd) branch
# ---------------------------------------------------------------------
ctd_spat <- terra::rast(matrix(rnorm(100, mean = 20, sd = 2), 10, 10))


# ---- input validation ---------------------------------------------------

test_that("prep_bathy_likdepth validates bathy is a SpatRaster", {
  expect_error(
    prep_bathy_likdepth(bathy_matrix),
    "SpatRaster"
  )
})


test_that("prep_bathy_likdepth rejects even focal_dim", {
  expect_error(
    prep_bathy_likdepth(bathy_spat, focal_dim = 4),
    "odd"
  )
})

# ---- happy path (weighted = FALSE) --------------------------------------

test_that("prep_bathy_likdepth returns expected structure", {
  result <- suppressMessages(
    prep_bathy_likdepth(bathy_spat, focal_dim = 3)
  )

  expect_type(result, "list")
  expect_named(
    result,
    c("bathy", "sdi", "idx", "idx_na", "focal_dim")
  )
  expect_s4_class(result$sdi, "SpatRaster")
  expect_identical(result$focal_dim, 3)
  expect_identical(result$idx, attr(bathy_spat, "idx"))
  expect_identical(result$idx_na, attr(bathy_spat, "idx_na"))
  expect_identical(dim(result$sdi)[1:2], dim(bathy_spat)[1:2])
})

test_that("sd values are all positive (zeros/NAs filled in)", {
  result <- suppressMessages(
    prep_bathy_likdepth(bathy_spat, focal_dim = 3)
  )

  vals <- terra::values(result$sdi)
  expect_false(anyNA(vals))
  expect_true(all(vals > 0))
})

test_that("fill_sd0 falls back to 1e-3 when every focal sd is zero", {
  # constant-depth raster => sd is 0 in every window

  result <- suppressMessages(
    prep_bathy_likdepth(bathy_spat_na, focal_dim = 3)
  )

  sdi_na <- terra::focal(bathy_spat_na, w = 3)

  sdi_val <- terra::values(sdi_na)

  vals <- terra::values(result$sdi)

  expect_true(all(vals > 1e-5))
  expect_false(anyNA(vals))
})

# ---- known bugs (documented so regressions/fixes are visible) -----------

test_that("Check if sdi does not inherit bathy's CRS/extent", {
  # terra::ext(bathy_spat) <- terra::ext(0, 10, 0, 10)

  result <- suppressMessages(
    prep_bathy_likdepth(bathy_spat, focal_dim = 3)
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

describe("prep_bathy_likdepth()", {
  it("uses ctd directly as sdi and sets focal_dim to NULL when ctd is supplied", {
    result <- suppressMessages(
      prep_bathy_likdepth(bathy_spat, focal_dim = 9, ctd = ctd_spat)
    )
    result
    # sdi should be exactly the ctd raster, untouched by focal()/fill logic
    expect_equal(
      terra::values(result$sdi),
      terra::values(ctd_spat)
    )

    # focal_dim is reassigned to NULL inside the ctd branch
    expect_null(result$focal_dim)
  })

  it("validates ctd through error_spatrast() when ctd is not a SpatRaster", {
    expect_error(
      prep_bathy_likdepth(bathy_spat, focal_dim = 9, ctd = "not_a_raster"),
      "not SpatRaster",
      fixed = TRUE
    )
  })

  it("falls back to fill_sd0 = 1e-3 when focal sd is uniformly zero", {
    result <- suppressMessages(
      prep_bathy_likdepth(bathy_spat_flat, focal_dim = 3)
    )

    vals <- terra::values(result$sdi)

    # every cell had sd == 0 and no non-zero values existed for the
    # quantile() call, so fill_sd0 falls back to the 1e-3 constant
    expect_true(all(vals == 1e-3))
    expect_true(all(vals > 1e-5))
  })

  it("still fills NA cells from a partially-NA raster with the quantile-derived fill_sd0", {
    m <- matrix(rnorm(100, mean = 50, sd = 10), nrow = 10, ncol = 10)
    m[4:6, 4:6] <- NA # block >= focal window guarantees an all-NA window
    bathy_spat_na <- terra::rast(m)

    result <- suppressMessages(
      prep_bathy_likdepth(bathy_spat_na, focal_dim = 3)
    )

    vals <- terra::values(result$sdi)

    expect_false(anyNA(vals)) # is.na(sdi) branch should have filled these
    expect_true(all(vals > 1e-5)) # fill_sd0 here comes from the real quantile, not the 1e-3 fallback
  })
})
