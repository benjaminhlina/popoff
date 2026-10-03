# tests/testthat/test-download_noaa_bathy.R

# helpers ----------------------------------------------------------------------

# marmap "bathy" objects: rows = lon, cols = lat, depth negative below sea level
fake_bathy <- function(
  lon = seq(-80, -79, by = 0.25),
  lat = seq(44, 45, by = 0.25)
) {
  z <- outer(seq_along(lon), seq_along(lat), \(i, j) -(100 * i + j))
  dimnames(z) <- list(as.character(lon), as.character(lat))
  class(z) <- "bathy"
  z
}

local_mock_noaa <- function(bathy = fake_bathy(), env = parent.frame()) {
  calls <- new.env()
  calls$n <- 0L
  calls$args <- NULL

  testthat::local_mocked_bindings(
    getNOAA.bathy = function(lon1, lon2, lat1, lat2, ...) {
      calls$n <- calls$n + 1L
      calls$args <- c(
        list(lon1 = lon1, lon2 = lon2, lat1 = lat1, lat2 = lat2),
        list(...)
      )
      bathy
    },
    .package = "marmap",
    .env = env
  )
  calls
}

xlim <- c(-80, -79)
ylim <- c(44, 45)

# input validation -------------------------------------------------------------

test_that("xlim must be a numeric vector of length 2", {
  calls <- local_mock_noaa()

  expect_error(download_noaa_bathy(c(-80, -79, -78), ylim), "xlim")
  expect_error(download_noaa_bathy(-80, ylim), "xlim")
  expect_error(download_noaa_bathy(c("a", "b"), ylim), "xlim")
  expect_error(download_noaa_bathy(c(-80, NA), ylim), "xlim")
  expect_error(download_noaa_bathy(NULL, ylim), "xlim")

  expect_equal(calls$n, 0L) # validation happens before any download
})

test_that("ylim must be a numeric vector of length 2", {
  calls <- local_mock_noaa()

  expect_error(download_noaa_bathy(xlim, c(44, 45, 46)), "ylim")
  expect_error(download_noaa_bathy(xlim, 44), "ylim")
  expect_error(download_noaa_bathy(xlim, c("a", "b")), "ylim")
  expect_error(download_noaa_bathy(xlim, c(NA, 45)), "ylim")

  expect_equal(calls$n, 0L)
})

test_that("land_value must be numeric", {
  calls <- local_mock_noaa()

  expect_error(download_noaa_bathy(xlim, ylim, land_value = "high"))
  expect_equal(calls$n, 0L)
})

# arguments passed to marmap ----------------------------------------------------

test_that("limits are passed to getNOAA.bathy as min/max", {
  calls <- local_mock_noaa()

  download_noaa_bathy(xlim, ylim)

  expect_equal(calls$n, 1L)
  expect_equal(calls$args$lon1, -80)
  expect_equal(calls$args$lon2, -79)
  expect_equal(calls$args$lat1, 44)
  expect_equal(calls$args$lat2, 45)
})

test_that("limits are order-independent", {
  calls <- local_mock_noaa()

  download_noaa_bathy(c(-79, -80), c(45, 44))

  expect_equal(calls$args$lon1, -80)
  expect_equal(calls$args$lon2, -79)
  expect_equal(calls$args$lat1, 44)
  expect_equal(calls$args$lat2, 45)
})

test_that("... is forwarded to getNOAA.bathy", {
  calls <- local_mock_noaa()

  download_noaa_bathy(xlim, ylim, resolution = 4, keep = TRUE)

  expect_equal(calls$args$resolution, 4)
  expect_true(calls$args$keep)
})

# output structure -------------------------------------------------------------

test_that("returns a plain numeric matrix with correct dimensions", {
  local_mock_noaa()

  res <- download_noaa_bathy(xlim, ylim)

  expect_true(is.matrix(res))
  expect_type(res, "double")
  expect_false(inherits(res, "bathy"))
  expect_equal(dim(res), c(5, 5))
})

test_that("longitude is ascending and latitude is descending", {
  local_mock_noaa()

  res <- download_noaa_bathy(xlim, ylim)

  lon <- as.numeric(colnames(res))
  lat <- as.numeric(rownames(res))

  expect_equal(lon, sort(lon))
  expect_equal(lat, sort(lat, decreasing = TRUE))
  expect_equal(range(lon), c(-80, -79))
  expect_equal(range(lat), c(44, 45))
})

test_that("values stay attached to the correct coordinates after reordering", {
  local_mock_noaa()

  res <- download_noaa_bathy(xlim, ylim)

  # fake depth = -(100 * lon_index + lat_index) in the original (ascending) grid
  # lon index 2 = -79.75, lat index 3 = 44.5
  expect_equal(res["44.5", "-79.75"], -203)
  # lon index 5 = -79, lat index 1 = 44
  expect_equal(res["44", "-79"], -501)
  # lon index 1 = -80, lat index 5 = 45
  expect_equal(res["45", "-80"], -105)
})

# land masking -----------------------------------------------------------------

test_that("cells above land_value are NA; default land_value is 0", {
  b <- fake_bathy()
  b[1, 1] <- 50 # lon -80, lat 44 -> land
  b[2, 2] <- 0 # at the threshold -> kept
  b[3, 3] <- -5 # shallow water -> kept
  local_mock_noaa(b)

  res <- download_noaa_bathy(xlim, ylim)

  expect_true(is.na(res["44", "-80"]))
  expect_equal(res["44.25", "-79.75"], 0)
  expect_equal(res["44.5", "-79.5"], -5)
  expect_equal(sum(is.na(res)), 1)
})

test_that("custom land_value changes the mask", {
  b <- fake_bathy()
  b[1, 1] <- 50
  b[2, 2] <- 0
  b[3, 3] <- -5
  local_mock_noaa(b)

  res <- download_noaa_bathy(xlim, ylim, land_value = -10)

  expect_true(is.na(res["44", "-80"]))
  expect_true(is.na(res["44.25", "-79.75"]))
  expect_true(is.na(res["44.5", "-79.5"]))
  expect_equal(sum(is.na(res)), 3)
})

test_that("land_value = NULL falls back to 0", {
  b <- fake_bathy()
  b[1, 1] <- 1
  local_mock_noaa(b)

  expect_equal(
    download_noaa_bathy(xlim, ylim, land_value = NULL),
    download_noaa_bathy(xlim, ylim, land_value = 0)
  )
})

# extent attributes ------------------------------------------------------------

test_that("xlim/ylim attributes are cell edges (centres +/- half a cell)", {
  local_mock_noaa() # 0.25 degree cells

  res <- download_noaa_bathy(xlim, ylim)

  expect_equal(attr(res, "xlim"), c(-80.50, -78.50))
  expect_equal(attr(res, "ylim"), c(43.50, 45.50))
})

test_that("attributes scale with cell size and handle unequal x/y resolution", {
  local_mock_noaa(fake_bathy(
    lon = seq(-80, -79, by = 0.5),
    lat = seq(44, 45, by = 0.1)
  ))

  res <- download_noaa_bathy(xlim, ylim)

  expect_equal(attr(res, "xlim"), c(-80.50, -78.50))
  expect_equal(attr(res, "ylim"), c(43.50, 45.50))
})

test_that("attributes round-trip into a correctly registered SpatRaster", {
  skip_if_not_installed("terra")
  local_mock_noaa()

  res <- download_noaa_bathy(xlim, ylim)
  r <- terra::rast(
    res,
    extent = terra::ext(c(attr(res, "xlim"), attr(res, "ylim")))
  )

  expect_equal(terra::res(r), c(0.40, 0.40))
  expect_equal(
    as.vector(terra::ext(r)),
    c(xmin = -80.50, xmax = -78.50, ymin = 43.50, ymax = 45.50),
    ignore_attr = TRUE
  )
})

# live download (optional) -----------------------------------------------------

test_that("live NOAA download returns a well-formed matrix", {
  skip_on_cran()
  skip_on_ci()
  skip_if_offline()
  skip_if_not_installed("marmap")

  res <- download_noaa_bathy(c(-80, -79.5), c(44.5, 45))

  expect_true(is.matrix(res))
  expect_gt(nrow(res), 1)
  expect_gt(ncol(res), 1)
  expect_length(attr(res, "xlim"), 2)
  expect_length(attr(res, "ylim"), 2)
  expect_lt(attr(res, "xlim")[1], attr(res, "xlim")[2])
})
