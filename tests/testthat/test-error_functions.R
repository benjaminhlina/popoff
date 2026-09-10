test_that("check if error_chracter works properly", {
  expect_null(error_chr("c"))

  expect_error(
    error_chr(6),
    "must be a",
    fixed = TRUE
  )
})


test_that("check if error_date works properly", {
  date <- as.Date("2026-09-10")
  expect_null(
    error_date(date, col = "date", expected_cls = "Date")
  )
  pos <- as.POSIXct("2026-09-10 14:16:00")
  expect_null(
    error_date(pos, col = "date_time", expected_cls = "POSIXct")
  )

  expect_error(
    error_date(date, col = "date", expected_cls = "POSIXct"),
    "must inherit from",
    fixed = TRUE
  )
  expect_error(
    error_date(pos, col = "date", expected_cls = "date"),
    "must inherit from",
    fixed = TRUE
  )
  expect_error(
    error_date(date, expected_cls = "date"),
    'argument "col" is missing, with no default',
    fixed = TRUE
  )
})


test_that("check if error_df works properly", {
  df <- data.frame(
    time = seq(1, 10, 1)
  )

  expect_null(error_df(df))

  list_df <- list(df = df)

  expect_error(
    error_df(
      list_df
    ),
    "with a data.frame",
    fixed = TRUE
  )
})

test_that("check if error_list works properly", {
  vec <- seq(1, 10, 1)
  list_df <- list(df = df)

  expect_null(error_list(list_df, "test"))

  expect_error(
    error_list(
      vec,
      "vec"
    ),
    "must be a",
    fixed = TRUE
  )
})

# test_that("check if error_raster works properly", {
#   bath <- raster::raster(bathy_matrix)

#   expect_null(
#     error_raster(bath)
#   )
#   expect_error(
#     error_raster(bathy_matrix),
#     "You supplied ",
#     fixed = TRUE
#   )
# })
