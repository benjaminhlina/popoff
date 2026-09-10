test_that("check if error_chracter works proplery", {
  expect_null(error_chr("c"))

  expect_error(
    error_chr(6),
    "must be a",
    fixed = TRUE
  )
})


test_that("check if error_date works proplery", {
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
