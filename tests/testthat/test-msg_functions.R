test_that("Check if msg_dates return what it is suppose to", {
  df <- data.frame(
    date = seq.Date(
      as.Date("2026-08-01"),
      as.Date("2026-08-31"),
      length.out = 31
    )
  )
  n_days <- nrow(df)

  expect_message(
    msg_dates(df, n_days),
    regexp = "Generating depth profile likelihood",
    fixed = TRUE
  )
})
