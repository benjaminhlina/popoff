test_that(".compute_depth_range returns NA,NA when depth is entirely NA", {
  expect_equal(
    .compute_depth_range(depth = c(NA, NA, NA), depth_run = c(NA, NA, NA)),
    c(NA_real_, NA_real_)
  )
})

test_that(".compute_depth_range returns range of depth_run when depth has values", {
  expect_equal(
    .compute_depth_range(depth = c(10, NA, 30), depth_run = c(10, 20, 30)),
    c(10, 30)
  )
})

test_that(".compute_temp_subset returns numeric(0) when depth is entirely NA", {
  expect_equal(
    .compute_temp_subset(
      depth = c(NA, NA),
      temperature = c(10, 12),
      depth_min = NA,
      depth_max = NA
    ),
    numeric(0)
  )
})

test_that(".compute_temp_subset returns numeric(0) when no valid temperatures exist", {
  expect_equal(
    .compute_temp_subset(
      depth = c(10, 20),
      temperature = c(NA, NA),
      depth_min = 10,
      depth_max = 20
    ),
    numeric(0)
  )
})

test_that(".compute_temp_subset filters temperatures within the depth range", {
  # depth range covers rows 1 and 3 only; row 2's depth (99) falls outside
  result <- .compute_temp_subset(
    depth = c(10, 99, 20),
    temperature = c(5, 100, 8),
    depth_min = 10,
    depth_max = 20
  )
  expect_equal(result, c(5, 8))
})

test_that(".compute_temp_subset handles reversed depth_min/depth_max", {
  # min/max passed in swapped order should still bound correctly
  result <- .compute_temp_subset(
    depth = c(10, 15, 20),
    temperature = c(5, 6, 7),
    depth_min = 20,
    depth_max = 10
  )
  expect_equal(result, c(5, 6, 7))
})

test_that(".compute_temp_bound uses temp_subset when non-empty", {
  expect_equal(
    .compute_temp_bound(
      temp_subset = c(5, 8, 12),
      valid_temps = c(1, 5, 8, 12, 99),
      type = "min"
    ),
    5
  )
  expect_equal(
    .compute_temp_bound(
      temp_subset = c(5, 8, 12),
      valid_temps = c(1, 5, 8, 12, 99),
      type = "max"
    ),
    12
  )
})

test_that(".compute_temp_bound falls back to .safe_range_val when temp_subset is empty", {
  # temp_subset empty but valid_temps has data outside the depth range -
  # bound should come from the full valid_temps pool, not temp_subset
  expect_equal(
    .compute_temp_bound(
      temp_subset = numeric(0),
      valid_temps = c(1, 5, 8),
      type = "min"
    ),
    1
  )
  expect_equal(
    .compute_temp_bound(
      temp_subset = numeric(0),
      valid_temps = c(1, 5, 8),
      type = "max"
    ),
    8
  )
})

test_that(".compute_temp_bound returns NA when both temp_subset and valid_temps are empty", {
  expect_equal(
    .compute_temp_bound(
      temp_subset = numeric(0),
      valid_temps = numeric(0),
      type = "min"
    ),
    NA_real_
  )
})
