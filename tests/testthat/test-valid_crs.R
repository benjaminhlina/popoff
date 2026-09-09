test_that(".valid_crs returns FALSE (not NULL) for invalid crs", {
  result <- .valid_crs("not_a_real_crs")
  expect_type(result, "logical")
  expect_false(result)
})

test_that(".valid_crs returns TRUE for a valid crs", {
  expect_true(.valid_crs("EPSG:4326"))
})
