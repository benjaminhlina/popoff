test_that(".has_valid_field returns FALSE when field is missing", {
  x <- list(a = 1)
  expect_false(.has_valid_field(x, "b", "numeric"))
})

test_that(".has_valid_field returns TRUE when field matches class check", {
  x <- list(bathy = terra::rast())
  expect_true(.has_valid_field(x, "bathy", "SpatRaster"))
})

test_that(".has_valid_field returns FALSE when field does not match class check", {
  x <- list(bathy = "not a raster")
  expect_false(.has_valid_field(x, "bathy", "SpatRaster"))
})

test_that(".has_valid_field works with a predicate function", {
  x <- list(idx = c(1, 2, 3))
  expect_true(.has_valid_field(x, "idx", is.numeric))

  x2 <- list(idx = "a")
  expect_false(.has_valid_field(x2, "idx", is.numeric))
})

test_that(".has_valid_field errors when check is neither a class string nor a function", {
  x <- list(a = 1)
  expect_error(
    .has_valid_field(x, "a", 42),
    "class string or a predicate function"
  )
})

test_that(".has_valid_field coerces non-scalar-TRUE predicate results to FALSE", {
  # isTRUE() guards against predicates that return a non-length-1 logical
  x <- list(a = c(1, 2))
  bad_predicate <- function(val) c(TRUE, TRUE)
  expect_false(.has_valid_field(x, "a", bad_predicate))
})
