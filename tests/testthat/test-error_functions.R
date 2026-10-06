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

  expect_null(error_list(list_df))

  expect_error(
    error_list(
      vec
    ),
    "must be a",
    fixed = TRUE
  )
})

test_that("check if error_matrix works properly", {
  expect_error(
    error_matrix(iniloc),
    "must be a matrix.",
    fixed = TRUE,
  )

  expect_error(
    error_matrix(
      matrix(c(2, 2, 2)),
      dim_check = c(3, 3)
    ),
    "contains index/indices out of range for a matrix",
    fixed = TRUE
  )

  error_matrix(
    matrix(c(2, 2, 2)),
    dim_check = c(1, 2)
  )
})

test_that("error_matrix: errors when a checked dimension is even", {
  even_matrix <- matrix(1:12, nrow = 4, ncol = 3)

  expect_error(
    error_matrix(even_matrix, dim_check = c(1, 2)),
    "must have odd dimensions",
    fixed = TRUE
  )
})

test_that("error_matrix: error message reports the actual offending dimensions", {
  even_matrix <- matrix(1:12, nrow = 4, ncol = 3)

  expect_error(
    error_matrix(even_matrix, dim_check = c(1, 2)),
    "4",
    fixed = TRUE
  )
})

test_that("error_matrix: passes silently when all checked dimensions are odd", {
  odd_matrix <- matrix(1:9, nrow = 3, ncol = 3)

  expect_no_error(
    error_matrix(odd_matrix, dim_check = c(1, 2))
  )
})

test_that("error_matrix: returns the matrix unchanged when validation passes", {
  odd_matrix <- matrix(1:9, nrow = 3, ncol = 3)

  result <- error_matrix(odd_matrix, dim_check = c(1, 2))

  expect_identical(result, odd_matrix)
})

test_that("error_matrix: errors if only one of several checked dimensions is even", {
  mixed_matrix <- matrix(1:15, nrow = 3, ncol = 5)
  # transpose-like case: check row count (3, odd) and a fabricated even one via dim_check reuse
  even_row_matrix <- matrix(1:20, nrow = 4, ncol = 5)

  # row dim (4) even, col dim (5) odd -> should still error due to `any()`
  expect_error(
    error_matrix(even_row_matrix, dim_check = c(1, 2)),
    "must have odd dimensions",
    fixed = TRUE
  )

  # both odd -> no error
  expect_no_error(
    error_matrix(mixed_matrix, dim_check = c(1, 2))
  )
})

test_that("error_matrix: dim_check on a single dimension only evaluates that dimension", {
  # 4 rows (even), 3 cols (odd); only check the odd one
  partial_matrix <- matrix(1:12, nrow = 4, ncol = 3)

  expect_no_error(
    error_matrix(partial_matrix, dim_check = 2)
  )

  expect_error(
    error_matrix(partial_matrix, dim_check = 1),
    "must have odd dimensions",
    fixed = TRUE
  )
})

test_that("error_matrix: no dimension check performed when dim_check is NULL", {
  even_matrix <- matrix(1:12, nrow = 4, ncol = 3)

  expect_no_error(
    error_matrix(even_matrix, dim_check = NULL)
  )
})


test_that("error_focal_dim: non-numeric matrix triggers error_numeric and does NOT fall through to error_matrix", {
  non_numeric_mat <- matrix(letters[1:4], nrow = 2)

  expect_error(
    error_focal_dim(non_numeric_mat),
    "needs to be a numeric",
    fixed = TRUE
  )
})

test_that("error_focal_dim: numeric matrix skips error_numeric and reaches error_matrix with correct dim_check", {
  numeric_mat <- matrix(1:4, nrow = 2)

  expect_error(
    error_focal_dim(numeric_mat),
    "The supplied matrix in",
    fixed = TRUE
  )

  expect_error(
    error_focal_dim(4),
    "has to be odd",
    fixed = TRUE
  )
})

test_that("error_matrix_likelihood errors properly when min is larger 
than max", {
  ma <- matrix(c(0, 1, 5, 2), nrow = 1, ncol = 4)

  expect_error(
    error_likelihood_matrix(ma, idx = 1),
    regexp = "Check for reversed or malformed",
    fixed = TRUE
  )
})


test_that("error_numeric_min works", {
  vec <- seq(1, 10, 1)

  expect_error(
    error_numeric_min(vec, col = "test", min_val = 2),
    regexp = "must contain values strictly greater than",
    fixed = TRUE
  )
})


test_that("error_logical errros properly", {
  val <- "FALSE"

  expect_error(
    error_logical(val),
    "needs to be a logical value",
    fixed = TRUE
  )
})


# test_that("error_out_dir works", {
describe("NULL input", {
  it("is treated as valid/skip and returns invisible NULL", {
    expect_null(error_out_dir(NULL))
  })
})

describe("invalid type or length", {
  it("aborts when out_dir is not a character", {
    err <- rlang::catch_cnd(error_out_dir(123))

    expect_true(grepl(
      "Invalid `out_dir` input.",
      conditionMessage(err),
      fixed = TRUE
    ))
    expect_true(grepl(
      "must be a single string",
      conditionMessage(err),
      fixed = TRUE
    ))
  })

  it("aborts when out_dir has length > 1", {
    err <- rlang::catch_cnd(error_out_dir(c("dir_1", "dir_2")))

    expect_true(grepl(
      "Invalid `out_dir` input.",
      conditionMessage(err),
      fixed = TRUE
    ))
    expect_true(grepl(
      "must be a single string",
      conditionMessage(err),
      fixed = TRUE
    ))
  })

  it("aborts when out_dir has length 0", {
    err <- rlang::catch_cnd(error_out_dir(character(0)))

    expect_true(grepl(
      "Invalid `out_dir` input.",
      conditionMessage(err),
      fixed = TRUE
    ))
  })
})

describe("non-existent directory", {
  it("aborts with a directory-not-found message", {
    fake_dir <- file.path(tempdir(), "definitely_does_not_exist_xyz")

    err <- rlang::catch_cnd(error_out_dir(fake_dir))

    expect_true(grepl(
      "Output directory not found.",
      conditionMessage(err),
      fixed = TRUE
    ))
    expect_true(grepl(
      "No such directory",
      conditionMessage(err),
      fixed = TRUE
    ))
    expect_true(grepl(
      "Please provide a valid, existing directory path.",
      conditionMessage(err),
      fixed = TRUE
    ))
  })
})

describe("valid, existing directory", {
  it("returns invisible NULL without error", {
    expect_null(error_out_dir(tempdir()))
  })
})
# })
