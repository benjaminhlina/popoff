# Tests for gaussian_kernel()
#
# Assumes the package already provides `error_numeric()` and `error_logical()`
# input-validators (used unchanged inside gaussian_kernel()); these tests
# exercise the behavior of gaussian_kernel() itself.

test_that("returns a square numeric matrix sized by round(kernal_size)", {
  k <- gaussian_kernel(5, 1)
  expect_true(is.matrix(k))
  expect_true(is.numeric(k))
  expect_equal(dim(k), c(5, 5))

  # non-integer kernal_size is rounded before building x
  expect_equal(dim(gaussian_kernel(4.6, 1, mask = FALSE)), c(5, 5))
  expect_equal(dim(gaussian_kernel(4.4, 1, mask = FALSE)), c(4, 4))
})

test_that("kernal_size below 1 is clamped to 1 (returns 1x1 matrix of 1)", {
  expect_equal(gaussian_kernel(0, 1, normalize = TRUE), matrix(1), )
  expect_equal(gaussian_kernel(-3, 1, normalize = TRUE), matrix(1))
  expect_equal(gaussian_kernel(1, 1, normalize = TRUE), matrix(1))

  # with the normalize = FALSE default, it's just the raw density value at distance 0
  expect_equal(dim(gaussian_kernel(0, 1)), c(1, 1))
  expect_equal(gaussian_kernel(0, 1)[1, 1], gaussian_kernel(1, 1)[1, 1])
})


test_that("mask = TRUE zeroes cells outside the circular radius, regardless of normalize", {
  for (norm_flag in c(FALSE, TRUE)) {
    k <- gaussian_kernel(5, 1, normalize = norm_flag)
    expect_equal(k[1, 1], 0)
    expect_equal(k[1, 5], 0)
    expect_equal(k[5, 1], 0)
    expect_equal(k[5, 5], 0)
    expect_true(k[3, 3] > 0)
    expect_true(k[2, 3] > 0)
  }
})
test_that("default arguments (norm = FALSE, mask = TRUE) sum to 1 and are symmetric", {
  k <- gaussian_kernel(5, 1)
  expect_equal(sum(k), 0.97, tolerance = 0.001)
  expect_equal(k, t(k), tolerance = 0.001)
})

test_that("mask = TRUE zeroes cells outside the circular radius, regardless of normalize", {
  for (norm_flag in c(FALSE, TRUE)) {
    k <- gaussian_kernel(5, 1, normalize = norm_flag)
    expect_equal(k[1, 1], 0)
    expect_equal(k[1, 5], 0)
    expect_equal(k[5, 1], 0)
    expect_equal(k[5, 5], 0)
    expect_true(k[3, 3] > 0)
    expect_true(k[2, 3] > 0)
  }
})

test_that("mask = TRUE when normalizatation = TRUE renormalizes so the kernel sums to 1", {
  expect_equal(
    sum(gaussian_kernel(5, 1, normalize = TRUE)),
    1,
    tolerance = 1e-8
  )
  expect_equal(
    sum(gaussian_kernel(4, 1, normalize = TRUE)),
    1,
    tolerance = 1e-8
  )
  expect_equal(
    sum(gaussian_kernel(5, 100, normalize = TRUE)),
    1,
    tolerance = 1e-8
  )
  expect_equal(
    sum(gaussian_kernel(5, 1, normalize = TRUE)),
    1,
    tolerance = 1e-8
  )
})

test_that("mask = FALSE, norm = FALSE returns the raw (unnormalized) product of 1D densities", {
  k <- gaussian_kernel(5, 1, mask = FALSE, norm = FALSE)
  expect_equal(sum(k), 0.9818148, tolerance = 1e-6)
  # unmasked kernel keeps positive values in the corners
  expect_true(k[1, 1] > 0)
})

test_that("mask = FALSE, norm = TRUE rescales the raw kernel to sum to 1", {
  k <- gaussian_kernel(5, 1, mask = FALSE, norm = TRUE)
  expect_equal(sum(k), 1, tolerance = 1e-8)
  expect_true(k[1, 1] > 0)
})

test_that("peak sits at the (shifted) mean and the kernel stays symmetric under muadv", {
  k0 <- gaussian_kernel(5, 1, mask = FALSE)
  expect_equal(which(k0 == max(k0), arr.ind = TRUE)[1, ], c(row = 3L, col = 3L))

  k_shift <- gaussian_kernel(5, 1, muadv = 1, mask = FALSE)
  expect_equal(
    which(k_shift == max(k_shift), arr.ind = TRUE)[1, ],
    c(row = 4L, col = 4L)
  )
  expect_equal(k_shift, t(k_shift), tolerance = 1e-8)
})

test_that("larger sigma produces a flatter (lower peak-to-edge contrast) kernel", {
  k_tight <- gaussian_kernel(5, 0.5, mask = FALSE)
  k_wide <- gaussian_kernel(5, 100, mask = FALSE)
  ratio_tight <- k_tight[3, 3] / k_tight[1, 1]
  ratio_wide <- k_wide[3, 3] / k_wide[1, 1]
  expect_true(ratio_wide < ratio_tight)
})

test_that("sign of sigma does not matter (only sigma^2 is used)", {
  expect_equal(gaussian_kernel(5, -1), gaussian_kernel(5, 1))
  expect_equal(
    gaussian_kernel(5, -2, mask = FALSE),
    gaussian_kernel(5, 2, mask = FALSE)
  )
})

test_that("sigma = 0 does not error and is handled via the is.finite/is.nan guards", {
  k_unmasked <- gaussian_kernel(5, 0, mask = FALSE)
  expect_true(all(k_unmasked == 0))

  # with mask = TRUE the subsequent 0/0 renormalization yields NaN, not an error
  expect_error(
    gaussian_kernel(5, 0, mask = TRUE),
    regexp = "Gaussian Kernel produces",
    fixed = TRUE
  )
})

test_that("kernal_size, sigma and muadv must be numeric", {
  expect_error(gaussian_kernel("5", 1))
  expect_error(gaussian_kernel(5, "1"))
  expect_error(gaussian_kernel(5, 1, muadv = "0"))
})

test_that("norm and mask must be logical", {
  expect_error(gaussian_kernel(5, 1, norm = "TRUE"))
  expect_error(gaussian_kernel(5, 1, mask = "TRUE"))
  expect_error(gaussian_kernel(5, 1, norm = 1))
  expect_error(gaussian_kernel(5, 1, mask = 1))
})


test_that("masked kernel is symmetric regardless of normalize", {
  expect_equal(
    gaussian_kernel(5, 1, normalize = FALSE),
    t(gaussian_kernel(5, 1, normalize = FALSE)),
    tolerance = 1e-8
  )
  expect_equal(
    gaussian_kernel(5, 1, normalize = TRUE),
    t(gaussian_kernel(5, 1, normalize = TRUE)),
    tolerance = 1e-8
  )
})


test_that("normalize = FALSE never rescales to sum to 1, masked or not (the fixed behavior)", {
  k_mask_false <- gaussian_kernel(5, 1, normalize = FALSE, mask = FALSE)
  expect_equal(sum(k_mask_false), 0.9818148, tolerance = 1e-6)

  # previously this was force-renormalized to 1 despite normalize = FALSE;
  # the switch means it no longer is
  k_mask_true <- gaussian_kernel(5, 1, normalize = FALSE, mask = TRUE)
  expect_equal(sum(k_mask_true), 0.9701547, tolerance = 1e-6)
  expect_false(isTRUE(all.equal(sum(k_mask_true), 1)))
})

test_that("masking + normalizing only rescales magnitude, not the relative shape", {
  k_unnorm <- gaussian_kernel(5, 1, normalize = FALSE, mask = TRUE)
  k_norm <- gaussian_kernel(5, 1, normalize = TRUE, mask = TRUE)
  expect_equal(k_unnorm / sum(k_unnorm), k_norm, tolerance = 1e-8)
})

test_that("default arguments now leave the kernel unnormalized (normalize defaults to FALSE)", {
  k <- gaussian_kernel(5, 1)
  expect_false(isTRUE(all.equal(sum(k), 1)))
  expect_equal(sum(k), 0.9701547, tolerance = 1e-6)
})


test_that("sigma = 0 with mask = FALSE does not error (is.finite guard zeroes it out)", {
  k_unmasked <- gaussian_kernel(5, 0, mask = FALSE)
  expect_true(all(k_unmasked == 0))
})

test_that("sigma = 0 with mask = TRUE errors instead of silently returning NaN", {
  expect_error(gaussian_kernel(5, 0, mask = TRUE))
})


test_that("error_kernel() errors only when mask = TRUE and sigma == 0", {
  expect_error(error_kernel(TRUE, 0))
  expect_silent(error_kernel(TRUE, 1))
  expect_silent(error_kernel(FALSE, 0))
  expect_silent(error_kernel(FALSE, 1))
})

test_that("error_kernel() does not fire for sigma values merely close to zero", {
  expect_silent(error_kernel(TRUE, 1e-8))
  expect_no_error(gaussian_kernel(5, 1e-8))
})
