describe("gaussian_kernel() kernel size", {
  # size = ceiling(sigma * 4), bumped to the next odd number
  it("returns a square matrix", {
    kern <- gaussian_kernel(sigma = 1)
    expect_true(is.matrix(kern))
    expect_equal(nrow(kern), ncol(kern))
  })

  it("derives size from sigma and forces it to be odd", {
    # sigma = 1   -> ceiling(4)   = 4  -> even -> 5
    expect_equal(dim(gaussian_kernel(sigma = 1)), c(5, 5))
    # sigma = 2   -> ceiling(8)   = 8  -> even -> 9
    expect_equal(dim(gaussian_kernel(sigma = 2)), c(9, 9))
    # sigma = 0.5 -> ceiling(2)   = 2  -> even -> 3
    expect_equal(dim(gaussian_kernel(sigma = 0.5)), c(3, 3))
    # sigma = 1.2 -> ceiling(4.8) = 5  -> already odd -> 5
    expect_equal(dim(gaussian_kernel(sigma = 1.2)), c(5, 5))
  })

  it("always returns an odd number of rows and columns", {
    for (s in c(0.3, 0.7, 1, 1.5, 2, 2.6, 3)) {
      kern <- gaussian_kernel(sigma = s)
      expect_equal(nrow(kern) %% 2, 1, info = paste("sigma =", s))
    }
  })

  it("returns a 1x1 kernel for very small sigma", {
    # sigma = 0.1 -> ceiling(0.4) = 1 -> odd -> 1
    kern <- gaussian_kernel(sigma = 0.1)
    expect_equal(dim(kern), c(1, 1))
    expect_gt(kern[1, 1], 0)
  })

  it("grows with sigma", {
    expect_lt(
      nrow(gaussian_kernel(sigma = 0.5)),
      nrow(gaussian_kernel(sigma = 2))
    )
  })
})

describe("gaussian_kernel() shape", {
  it("peaks in the centre cell when muadv = 0", {
    kern <- gaussian_kernel(sigma = 2, mask = FALSE)
    centre <- (nrow(kern) + 1) / 2
    expect_equal(which.max(kern), (centre - 1) * nrow(kern) + centre)
  })

  it("is symmetric about the centre when muadv = 0", {
    kern <- gaussian_kernel(sigma = 2, mask = FALSE)
    expect_equal(kern, t(kern))
    expect_equal(kern, kern[nrow(kern):1, ])
    expect_equal(kern, kern[, ncol(kern):1])
  })

  it("has no negative, NA or non-finite values", {
    kern <- gaussian_kernel(sigma = 1.5)
    expect_true(all(is.finite(kern)))
    expect_true(all(kern >= 0))
  })

  it("decreases monotonically away from the centre along the middle row", {
    kern <- gaussian_kernel(sigma = 2, mask = FALSE)
    mid <- (nrow(kern) + 1) / 2
    right_half <- kern[mid, mid:ncol(kern)]
    expect_true(all(diff(right_half) < 0))
  })
})

describe("gaussian_kernel() muadv", {
  it("shifts the peak diagonally when muadv is a single value", {
    # 5x5 kernel, centre is (3, 3); muadv = 1 moves the peak to (4, 4)
    kern <- gaussian_kernel(sigma = 1, muadv = 1, mask = FALSE)
    expect_equal(
      unname(which(kern == max(kern), arr.ind = TRUE)[1, ]),
      c(4, 4)
    )
  })

  it("shifts the peak the other way for negative muadv", {
    kern <- gaussian_kernel(sigma = 1, muadv = -1, mask = FALSE)
    expect_equal(
      unname(which(kern == max(kern), arr.ind = TRUE)[1, ]),
      c(2, 2)
    )
  })

  it("accepts a length-2 muadv for independent x/y drift", {
    # NOTE: assumes error_numeric() permits length > 1. If it does not,
    # this is the place to document that muadv is effectively scalar.
    # fx (rows) uses mu[1], fy (columns) uses mu[2].
    kern <- gaussian_kernel(sigma = 1, muadv = c(1, 0), mask = FALSE)
    expect_equal(
      unname(which(kern == max(kern), arr.ind = TRUE)[1, ]),
      c(4, 3)
    )
  })

  it("does not change the kernel dimensions", {
    expect_equal(
      dim(gaussian_kernel(sigma = 1, muadv = 1)),
      dim(gaussian_kernel(sigma = 1, muadv = 0))
    )
  })
})

describe("gaussian_kernel() normalize", {
  it("sums to 1 when normalize = TRUE and mask = FALSE", {
    kern <- gaussian_kernel(sigma = 2, normalize = TRUE, mask = FALSE)
    expect_equal(sum(kern), 1)
  })

  it("sums to 1 when normalize = TRUE and mask = TRUE", {
    # masking removes mass, so the second normalization is what makes this 1
    kern <- gaussian_kernel(sigma = 2, normalize = TRUE, mask = TRUE)
    expect_equal(sum(kern), 1)
  })

  it("does not sum to 1 when normalize = FALSE", {
    kern <- gaussian_kernel(sigma = 2, normalize = FALSE, mask = FALSE)
    expect_false(isTRUE(all.equal(sum(kern), 1)))
  })

  it("gives the same result whether or not normalize runs before masking", {
    # normalizing is a rescale, so normalize + mask must equal
    # mask then rescale by the remaining sum
    masked <- gaussian_kernel(sigma = 2, normalize = FALSE, mask = TRUE)
    both <- gaussian_kernel(sigma = 2, normalize = TRUE, mask = TRUE)
    expect_equal(both, masked / sum(masked))
  })

  it("preserves the relative shape of the kernel", {
    raw <- gaussian_kernel(sigma = 2, normalize = FALSE, mask = FALSE)
    norm <- gaussian_kernel(sigma = 2, normalize = TRUE, mask = FALSE)
    expect_equal(norm, raw / sum(raw))
  })

  it("normalizes a 1x1 kernel to exactly 1", {
    kern <- gaussian_kernel(sigma = 0.1, normalize = TRUE)
    expect_equal(kern[1, 1], 1)
  })
})

describe("gaussian_kernel() mask", {
  it("zeroes the four corners of a 5x5 kernel", {
    # corner distance sqrt(8) ~ 2.83 rounds to 3, which is > radius 2
    kern <- gaussian_kernel(sigma = 1, mask = TRUE)
    corners <- c(kern[1, 1], kern[1, 5], kern[5, 1], kern[5, 5])
    expect_equal(corners, rep(0, 4))
  })

  it("keeps the corners when mask = FALSE", {
    kern <- gaussian_kernel(sigma = 1, mask = FALSE)
    expect_true(all(c(kern[1, 1], kern[1, 5], kern[5, 1], kern[5, 5]) > 0))
  })

  it("keeps cells at (+-2, +-1) from the centre in a 5x5 kernel", {
    # distance sqrt(5) ~ 2.24 rounds to 2, which is not > radius 2
    kern <- gaussian_kernel(sigma = 1, mask = TRUE)
    expect_gt(kern[1, 2], 0)
    expect_gt(kern[2, 1], 0)
  })

  it("does not mask anything in a 3x3 kernel", {
    # corner distance sqrt(2) ~ 1.41 rounds to 1, which is not > radius 1
    expect_equal(
      gaussian_kernel(sigma = 0.5, mask = TRUE),
      gaussian_kernel(sigma = 0.5, mask = FALSE)
    )
  })

  it("leaves the unmasked cells unchanged when normalize = FALSE", {
    masked <- gaussian_kernel(sigma = 2, mask = TRUE)
    unmasked <- gaussian_kernel(sigma = 2, mask = FALSE)
    keep <- masked != 0
    expect_equal(masked[keep], unmasked[keep])
  })

  it("removes mass relative to the unmasked kernel", {
    masked <- gaussian_kernel(sigma = 2, mask = TRUE)
    unmasked <- gaussian_kernel(sigma = 2, mask = FALSE)
    expect_lt(sum(masked), sum(unmasked))
  })

  it("produces a circular footprint: every non-zero cell is within the radius", {
    kern <- gaussian_kernel(sigma = 2, mask = TRUE)
    n <- nrow(kern)
    centre <- (n + 1) / 2
    radius <- trunc(n / 2)
    idx <- which(kern != 0, arr.ind = TRUE)
    dist <- sqrt((idx[, 1] - centre)^2 + (idx[, 2] - centre)^2)
    expect_true(all(round(dist) <= radius))
  })

  it("stays symmetric after masking", {
    kern <- gaussian_kernel(sigma = 2, mask = TRUE)
    expect_equal(kern, t(kern))
    expect_equal(kern, kern[nrow(kern):1, ])
  })
})

describe("gaussian_kernel() back_transform", {
  it("back-transforms sigma before computing the kernel size", {
    # raw sigma = 0.5 would give a 3x3 kernel; the mocked back-transform
    # returns 2, which should give 9x9. This is the regression test for
    # the old bug where the result was assigned to `parameter` and ignored.
    local_mocked_bindings(
      error_transform = function(...) invisible(NULL),
      .transpose_back = function(parameter_real, parameter_range) 2
    )
    kern <- gaussian_kernel(
      sigma = 0.5,
      back_transform = TRUE,
      sigma_range = c(0.1, 3)
    )
    expect_equal(dim(kern), c(9, 9))
  })

  it("matches a direct call using the back-transformed sigma", {
    local_mocked_bindings(
      error_transform = function(...) invisible(NULL),
      .transpose_back = function(parameter_real, parameter_range) 2
    )
    via_transform <- gaussian_kernel(
      sigma = 0.5,
      back_transform = TRUE,
      sigma_range = c(0.1, 3)
    )
    direct <- gaussian_kernel(sigma = 2)
    expect_equal(via_transform, direct)
  })

  it("passes sigma and sigma_range through to .transpose_back()", {
    captured <- NULL
    local_mocked_bindings(
      error_transform = function(...) invisible(NULL),
      .transpose_back = function(parameter_real, parameter_range) {
        captured <<- list(real = parameter_real, range = parameter_range)
        1
      }
    )
    gaussian_kernel(
      sigma = 0.25,
      back_transform = TRUE,
      sigma_range = c(0.1, 3)
    )
    expect_equal(captured$real, 0.25)
    expect_equal(captured$range, c(0.1, 3))
  })

  it("does not call .transpose_back() when back_transform = FALSE", {
    local_mocked_bindings(
      .transpose_back = function(...) stop("should not be called")
    )
    expect_no_error(gaussian_kernel(sigma = 1, back_transform = FALSE))
  })
})

describe("gaussian_kernel() input validation", {
  it("errors on a non-numeric sigma", {
    expect_error(gaussian_kernel(sigma = "a"))
  })

  it("errors on a non-numeric muadv", {
    expect_error(gaussian_kernel(sigma = 1, muadv = "a"))
  })

  it("errors on a non-logical normalize", {
    expect_error(gaussian_kernel(sigma = 1, normalize = "yes"))
  })

  it("errors on a non-logical mask", {
    expect_error(gaussian_kernel(sigma = 1, mask = "yes"))
  })

  it("errors on a non-logical back_transform", {
    expect_error(gaussian_kernel(sigma = 1, back_transform = "yes"))
  })

  it("errors when back_transform = TRUE without a sigma_range", {
    expect_error(
      gaussian_kernel(sigma = 1, back_transform = TRUE, sigma_range = NULL)
    )
  })

  it("no longer accepts kernal_size", {
    # the argument was removed because it was always overwritten
    expect_error(gaussian_kernel(kernal_size = 5, sigma = 1), "unused argument")
  })
})
