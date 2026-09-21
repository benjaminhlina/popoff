# Tests for .transpose_real() / .transpose_back()
#
# These implement a logit-style reparameterization used to keep a bounded
# parameter unconstrained during optimization:
#   real  = log((x - lo) / (hi - x))
#   x     = exp(real) * (hi - lo) / (1 + exp(real)) + lo
#
# .transpose_back() is the mathematical inverse of .transpose_real() over
# the open interval (lo, hi).

describe(".transpose_real()", {
  it("returns 0 at the midpoint of the range", {
    expect_equal(.transpose_real(5, c(0, 10)), 0)
  })

  it("returns a negative value below the midpoint", {
    expect_lt(.transpose_real(2, c(0, 10)), 0)
  })

  it("returns a positive value above the midpoint", {
    expect_gt(.transpose_real(8, c(0, 10)), 0)
  })

  it("matches a hand-computed value", {
    # log((3 - 0) / (10 - 3)) = log(3 / 7)
    expect_equal(.transpose_real(3, c(0, 10)), log(3 / 7))
  })

  it("is vectorized over parameter", {
    result <- .transpose_real(c(2, 5, 8), c(0, 10))
    expect_equal(result, log(c(2, 5, 8) / (10 - c(2, 5, 8))))
  })

  it("returns -Inf / Inf at the lower / upper bound", {
    expect_equal(.transpose_real(0, c(0, 10)), -Inf)
    expect_equal(.transpose_real(10, c(0, 10)), Inf)
  })
})

describe(".transpose_back()", {
  it("returns the midpoint of the range at 0", {
    expect_equal(.transpose_back(0, c(0, 10)), 5)
  })

  it("stays within the parameter range for finite input", {
    result <- .transpose_back(c(-5, -1, 0, 1, 5), c(0, 10))
    expect_true(all(result > 0 & result < 10))
  })

  it("is vectorized over parameter_real", {
    real <- c(-1, 0, 1)
    result <- .transpose_back(real, c(0, 10))
    expected <- exp(real) * 10 / (1 + exp(real)) + 0
    expect_equal(result, expected)
  })

  it("approaches the bounds as parameter_real goes to +/-Inf", {
    expect_equal(.transpose_back(-Inf, c(0, 10)), 0)
    expect_equal(.transpose_back(Inf, c(0, NaN)), NaN)
  })

  it("respects a non-zero lower bound", {
    expect_equal(.transpose_back(0, c(2, 12)), 7)
  })
})

describe(".transpose_real() and .transpose_back() round-trip", {
  it("recovers the original parameter within the open interval", {
    values <- c(0.1, 1, 2.5, 5, 7.5, 9.9)
    range <- c(0, 10)

    recovered <- .transpose_back(.transpose_real(values, range), range)

    expect_equal(recovered, values, tolerance = 1e-8)
  })

  it("recovers the original parameter_real for a range of reals", {
    reals <- c(-3, -1, 0, 1, 3)
    range <- c(0, 10)

    recovered <- .transpose_real(.transpose_back(reals, range), range)

    expect_equal(recovered, reals, tolerance = 1e-8)
  })
})
