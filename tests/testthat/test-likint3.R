# ----- likint3() -----

describe("likint3()", {
  it("preserves NA positions and fills non-NA positions with the integration result", {
    w <- c(1, NA, 3, NA, 5)
    wsd <- c(0.5, 0.5, 0.5, 0.5, 0.5)

    local_mocked_bindings(
      .integrate_dnorm_rows = function(w, wsd, minT, maxT) rep(0.42, length(w))
    )

    result <- likint3(w, wsd, minT = 0, maxT = 10)

    expect_equal(result, c(0.42, NA, 0.42, NA, 0.42))
  })

  it("calls .integrate_dnorm_rows with only the non-NA subset of w and wsd, and forwards minT/maxT unchanged", {
    w <- c(1, NA, 3)
    wsd <- c(0.1, 0.2, 0.3)
    captured <- NULL

    local_mocked_bindings(
      .integrate_dnorm_rows = function(w, wsd, minT, maxT) {
        captured <<- list(w = w, wsd = wsd, minT = minT, maxT = maxT)
        rep(1, length(w))
      }
    )

    likint3(w, wsd, minT = -5, maxT = 5)

    expect_equal(captured$w, c(1, 3))
    expect_equal(captured$wsd, c(0.1, 0.3))
    expect_equal(captured$minT, -5)
    expect_equal(captured$maxT, 5)
  })

  it("returns an all-NA result, without erroring, when every element of w is NA", {
    w <- c(NA_real_, NA_real_)
    wsd <- c(0.5, 0.5)

    called_with_empty <- FALSE
    local_mocked_bindings(
      .integrate_dnorm_rows = function(w, wsd, minT, maxT) {
        called_with_empty <<- length(w) == 0 && length(wsd) == 0
        numeric(0)
      }
    )

    result <- likint3(w, wsd, minT = 0, maxT = 1)

    expect_true(all(is.na(result)))
    expect_length(result, 2)
    expect_true(called_with_empty)
  })

  it("preserves the dimensions of a matrix input", {
    w <- matrix(c(1, NA, 3, 4), nrow = 2)
    wsd <- matrix(c(0.5, 0.5, 0.5, 0.5), nrow = 2)

    local_mocked_bindings(
      .integrate_dnorm_rows = function(w, wsd, minT, maxT) rep(0.9, length(w))
    )

    result <- likint3(w, wsd, minT = 0, maxT = 1)

    expect_equal(dim(result), c(2, 2))
    expect_true(is.na(result[2, 1]))
    expect_equal(result[1, 1], 0.9)
  })
})
