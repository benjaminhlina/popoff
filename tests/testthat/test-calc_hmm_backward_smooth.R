# Tests for hmm_backward_smooth()
#
# Strategy
# - Helpers (validators, msg_*) are mocked with local_mocked_bindings() so each
#   test isolates one branch of the smoother.
# - imagine::convolution2D() is mocked to the identity for the algebra tests
#   (so expected values can be computed by hand), and left real for one
#   integration test that checks invariants only.
# - Fixtures are small random arrays; edge cases are built by overwriting
#   specific cells with point masses.

# ---- fixtures and helpers ---------------------------------------------------

make_fwd <- function(nx = 4, ny = 5, n_time = 3, land = NULL, seed = 42) {
  set.seed(seed)
  if (is.null(land)) {
    land <- matrix(FALSE, nx, ny)
  }

  norm_array <- function(n) {
    out <- array(NA_real_, dim = c(nx, ny, n))
    for (i in seq_len(n)) {
      m <- matrix(stats::runif(nx * ny, 0.1, 1), nx, ny)
      m[land] <- NA_real_
      out[,, i] <- m / sum(m, na.rm = TRUE)
    }
    out
  }

  dates <- format(seq(as.Date("2024-01-01"), by = "day", length.out = n_time))

  alpha <- norm_array(n_time)
  dimnames(alpha) <- list(NULL, NULL, dates)

  pred <- norm_array(n_time)
  dimnames(pred) <- dimnames(alpha)

  # likelihood_array has n_time + 1 slices; the last is the pop-off slice
  lik <- array(
    stats::runif(nx * ny * (n_time + 1), 0.1, 1),
    dim = c(nx, ny, n_time + 1)
  )
  lik[rep(as.vector(land), n_time + 1)] <- NA_real_

  k <- matrix(c(1, 2, 1, 2, 4, 2, 1, 2, 1), 3, 3)

  list(
    normalised_posterior = alpha,
    pc_predicted = pred,
    likelihood_array = lik,
    kernel = k / sum(k)
  )
}

point_mass <- function(nx, ny, i, j) {
  m <- matrix(0, nx, ny)
  m[i, j] <- 1
  m
}

# One backward step with an identity convolution, computed independently of
# the function under test
backward_step <- function(smooth_t, pred_t, alpha_prev) {
  rat <- smooth_t / pred_t
  rat[!is.finite(rat)] <- 0
  post <- rat * alpha_prev
  post / sum(post, na.rm = TRUE)
}

# Mock the helpers hmm_backward_smooth() calls. Returns an environment that
# records how they were called.
stub_helpers <- function(
  env = parent.frame(),
  mock_error_fwd = TRUE,
  mock_error_logical = TRUE
) {
  spy <- new.env(parent = emptyenv())

  if (mock_error_fwd) {
    testthat::local_mocked_bindings(
      error_forward_filter_obj = function(...) {
        spy$fwd_obj <- list(...)[[1]]
        invisible(NULL)
      },
      .env = env
    )
  }
  if (mock_error_logical) {
    testthat::local_mocked_bindings(
      error_logical = function(...) {
        spy$logical <- list(...)[[1]]
        invisible(NULL)
      },
      .env = env
    )
  }

  testthat::local_mocked_bindings(
    msg_start = function(...) {
      spy$start <- list(...)[[1]]
      "t0-token"
    },
    msg_end = function(...) {
      spy$end <- list(...)
      invisible(NULL)
    },
    msg_progress = function(...) {
      spy$progress <- list(...)
      # cli_progress_update()/_done() need a live bar in the caller's frame
      cli::cli_progress_bar(
        name = "test",
        total = spy$progress$n_time,
        .envir = parent.frame()
      )
    },
    .env = env
  )

  spy
}

stub_convolution <- function(
  fn = function(x, kernel, ...) x,
  env = parent.frame()
) {
  testthat::local_mocked_bindings(
    convolution2D = fn,
    .package = "imagine",
    .env = env
  )
}

# ---- tests ------------------------------------------------------------------

describe("hmm_backward_smooth() input validation", {
  it("validates the filter object and progress before doing any work", {
    fwd <- make_fwd(n_time = 2)
    spy <- stub_helpers()
    stub_convolution()

    hmm_backward_smooth(fwd, progress = FALSE)

    expect_identical(spy$fwd_obj, fwd)
    expect_false(spy$logical)
  })

  it("aborts, before any work, when the filter object is invalid", {
    fwd <- make_fwd(n_time = 2)
    spy <- stub_helpers(mock_error_fwd = FALSE)
    testthat::local_mocked_bindings(
      error_forward_filter_obj = function(...) {
        cli::cli_abort("bad filter object")
      }
    )

    err <- rlang::catch_cnd(hmm_backward_smooth(fwd), classes = "error")

    expect_false(is.null(err))
    expect_true(grepl("bad filter object", conditionMessage(err), fixed = TRUE))
    expect_null(spy$start)
  })

  it("aborts when progress is not a single logical", {
    fwd <- make_fwd(n_time = 2)
    stub_helpers(mock_error_logical = FALSE)

    expect_error(hmm_backward_smooth(fwd, progress = "yes"))
  })

  it("aborts on an object that is not a forward filter result", {
    expect_error(hmm_backward_smooth(list()))
  })
})

describe("hmm_backward_smooth() output structure", {
  it("returns an array with the dim and dimnames of the filtered posterior", {
    fwd <- make_fwd(n_time = 3)
    stub_helpers()
    stub_convolution()

    out <- hmm_backward_smooth(fwd, progress = FALSE)

    expect_true(is.array(out))
    expect_type(out, "double")
    expect_identical(dim(out), dim(fwd$normalised_posterior))
    expect_identical(dimnames(out), dimnames(fwd$normalised_posterior))
  })

  it("keeps land NA in every slice and normalises every slice over water", {
    land <- matrix(FALSE, 4, 5)
    land[1, 1] <- TRUE
    land[4, 5] <- TRUE
    fwd <- make_fwd(n_time = 3, land = land)
    stub_helpers()
    stub_convolution()

    out <- hmm_backward_smooth(fwd, progress = FALSE)

    for (t in seq_len(dim(out)[3])) {
      s <- out[,, t]
      expect_true(all(is.na(s[land])))
      expect_true(all(is.finite(s[!land])))
      expect_true(all(s[!land] >= 0))
      expect_equal(sum(s, na.rm = TRUE), 1)
    }
  })

  it("takes the land mask from the first slice of the filtered posterior only", {
    fwd <- make_fwd(n_time = 3)
    fwd$normalised_posterior[2, 3, 1] <- NA_real_ # NA in slice 1 only
    stub_helpers()
    stub_convolution()

    out <- hmm_backward_smooth(fwd, progress = FALSE)

    expect_true(all(is.na(out[2, 3, ])))
    expect_equal(sum(out[,, 1], na.rm = TRUE), 1)
  })
})

describe("hmm_backward_smooth() pop-off boundary condition", {
  it("sets the last slice to the normalised pop-off likelihood, ignoring land", {
    land <- matrix(FALSE, 4, 5)
    land[1, 1] <- TRUE
    fwd <- make_fwd(n_time = 3, land = land)
    fwd$likelihood_array[1, 1, 4] <- 1e6 # land mass must not count
    stub_helpers()
    stub_convolution()

    out <- hmm_backward_smooth(fwd, progress = FALSE)

    pop <- fwd$likelihood_array[,, 4]
    pop[land] <- 0
    expected <- pop / sum(pop)

    expect_equal(out[,, 3][!land], expected[!land])
  })

  it("falls back to the filtered posterior when the pop-off slice has no mass", {
    land <- matrix(FALSE, 4, 5)
    land[1, 1] <- TRUE

    for (empty in list(0, NA_real_)) {
      fwd <- make_fwd(n_time = 3, land = land)
      fwd$likelihood_array[,, 4] <- empty
      stub_helpers()
      stub_convolution()

      expect_no_warning(out <- hmm_backward_smooth(fwd, progress = FALSE))

      expect_equal(
        out[,, 3][!land],
        fwd$normalised_posterior[,, 3][!land]
      )
    }
  })
})

describe("hmm_backward_smooth() backward recursion", {
  it("matches a hand-computed single step (n_time = 2)", {
    land <- matrix(FALSE, 4, 5)
    land[1, 1] <- TRUE
    land[4, 5] <- TRUE
    fwd <- make_fwd(n_time = 2, land = land)
    stub_helpers()
    stub_convolution()

    out <- hmm_backward_smooth(fwd, progress = FALSE)

    pop <- fwd$likelihood_array[,, 3]
    pop[land] <- 0
    s2 <- pop / sum(pop)
    s1 <- backward_step(
      s2,
      fwd$pc_predicted[,, 2],
      fwd$normalised_posterior[,, 1]
    )

    expect_equal(out[,, 2][!land], s2[!land])
    expect_equal(out[,, 1][!land], s1[!land])
  })

  it("chains steps correctly over several time steps (n_time = 4)", {
    land <- matrix(FALSE, 4, 5)
    land[2, 2] <- TRUE
    fwd <- make_fwd(n_time = 4, land = land)
    stub_helpers()
    stub_convolution()

    out <- hmm_backward_smooth(fwd, progress = FALSE)

    pop <- fwd$likelihood_array[,, 5]
    pop[land] <- 0
    expected <- array(NA_real_, dim = dim(out))
    expected[,, 4] <- pop / sum(pop)
    for (t in 4:2) {
      expected[,, t - 1] <- backward_step(
        expected[,, t],
        fwd$pc_predicted[,, t],
        fwd$normalised_posterior[,, t - 1]
      )
    }

    for (t in 1:4) {
      expect_equal(out[,, t][!land], expected[,, t][!land])
    }
  })

  it("convolves once per backward step, with the kernel from the filter object", {
    fwd <- make_fwd(n_time = 4)
    stub_helpers()
    seen <- list()
    stub_convolution(function(x, kernel, ...) {
      seen[[length(seen) + 1]] <<- kernel
      x
    })

    hmm_backward_smooth(fwd, progress = FALSE)

    expect_length(seen, 3)
    for (k in seen) {
      expect_identical(k, fwd$kernel)
    }
  })
})

describe("hmm_backward_smooth() non-finite handling", {
  it("turns x/0 and 0/0 in the smoothed/predicted ratio into zero mass", {
    fwd <- make_fwd(n_time = 2)
    fwd$pc_predicted[2, 2, 2] <- 0 # smooth > 0, pred = 0  -> Inf
    fwd$likelihood_array[3, 3, 3] <- 0 # smooth = 0 ...
    fwd$pc_predicted[3, 3, 2] <- 0 # ... and pred = 0     -> NaN
    stub_helpers()
    stub_convolution()

    expect_no_warning(out <- hmm_backward_smooth(fwd, progress = FALSE))

    out <- unname(out)
    expect_false(all(is.nan(out)))
    expect_true(all(is.finite(out)))
    expect_equal(out[2, 2, 1], 0)
    expect_equal(out[3, 3, 1], 0)
    expect_equal(sum(out[,, 1]), 1)
  })

  it("zeroes non-finite convolution output on water cells", {
    fwd <- make_fwd(n_time = 2)
    stub_helpers()
    stub_convolution(function(x, kernel, ...) {
      x[2, 2] <- NaN
      x
    })

    out <- hmm_backward_smooth(fwd, progress = FALSE)
    out <- unname(out)
    expect_false(all(is.nan(out)))
    expect_equal(out[2, 2, 1], 0)
    expect_equal(sum(out[,, 1]), 1)
  })
})

describe("hmm_backward_smooth() collapsed steps", {
  it("uses the filtered posterior and warns (singular) when one step collapses", {
    fwd <- make_fwd(nx = 3, ny = 3, n_time = 2)
    fwd$likelihood_array[,, 3] <- point_mass(3, 3, 1, 1)
    fwd$normalised_posterior[,, 1] <- point_mass(3, 3, 3, 3)
    stub_helpers()
    stub_convolution()

    expect_warning(
      out <- hmm_backward_smooth(fwd, progress = FALSE),
      "1 step collapsed during smoothing",
      fixed = TRUE
    )

    expect_equal(out[,, 1], point_mass(3, 3, 3, 3))
  })

  it("counts every collapsed step and pluralises the warning", {
    fwd <- make_fwd(nx = 3, ny = 3, n_time = 3)
    fwd$likelihood_array[,, 4] <- point_mass(3, 3, 1, 1)
    fwd$normalised_posterior[,, 2] <- point_mass(3, 3, 3, 3)
    fwd$normalised_posterior[,, 1] <- point_mass(3, 3, 2, 2)
    stub_helpers()
    stub_convolution()

    expect_warning(
      out <- hmm_backward_smooth(fwd, progress = FALSE),
      "2 steps collapsed during smoothing",
      fixed = TRUE
    )

    expect_equal(out[,, 2], point_mass(3, 3, 3, 3))
    expect_equal(out[,, 1], point_mass(3, 3, 2, 2))
  })

  it("aborts, naming the offending time step, when the fallback also has no mass", {
    fwd <- make_fwd(n_time = 3)
    fwd$normalised_posterior[,, 2] <- 0
    stub_helpers()
    stub_convolution()

    err <- rlang::catch_cnd(
      hmm_backward_smooth(fwd, progress = FALSE),
      classes = "error"
    )
    msg <- conditionMessage(err)

    expect_false(is.null(err))
    expect_true(grepl("Smoothing collapsed at time step 2", msg, fixed = TRUE))
    expect_true(grepl("filtered posterior also sums to 0", msg, fixed = TRUE))
  })
})

describe("hmm_backward_smooth() messaging and progress", {
  it("brackets the run with msg_start()/msg_end() using the same timer", {
    fwd <- make_fwd(n_time = 2)
    spy <- stub_helpers()
    stub_convolution()

    hmm_backward_smooth(fwd, progress = FALSE)

    expect_identical(spy$start, "HMM backwards smoothing")
    expect_identical(spy$end[[1]], "HMM backwards smoothing")
    expect_identical(spy$end[[2]], "t0-token")
  })

  it("starts a progress bar sized to n_time only when progress = TRUE", {
    fwd <- make_fwd(n_time = 3)
    spy <- stub_helpers()
    stub_convolution()

    hmm_backward_smooth(fwd, progress = TRUE)

    expect_identical(spy$progress$n_time, 3L)
    expect_identical(spy$progress$name, "HMM backwards smoothing")
    expect_identical(spy$progress$type, "Smoothing")
  })

  it("does not touch the progress machinery when progress = FALSE", {
    fwd <- make_fwd(n_time = 3)
    spy <- stub_helpers()
    stub_convolution()

    hmm_backward_smooth(fwd, progress = FALSE)

    expect_null(spy$progress)
  })
})

describe("hmm_backward_smooth() with the real convolution", {
  it("returns a valid, land-masked posterior for every slice", {
    land <- matrix(FALSE, 6, 7)
    land[1, 1] <- TRUE
    land[6, 7] <- TRUE
    land[3, 4] <- TRUE
    fwd <- make_fwd(nx = 6, ny = 7, n_time = 5, land = land)
    stub_helpers()

    expect_no_warning(out <- hmm_backward_smooth(fwd, progress = FALSE))

    for (t in seq_len(dim(out)[3])) {
      s <- out[,, t]
      expect_true(all(is.na(s[land])))
      expect_true(all(is.finite(s[!land])))
      expect_true(all(s[!land] >= 0))
      expect_equal(sum(s, na.rm = TRUE), 1)
    }
  })
})
