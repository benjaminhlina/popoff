# Tests for hmm_backward_smooth()
#
# Strategy
# - Helpers (validators, msg_*, error_smoothing_collapsed) are mocked with
#   local_mocked_bindings() so each test isolates one branch of the smoother.
#   The warning about collapsed steps now lives in msg_fallback() and the abort
#   for an unrecoverable collapse lives in error_smoothing_collapsed(); here we
#   only check that the smoother calls them with the right arguments. Their
#   messages belong in the helpers' own tests.
# - imagine::convolution2D() is mocked to the identity for the algebra tests
#   (so expected values can be computed by hand), and left real for one
#   integration test that checks invariants only.
# - Fixtures are small random arrays; edge cases are built by overwriting
#   specific cells with point masses.
# - The pop-off boundary condition is the LAST slice of likelihood_array, i.e.
#   slice n_time (the likelihood array has the same number of slices as the
#   filtered posterior).

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

  # likelihood_array has n_time slices; the last (slice n_time) is the
  # pop-off slice
  lik <- array(
    stats::runif(nx * ny * n_time, 0.1, 1),
    dim = c(nx, ny, n_time)
  )
  lik[rep(as.vector(land), n_time)] <- NA_real_

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
# records how they were called:
#   fwd_obj, logical  first argument of the two validators
#   start, end        msg_start() / msg_end() arguments
#   progress          list of everything passed to msg_progress()
#   fallback          list of everything passed to msg_fallback()
#   collapsed         one list(fac, time_step) per error_smoothing_collapsed()
#                     call, in call order
stub_helpers <- function(
  env = parent.frame(),
  mock_error_fwd = TRUE,
  mock_error_logical = TRUE,
  mock_error_collapsed = TRUE
) {
  spy <- new.env(parent = emptyenv())
  spy$collapsed <- list()

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
  if (mock_error_collapsed) {
    testthat::local_mocked_bindings(
      error_smoothing_collapsed = function(...) {
        spy$collapsed[[length(spy$collapsed) + 1L]] <- list(...)
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
      args <- list(...)
      spy$progress <- args
      # msg_progress() is now always called; the first argument is the
      # progress flag. cli_progress_update()/_done() need a live bar in the
      # caller's frame, so only start one when progress is on.
      if (isTRUE(args[[1]])) {
        cli::cli_progress_bar(
          name = "test",
          total = args$n_time,
          .envir = parent.frame()
        )
      }
      invisible(NULL)
    },
    msg_fallback = function(...) {
      spy$fallback <- list(...)
      invisible(NULL)
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

  it("forces land to NA even where later filtered slices carry mass on it", {
    land <- matrix(FALSE, 4, 5)
    land[1, 1] <- TRUE
    fwd <- make_fwd(n_time = 3, land = land)
    fwd$normalised_posterior[1, 1, 2:3] <- 0.3 # land is NA only in slice 1
    stub_helpers()
    stub_convolution()

    out <- hmm_backward_smooth(fwd, progress = FALSE)

    expect_true(all(is.na(out[1, 1, ])))
    for (t in seq_len(dim(out)[3])) {
      expect_equal(sum(out[,, t], na.rm = TRUE), 1)
    }
  })
})

describe("hmm_backward_smooth() pop-off boundary condition", {
  it("sets the last slice to the normalised pop-off likelihood, ignoring land", {
    land <- matrix(FALSE, 4, 5)
    land[1, 1] <- TRUE
    fwd <- make_fwd(n_time = 3, land = land)
    fwd$likelihood_array[1, 1, 3] <- 1e6 # land mass must not count
    stub_helpers()
    stub_convolution()

    out <- hmm_backward_smooth(fwd, progress = FALSE)

    pop <- fwd$likelihood_array[,, 3]
    pop[land] <- 0
    expected <- pop / sum(pop)

    expect_equal(out[,, 3][!land], expected[!land])
  })

  it("falls back to the filtered posterior when the pop-off slice has no mass", {
    land <- matrix(FALSE, 4, 5)
    land[1, 1] <- TRUE

    for (empty in list(0, NA_real_)) {
      fwd <- make_fwd(n_time = 3, land = land)
      fwd$likelihood_array[,, 3] <- empty
      spy <- stub_helpers()
      stub_convolution()

      expect_no_warning(out <- hmm_backward_smooth(fwd, progress = FALSE))

      expect_equal(
        out[,, 3][!land],
        fwd$normalised_posterior[,, 3][!land]
      )
      # the pop-off fallback is not a collapsed step
      expect_identical(spy$fallback$n_fallback, 0L)
      expect_length(spy$collapsed, 0)
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

    pop <- fwd$likelihood_array[,, 2]
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

    pop <- fwd$likelihood_array[,, 4]
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
    fwd$likelihood_array[3, 3, 2] <- 0 # smooth = 0 ...
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
  it("uses the filtered posterior and reports one collapsed step", {
    fwd <- make_fwd(nx = 3, ny = 3, n_time = 2)
    fwd$likelihood_array[,, 2] <- point_mass(3, 3, 1, 1)
    fwd$normalised_posterior[,, 1] <- point_mass(3, 3, 3, 3)
    spy <- stub_helpers()
    stub_convolution()

    expect_no_warning(out <- hmm_backward_smooth(fwd, progress = FALSE))

    expect_equal(out[,, 1], point_mass(3, 3, 3, 3))
    expect_identical(spy$fallback$n_fallback, 1L)
    expect_length(spy$collapsed, 1)
    expect_equal(spy$collapsed[[1]]$fac, 1)
    expect_equal(spy$collapsed[[1]]$time_step, 2)
  })

  it("counts every collapsed step and passes each one to the helper", {
    fwd <- make_fwd(nx = 3, ny = 3, n_time = 3)
    fwd$likelihood_array[,, 3] <- point_mass(3, 3, 1, 1)
    fwd$normalised_posterior[,, 2] <- point_mass(3, 3, 3, 3)
    fwd$normalised_posterior[,, 1] <- point_mass(3, 3, 2, 2)
    spy <- stub_helpers()
    stub_convolution()

    out <- hmm_backward_smooth(fwd, progress = FALSE)

    expect_equal(out[,, 2], point_mass(3, 3, 3, 3))
    expect_equal(out[,, 1], point_mass(3, 3, 2, 2))
    expect_identical(spy$fallback$n_fallback, 2L)
    expect_length(spy$collapsed, 2)
    expect_equal(
      vapply(spy$collapsed, function(x) x$time_step, numeric(1)),
      c(3, 2)
    )
    expect_equal(
      vapply(spy$collapsed, function(x) x$fac, numeric(1)),
      c(1, 1)
    )
  })

  it("reports zero collapsed steps when nothing collapses", {
    fwd <- make_fwd(n_time = 3)
    spy <- stub_helpers()
    stub_convolution()

    hmm_backward_smooth(fwd, progress = FALSE)

    expect_identical(spy$fallback$n_fallback, 0L)
    expect_length(spy$collapsed, 0)
  })

  it("stops when error_smoothing_collapsed() aborts because the fallback also has no mass", {
    fwd <- make_fwd(n_time = 3)
    fwd$normalised_posterior[,, 2] <- 0
    spy <- stub_helpers(mock_error_collapsed = FALSE)
    seen <- NULL
    testthat::local_mocked_bindings(
      error_smoothing_collapsed = function(fac, time_step) {
        seen <<- list(fac = fac, time_step = time_step)
        cli::cli_abort("smoothing collapsed")
      }
    )
    stub_convolution()

    err <- rlang::catch_cnd(
      hmm_backward_smooth(fwd, progress = FALSE),
      classes = "error"
    )

    expect_false(is.null(err))
    expect_true(grepl(
      "smoothing collapsed",
      conditionMessage(err),
      fixed = TRUE
    ))
    expect_equal(seen$fac, 0)
    expect_equal(seen$time_step, 3)
    # aborted before the closing messages
    expect_null(spy$fallback)
    expect_null(spy$end)
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

  it("passes progress = TRUE and a bar sized to n_time to msg_progress()", {
    fwd <- make_fwd(n_time = 3)
    spy <- stub_helpers()
    stub_convolution()

    hmm_backward_smooth(fwd, progress = TRUE)

    # first argument is the progress flag
    expect_true(spy$progress[[1]])
    expect_identical(spy$progress$n_time, 3L)
    expect_identical(spy$progress$name, "HMM backwards smoothing")
    expect_identical(spy$progress$type, "Smoothing")
  })

  it("passes progress = FALSE to msg_progress() and still runs", {
    fwd <- make_fwd(n_time = 3)
    spy <- stub_helpers()
    stub_convolution()

    expect_no_error(out <- hmm_backward_smooth(fwd, progress = FALSE))

    # msg_progress() is always called; it decides whether to draw a bar
    expect_false(spy$progress[[1]])
    expect_identical(dim(out), dim(fwd$normalised_posterior))
  })
})

describe("hmm_backward_smooth() with the real convolution", {
  it("returns a valid, land-masked posterior for every slice", {
    land <- matrix(FALSE, 6, 7)
    land[1, 1] <- TRUE
    land[6, 7] <- TRUE
    land[3, 4] <- TRUE
    fwd <- make_fwd(nx = 6, ny = 7, n_time = 5, land = land)
    spy <- stub_helpers()

    expect_no_warning(out <- hmm_backward_smooth(fwd, progress = FALSE))

    for (t in seq_len(dim(out)[3])) {
      s <- out[,, t]
      expect_true(all(is.na(s[land])))
      expect_true(all(is.finite(s[!land])))
      expect_true(all(s[!land] >= 0))
      expect_equal(sum(s, na.rm = TRUE), 1)
    }
    expect_identical(spy$fallback$n_fallback, 0L)
  })
})
