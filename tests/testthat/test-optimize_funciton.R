# Tests for optimize_sigma()
#
# eval_sigma_iterate() is mocked with a quadratic profile whose minimum sits at
# `target`, so the refinement window logic can be checked without running the
# real (parallel) likelihood evaluation. error_*() validators are mocked as
# no-ops except where validation itself is under test.

# ---- helpers ---------------------------------------------------------------

lik <- array(0, dim = c(3, 3, 2))

make_input <- function(sigma = seq(0.1, 3, by = 0.1), target = 1.234) {
  data.frame(sigma = sigma, neg_log_lik = (sigma - target)^2)
}

# Records every call to eval_sigma_iterate() and returns a quadratic profile
make_eval_mock <- function(target = 1.234) {
  rec <- new.env(parent = emptyenv())
  rec$calls <- list()
  rec$fn <- function(sigma_vec, ...) {
    rec$calls[[length(rec$calls) + 1]] <- list(
      sigma_vec = sigma_vec,
      dots = list(...)
    )
    data.frame(sigma = sigma_vec, neg_log_lik = (sigma_vec - target)^2)
  }
  rec
}

mock_validators <- function(.env = parent.frame()) {
  testthat::local_mocked_bindings(
    error_numeric = function(...) invisible(NULL),
    error_logical = function(...) invisible(NULL),
    .env = .env
  )
}

# ---- output structure ------------------------------------------------------

describe("optimize_sigma() output", {
  it("returns a one-row data.frame with optimal_sigma and neg_log_lik", {
    m <- make_eval_mock()
    mock_validators()
    local_mocked_bindings(eval_sigma_iterate = m$fn)

    out <- optimize_sigma(make_input(), likelihood_array = lik)

    expect_s3_class(out, "data.frame")
    expect_equal(nrow(out), 1L)
    expect_named(out, c("optimal_sigma", "neg_log_lik"))
  })

  it("reports the minimum of the final iteration, with a consistent neg_log_lik", {
    m <- make_eval_mock(target = 1.234)
    mock_validators()
    local_mocked_bindings(eval_sigma_iterate = m$fn)

    out <- optimize_sigma(make_input(), likelihood_array = lik)

    final <- m$fn(m$calls[[length(m$calls)]]$sigma_vec)
    expect_equal(out$neg_log_lik, min(final$neg_log_lik))
    expect_equal(out$optimal_sigma, final$sigma[which.min(final$neg_log_lik)])
    expect_equal(out$neg_log_lik, (out$optimal_sigma - 1.234)^2)
  })

  it("converges close to the true optimum with default settings", {
    m <- make_eval_mock(target = 1.234)
    mock_validators()
    local_mocked_bindings(eval_sigma_iterate = m$fn)

    out <- optimize_sigma(make_input(), likelihood_array = lik)

    expect_equal(out$optimal_sigma, 1.234, tolerance = 1e-3)
  })

  it("returns the best row of the input when n_interval = 0", {
    m <- make_eval_mock()
    mock_validators()
    local_mocked_bindings(eval_sigma_iterate = m$fn)

    out <- optimize_sigma(make_input(), likelihood_array = lik, n_interval = 0)

    expect_length(m$calls, 0L)
    expect_equal(out$optimal_sigma, 1.2)
  })
})

# ---- refinement loop -------------------------------------------------------

describe("optimize_sigma() refinement loop", {
  it("calls eval_sigma_iterate() once per interval", {
    for (n in c(1L, 2L, 4L)) {
      m <- make_eval_mock()
      mock_validators()
      local_mocked_bindings(eval_sigma_iterate = m$fn)

      optimize_sigma(make_input(), likelihood_array = lik, n_interval = n)

      expect_length(m$calls, n)
    }
  })

  it("evaluates n_points sigmas on every iteration", {
    m <- make_eval_mock()
    mock_validators()
    local_mocked_bindings(eval_sigma_iterate = m$fn)

    optimize_sigma(
      make_input(),
      likelihood_array = lik,
      n_points = 9,
      n_interval = 3
    )

    expect_equal(
      vapply(m$calls, function(x) length(x$sigma_vec), integer(1)),
      rep(9L, 3)
    )
  })

  it("centres the first window on the best input sigma, +/- one grid step", {
    m <- make_eval_mock(target = 1.234)
    mock_validators()
    local_mocked_bindings(eval_sigma_iterate = m$fn)

    optimize_sigma(make_input(), likelihood_array = lik)

    # input best = 1.2, input spacing = 0.1
    expect_equal(range(m$calls[[1]]$sigma_vec), c(1.1, 1.3))
  })

  it("shrinks the window using the spacing of the previous result", {
    m <- make_eval_mock(target = 1.234)
    mock_validators()
    local_mocked_bindings(eval_sigma_iterate = m$fn)

    optimize_sigma(make_input(), likelihood_array = lik)

    # iteration 1: 17 points over [1.1, 1.3] -> spacing 0.0125, best = 1.2375
    # iteration 2: best +/- 0.0125
    expect_equal(range(m$calls[[2]]$sigma_vec), c(1.225, 1.25))
    width <- vapply(m$calls, function(x) diff(range(x$sigma_vec)), numeric(1))
    expect_true(all(diff(width) < 0))
  })

  it("clips the window at 0 on the lower side", {
    m <- make_eval_mock(target = -1)
    mock_validators()
    local_mocked_bindings(eval_sigma_iterate = m$fn)

    optimize_sigma(
      make_input(sigma = seq(0, 3, by = 0.5), target = -1),
      likelihood_array = lik,
      n_interval = 1
    )

    expect_equal(range(m$calls[[1]]$sigma_vec), c(0, 0.5))
    expect_true(all(m$calls[[1]]$sigma_vec >= 0))
  })

  it("clips the window at the maximum input sigma on the upper side", {
    m <- make_eval_mock(target = 10)
    mock_validators()
    local_mocked_bindings(eval_sigma_iterate = m$fn)

    optimize_sigma(
      make_input(sigma = seq(0.5, 3, by = 0.5), target = 10),
      likelihood_array = lik,
      n_interval = 3
    )

    all_sigma <- unlist(lapply(m$calls, `[[`, "sigma_vec"))
    expect_lte(max(all_sigma), 3)
    expect_equal(max(m$calls[[1]]$sigma_vec), 3)
  })

  it("picks the first row when neg_log_lik is tied", {
    m <- make_eval_mock()
    mock_validators()
    local_mocked_bindings(eval_sigma_iterate = m$fn)

    flat <- data.frame(sigma = seq(1, 2, by = 0.25), neg_log_lik = 5)
    optimize_sigma(flat, likelihood_array = lik, n_interval = 1)

    # best = 1 (first row), step = 0.25 -> window [0.75, 1.25]
    expect_equal(range(m$calls[[1]]$sigma_vec), c(0.75, 1.25))
  })

  it("ignores NA neg_log_lik values when choosing the best sigma", {
    m <- make_eval_mock()
    mock_validators()
    local_mocked_bindings(eval_sigma_iterate = m$fn)

    input <- make_input()
    input$neg_log_lik[input$sigma == 1.2] <- NA

    optimize_sigma(input, likelihood_array = lik, n_interval = 1)

    # 1.3 is now the lowest non-NA value
    expect_equal(range(m$calls[[1]]$sigma_vec), c(1.1, 1.3))
  })
})

# ---- argument forwarding ---------------------------------------------------

describe("optimize_sigma() argument forwarding", {
  it("passes evaluation settings through to every eval_sigma_iterate() call", {
    m <- make_eval_mock()
    mock_validators()
    local_mocked_bindings(eval_sigma_iterate = m$fn)

    optimize_sigma(
      make_input(),
      likelihood_array = lik,
      n_interval = 2,
      normalize = TRUE,
      mask_likelihood = TRUE,
      workers = 3,
      parallel = FALSE,
      .progress = FALSE
    )

    for (cl in m$calls) {
      expect_identical(cl$dots$likelihood_array, lik)
      expect_true(cl$dots$normalize)
      expect_true(cl$dots$mask_likelihood)
      expect_equal(cl$dots$workers, 3)
      expect_false(cl$dots$parallel)
      expect_false(cl$dots$.progress)
    }
  })
})

# ---- validation ------------------------------------------------------------

describe("optimize_sigma() validation", {
  it("validates each argument with the matching error_*() helper", {
    m <- make_eval_mock()
    numeric_seen <- character()
    logical_seen <- character()
    local_mocked_bindings(
      error_numeric = function(x) {
        numeric_seen <<- c(numeric_seen, rlang::as_label(rlang::enexpr(x)))
        invisible(NULL)
      },
      error_logical = function(x) {
        logical_seen <<- c(logical_seen, rlang::as_label(rlang::enexpr(x)))
        invisible(NULL)
      },
      eval_sigma_iterate = m$fn
    )

    optimize_sigma(make_input(), likelihood_array = lik)

    expect_setequal(numeric_seen, c("n_points", "n_interval"))
    expect_setequal(
      logical_seen,
      c("normalize", "mask_likelihood", "parallel", ".progress")
    )
  })

  it("errors on invalid arguments before any likelihood evaluation", {
    # sentinel: if validation lets bad input through, eval is reached
    local_mocked_bindings(
      eval_sigma_iterate = function(...) {
        rlang::abort("eval reached", class = "eval_reached")
      }
    )

    bad <- list(
      n_points = list(n_points = "a"),
      n_interval = list(n_interval = "a"),
      normalize = list(normalize = "yes"),
      mask_likelihood = list(mask_likelihood = "yes"),
      parallel = list(parallel = "yes"),
      .progress = list(.progress = "yes")
    )

    for (nm in names(bad)) {
      err <- rlang::catch_cnd(
        do.call(
          optimize_sigma,
          c(
            list(sigma_iterate_obj = make_input(), likelihood_array = lik),
            bad[[nm]]
          )
        )
      )
      expect_false(is.null(err), info = nm)
      expect_false(inherits(err, "eval_reached"), info = nm)
    }
  })
})

# ---- edge cases not guarded in the function --------------------------------

describe("optimize_sigma() unguarded edge cases", {
  it("[KNOWN BUG] single-row input fails with an opaque seq() error", {
    # step = (max - min) / (nrow - 1) = 0 / 0 = NaN, so seq() receives NaN.
    # TODO: should fail early with an informative cli_abort() instead
    m <- make_eval_mock()
    mock_validators()
    local_mocked_bindings(eval_sigma_iterate = m$fn)

    err <- rlang::catch_cnd(
      optimize_sigma(
        data.frame(sigma = 1, neg_log_lik = 0.5),
        likelihood_array = lik
      )
    )

    expect_false(is.null(err))
    expect_length(m$calls, 0L)
    # post-fix: expect_match(conditionMessage(err), "at least two", fixed = TRUE)
  })

  it("[KNOWN BUG] all-NA neg_log_lik fails with an opaque seq() error", {
    # which.min() returns integer(0), so `best` has length 0.
    # TODO: should fail early with an informative cli_abort() instead
    m <- make_eval_mock()
    mock_validators()
    local_mocked_bindings(eval_sigma_iterate = m$fn)

    err <- rlang::catch_cnd(
      optimize_sigma(
        data.frame(sigma = c(0.5, 1, 1.5), neg_log_lik = NA_real_),
        likelihood_array = lik
      )
    )

    expect_false(is.null(err))
    expect_length(m$calls, 0L)
    # post-fix: expect_match(conditionMessage(err), "non-missing", fixed = TRUE)
  })
})
