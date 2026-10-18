# tests/testthat/test-eval_sigma_iterate.R
#
# Scope: the control flow of eval_sigma_iterate() only (daemon handling, arg
# assembly, sequential vs. parallel dispatch, result binding).
# eval_sigma(), msg_start()/msg_end(), mori and mirai are mocked so no real
# daemons are started and no likelihood maths is exercised here.

# ---- helpers ---------------------------------------------------------------

# Stand-in for eval_sigma(): one row per sigma, records what it was called with.
make_fake_eval_sigma <- function(seen = new.env()) {
  seen$calls <- list()
  fn <- function(sigma, likelihood_array, ..., progress) {
    seen$calls[[length(seen$calls) + 1]] <- list(
      sigma = sigma,
      likelihood_array = likelihood_array,
      progress = progress,
      dots = list(...)
    )
    dplyr::tibble(sigma = sigma, nll = sigma^2)
  }
  list(fn = fn, seen = seen)
}

# Eager, in-process stand-in for mirai::mirai_map(): returns a plain list so
# both `results[]` and `results[TRUE]` work, and records the .args it received.
make_fake_mirai_map <- function(rec = new.env()) {
  fn <- function(.x, .f, ..., .args = list()) {
    rec$x <- .x
    rec$args <- .args
    lapply(.x, function(s) do.call(.f, c(list(s), .args)))
  }
  list(fn = fn, rec = rec)
}

lik <- array(runif(2 * 2 * 3), dim = c(2, 2, 3))
sigmas <- c(0.1, 0.5, 1)

describe("eval_sigma_iterate() sequential path (parallel = FALSE)", {
  it("returns one row per sigma, in the order of sigma_vec", {
    fake <- make_fake_eval_sigma()
    local_mocked_bindings(
      eval_sigma = fake$fn,
      msg_start = function(...) Sys.time(),
      msg_end = function(...) invisible(NULL)
    )

    out <- eval_sigma_iterate(
      sigma_vec = sigmas,
      likelihood_array = lik,
      parallel = FALSE,
      .progress = FALSE
    )

    expect_s3_class(out, "data.frame")
    expect_equal(nrow(out), length(sigmas))
    expect_equal(out$sigma, sigmas)
    expect_equal(out$nll, sigmas^2)
  })

  it("passes likelihood_array, progress = FALSE and extra args to eval_sigma", {
    fake <- make_fake_eval_sigma()
    local_mocked_bindings(
      eval_sigma = fake$fn,
      msg_start = function(...) Sys.time(),
      msg_end = function(...) invisible(NULL)
    )

    eval_sigma_iterate(
      sigma_vec = sigmas,
      likelihood_array = lik,
      parallel = FALSE,
      extra_opt = "abc",
      .progress = FALSE
    )

    calls <- fake$seen$calls
    expect_length(calls, length(sigmas))
    expect_identical(calls[[1]]$likelihood_array, lik)
    expect_false(calls[[1]]$progress)
    expect_identical(calls[[1]]$dots$extra_opt, "abc")
  })

  it("does not touch mirai or mori", {
    fake <- make_fake_eval_sigma()
    local_mocked_bindings(
      eval_sigma = fake$fn,
      msg_start = function(...) Sys.time(),
      msg_end = function(...) invisible(NULL)
    )
    local_mocked_bindings(
      daemons_set = function(...) stop("mirai should not be used"),
      daemons = function(...) stop("mirai should not be used"),
      mirai_map = function(...) stop("mirai should not be used"),
      .package = "mirai"
    )
    local_mocked_bindings(
      share = function(...) stop("mori should not be used"),
      .package = "mori"
    )

    expect_no_error(
      eval_sigma_iterate(sigmas, lik, parallel = FALSE, .progress = FALSE)
    )
  })

  it("returns an empty data frame for an empty sigma_vec", {
    fake <- make_fake_eval_sigma()
    local_mocked_bindings(
      eval_sigma = fake$fn,
      msg_start = function(...) Sys.time(),
      msg_end = function(...) invisible(NULL)
    )

    out <- eval_sigma_iterate(
      numeric(0),
      lik,
      parallel = FALSE,
      .progress = FALSE
    )

    expect_s3_class(out, "data.frame")
    expect_equal(nrow(out), 0)
  })

  it("calls msg_start() once and msg_end() once with the t0 it returned", {
    fake <- make_fake_eval_sigma()
    sentinel <- as.POSIXct("2026-01-01 00:00:00", tz = "UTC")
    msg_log <- new.env()
    msg_log$start <- 0
    msg_log$end_t0 <- NULL
    msg_log$end <- 0

    local_mocked_bindings(
      eval_sigma = fake$fn,
      msg_start = function(...) {
        msg_log$start <- msg_log$start + 1
        sentinel
      },
      msg_end = function(..., t0) {
        msg_log$end <- msg_log$end + 1
        msg_log$end_t0 <- t0
        invisible(NULL)
      }
    )

    eval_sigma_iterate(sigmas, lik, parallel = FALSE, .progress = FALSE)

    expect_equal(msg_log$start, 1)
    expect_equal(msg_log$end, 1)
    expect_identical(msg_log$end_t0, sentinel)
  })
})

describe("eval_sigma_iterate() daemon handling (parallel = TRUE)", {
  it("aborts when no daemons are set and workers is NULL", {
    local_mocked_bindings(
      msg_start = function(...) Sys.time(),
      msg_end = function(...) invisible(NULL)
    )
    local_mocked_bindings(daemons_set = function(...) FALSE, .package = "mirai")

    err <- rlang::catch_cnd(
      eval_sigma_iterate(sigmas, lik, parallel = TRUE, workers = NULL),
      classes = "rlang_error"
    )

    expect_s3_class(err, "rlang_error")
    msg <- conditionMessage(err)
    expect_true(grepl("no mirai daemons are set", msg, fixed = TRUE))
    expect_true(grepl("workers", msg, fixed = TRUE))
  })

  it("aborts when workers < 2", {
    local_mocked_bindings(
      msg_start = function(...) Sys.time(),
      msg_end = function(...) invisible(NULL)
    )
    local_mocked_bindings(daemons_set = function(...) FALSE, .package = "mirai")

    err <- rlang::catch_cnd(
      eval_sigma_iterate(sigmas, lik, parallel = TRUE, workers = 1),
      classes = "rlang_error"
    )

    expect_s3_class(err, "rlang_error")
    expect_true(grepl(
      "no mirai daemons are set",
      conditionMessage(err),
      fixed = TRUE
    ))
  })

  it("starts daemons with `workers` and shuts them down (daemons(0)) on exit", {
    fake <- make_fake_eval_sigma()
    fake_map <- make_fake_mirai_map()
    d_calls <- list()

    local_mocked_bindings(
      eval_sigma = fake$fn,
      msg_start = function(...) Sys.time(),
      msg_end = function(...) invisible(NULL)
    )
    local_mocked_bindings(
      daemons_set = function(...) FALSE,
      daemons = function(n, ...) {
        d_calls[[length(d_calls) + 1]] <<- n
        invisible(NULL)
      },
      mirai_map = fake_map$fn,
      .package = "mirai"
    )
    local_mocked_bindings(share = function(x, ...) x, .package = "mori")

    eval_sigma_iterate(
      sigmas,
      lik,
      parallel = TRUE,
      workers = 3,
      .progress = FALSE
    )

    expect_equal(d_calls, list(3, 0))
  })

  it("reuses existing daemons: informs, ignores `workers`, never starts/stops daemons", {
    fake <- make_fake_eval_sigma()
    fake_map <- make_fake_mirai_map()
    d_calls <- 0

    local_mocked_bindings(
      eval_sigma = fake$fn,
      msg_start = function(...) Sys.time(),
      msg_end = function(...) invisible(NULL)
    )
    local_mocked_bindings(
      daemons_set = function(...) TRUE,
      status = function(...) list(connections = 4L),
      daemons = function(...) {
        d_calls <<- d_calls + 1
        invisible(NULL)
      },
      mirai_map = fake_map$fn,
      .package = "mirai"
    )
    local_mocked_bindings(share = function(x, ...) x, .package = "mori")

    expect_message(
      out <- eval_sigma_iterate(
        sigmas,
        lik,
        parallel = TRUE,
        workers = 99,
        .progress = FALSE
      ),
      "Daemons already running",
      fixed = TRUE
    )

    expect_equal(d_calls, 0)
    expect_equal(nrow(out), length(sigmas))
  })

  it("does not require `workers` when daemons are already set", {
    fake <- make_fake_eval_sigma()
    fake_map <- make_fake_mirai_map()

    local_mocked_bindings(
      eval_sigma = fake$fn,
      msg_start = function(...) Sys.time(),
      msg_end = function(...) invisible(NULL)
    )
    local_mocked_bindings(
      daemons_set = function(...) TRUE,
      status = function(...) list(connections = 2L),
      mirai_map = fake_map$fn,
      .package = "mirai"
    )
    local_mocked_bindings(share = function(x, ...) x, .package = "mori")

    expect_no_error(
      suppressMessages(
        eval_sigma_iterate(
          sigmas,
          lik,
          parallel = TRUE,
          workers = NULL,
          .progress = FALSE
        )
      )
    )
  })
})

describe("eval_sigma_iterate() parallel dispatch and arg assembly", {
  # Shared setup: daemons "already set" so no start/stop logic is involved.
  setup_parallel <- function(env = parent.frame()) {
    fake <- make_fake_eval_sigma()
    fake_map <- make_fake_mirai_map()
    shared_sentinel <- structure(list(), class = "fake_shared")

    local_mocked_bindings(
      eval_sigma = fake$fn,
      msg_start = function(...) Sys.time(),
      msg_end = function(...) invisible(NULL),
      .env = env
    )
    local_mocked_bindings(
      daemons_set = function(...) TRUE,
      status = function(...) list(connections = 2L),
      mirai_map = fake_map$fn,
      .package = "mirai",
      .env = env
    )
    local_mocked_bindings(
      share = function(x, ...) shared_sentinel,
      .package = "mori",
      .env = env
    )
    list(fake = fake, fake_map = fake_map, shared = shared_sentinel)
  }

  it("shares the likelihood array via mori and passes the shared object to the workers", {
    s <- setup_parallel()

    suppressMessages(
      eval_sigma_iterate(sigmas, lik, parallel = TRUE, .progress = FALSE)
    )

    expect_identical(s$fake_map$rec$args$likelihood_array, s$shared)
    expect_false(identical(s$fake_map$rec$args$likelihood_array, lik))
  })

  it("maps over sigma_vec and forces progress = FALSE on the workers", {
    s <- setup_parallel()

    suppressMessages(
      eval_sigma_iterate(sigmas, lik, parallel = TRUE, .progress = FALSE)
    )

    expect_equal(s$fake_map$rec$x, sigmas)
    expect_false(s$fake_map$rec$args$progress)
  })

  it("forwards extra arguments in ... to the workers", {
    s <- setup_parallel()

    suppressMessages(
      eval_sigma_iterate(
        sigmas,
        lik,
        parallel = TRUE,
        extra_opt = "abc",
        another = 42,
        .progress = FALSE
      )
    )

    expect_identical(s$fake_map$rec$args$extra_opt, "abc")
    expect_identical(s$fake_map$rec$args$another, 42)
  })

  it("binds the per-sigma results into one data frame in sigma_vec order", {
    s <- setup_parallel()

    out <- suppressMessages(
      eval_sigma_iterate(sigmas, lik, parallel = TRUE, .progress = FALSE)
    )

    expect_s3_class(out, "data.frame")
    expect_equal(out$sigma, sigmas)
    expect_equal(out$nll, sigmas^2)
  })

  it("collects results for both .progress = TRUE and .progress = FALSE", {
    # NOTE: mirai_map is mocked to return a plain list, so this checks that
    # both branches collect and bind; it does not exercise mirai's own
    # progress bar. See the .progress note in the review comments.
    s <- setup_parallel()

    out_prog <- suppressMessages(
      eval_sigma_iterate(sigmas, lik, parallel = TRUE, .progress = TRUE)
    )
    out_quiet <- suppressMessages(
      eval_sigma_iterate(sigmas, lik, parallel = TRUE, .progress = FALSE)
    )

    expect_equal(out_prog, out_quiet)
  })
})

# ---- optional follow-up (not included to keep scope tight) -----------------
# A real-daemon integration test: mirai::daemons(2); call eval_sigma_iterate()
# on a tiny real likelihood array with the real eval_sigma(); compare against
# the parallel = FALSE result. Needs skip_on_cran() and skip_if_not_installed()
# for mirai/mori, and popr installed so the daemons can resolve eval_sigma.
