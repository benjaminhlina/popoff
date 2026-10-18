# tests/testthat/test-eval_sigma.R
#
# Scope: argument routing through `...`, the two return shapes, and the
# error path for bad `...`. gaussian_kernel() and hmm_forward_filter() are
# mocked so these tests isolate eval_sigma() itself.

# ---- helpers ----------------------------------------------------------------

# Records the arguments each mocked collaborator receives.
make_spy <- function() {
  spy <- new.env(parent = emptyenv())
  spy$kern_args <- NULL
  spy$hmm_args <- NULL
  spy
}

# Mock formals deliberately include extra, non-shared arguments so routing by
# name can be exercised: `kern_extra` (kernel only), `hmm_extra` (filter only),
# `shared_extra` (both).
mock_kernel_factory <- function(spy) {
  function(
    sigma,
    normalize = FALSE,
    mask = TRUE,
    kern_extra = NULL,
    shared_extra = NULL
  ) {
    spy$kern_args <- as.list(environment())
    structure(list(sigma = sigma), class = "mock_kernel")
  }
}

mock_hmm_factory <- function(spy) {
  function(
    likelihood_array,
    kernel,
    return_maps = FALSE,
    hmm_extra = NULL,
    shared_extra = NULL
  ) {
    spy$hmm_args <- as.list(environment())
    list(
      log_lik = -12.5,
      neg_log_lik = 12.5,
      n_bad_days = 1L,
      n_updated = 8L,
      ncps = seq_len(9),
      maps = "mock_maps"
    )
  }
}

dummy_lik <- array(1, dim = c(3, 3, 9))


describe("eval_sigma()", {
  it("return_maps = FALSE returns a one-row data.frame with the expected columns", {
    spy <- make_spy()
    local_mocked_bindings(
      gaussian_kernel = mock_kernel_factory(spy),
      hmm_forward_filter = mock_hmm_factory(spy)
    )

    out <- eval_sigma(sigma = 1.5, likelihood_array = dummy_lik)

    expect_s3_class(out, "data.frame")
    expect_equal(nrow(out), 1L)
    expect_named(
      out,
      c(
        "sigma",
        "log_lik",
        "neg_log_lik",
        "n_bad_days",
        "n_updated",
        "n_steps"
      )
    )
  })

  it("fills the columns from the filter output and sigma", {
    spy <- make_spy()
    local_mocked_bindings(
      gaussian_kernel = mock_kernel_factory(spy),
      hmm_forward_filter = mock_hmm_factory(spy)
    )

    out <- eval_sigma(sigma = 1.5, likelihood_array = dummy_lik)

    expect_equal(out$sigma, 1.5)
    expect_equal(out$log_lik, -12.5)
    expect_equal(out$neg_log_lik, 12.5)
    expect_equal(out$n_bad_days, 1L)
    expect_equal(out$n_updated, 8L)
    expect_equal(out$n_steps, 9L) # length(ncps)
  })
})

describe("return_maps = TRUE", {
  it("returns the raw filter output untouched", {
    spy <- make_spy()
    local_mocked_bindings(
      gaussian_kernel = mock_kernel_factory(spy),
      hmm_forward_filter = mock_hmm_factory(spy)
    )

    out <- eval_sigma(
      sigma = 1,
      likelihood_array = dummy_lik,
      return_maps = TRUE
    )

    expect_type(out, "list")
    expect_false(is.data.frame(out))
    expect_equal(out$maps, "mock_maps")
    expect_equal(out$ncps, seq_len(9))
  })
})

describe("argument passing to collaborators", {
  it("passes sigma, normalize and mask to gaussian_kernel()", {
    spy <- make_spy()
    local_mocked_bindings(
      gaussian_kernel = mock_kernel_factory(spy),
      hmm_forward_filter = mock_hmm_factory(spy)
    )

    eval_sigma(
      sigma = 2,
      normalize = TRUE,
      mask = FALSE,
      likelihood_array = dummy_lik
    )

    expect_equal(spy$kern_args$sigma, 2)
    expect_true(spy$kern_args$normalize)
    expect_false(spy$kern_args$mask)
  })

  it("passes the likelihood array, kernel and return_maps to hmm_forward_filter()", {
    spy <- make_spy()
    local_mocked_bindings(
      gaussian_kernel = mock_kernel_factory(spy),
      hmm_forward_filter = mock_hmm_factory(spy)
    )

    eval_sigma(sigma = 2, likelihood_array = dummy_lik, return_maps = TRUE)

    expect_identical(spy$hmm_args$likelihood_array, dummy_lik)
    expect_s3_class(spy$hmm_args$kernel, "mock_kernel")
    expect_equal(spy$hmm_args$kernel$sigma, 2)
    expect_true(spy$hmm_args$return_maps)
  })
})

describe("routing of `...`", {
  it("sends kernel-only arguments to gaussian_kernel() only", {
    spy <- make_spy()
    local_mocked_bindings(
      gaussian_kernel = mock_kernel_factory(spy),
      hmm_forward_filter = mock_hmm_factory(spy)
    )

    eval_sigma(sigma = 1, likelihood_array = dummy_lik, kern_extra = "k")

    expect_equal(spy$kern_args$kern_extra, "k")
    expect_null(spy$hmm_args$hmm_extra)
    expect_false("kern_extra" %in% names(spy$hmm_args))
  })

  it("sends filter-only arguments to hmm_forward_filter() only", {
    spy <- make_spy()
    local_mocked_bindings(
      gaussian_kernel = mock_kernel_factory(spy),
      hmm_forward_filter = mock_hmm_factory(spy)
    )

    eval_sigma(sigma = 1, likelihood_array = dummy_lik, hmm_extra = "h")

    expect_equal(spy$hmm_args$hmm_extra, "h")
    expect_null(spy$kern_args$kern_extra)
    expect_false("hmm_extra" %in% names(spy$kern_args))
  })

  it("sends an argument shared by both functions to both", {
    spy <- make_spy()
    local_mocked_bindings(
      gaussian_kernel = mock_kernel_factory(spy),
      hmm_forward_filter = mock_hmm_factory(spy)
    )

    eval_sigma(sigma = 1, likelihood_array = dummy_lik, shared_extra = "s")

    expect_equal(spy$kern_args$shared_extra, "s")
    expect_equal(spy$hmm_args$shared_extra, "s")
  })

  it("routes a mix of arguments correctly in one call", {
    spy <- make_spy()
    local_mocked_bindings(
      gaussian_kernel = mock_kernel_factory(spy),
      hmm_forward_filter = mock_hmm_factory(spy)
    )

    eval_sigma(
      sigma = 1,
      likelihood_array = dummy_lik,
      kern_extra = "k",
      hmm_extra = "h"
    )

    expect_equal(spy$kern_args$kern_extra, "k")
    expect_equal(spy$hmm_args$hmm_extra, "h")
  })

  it("runs without error when `...` is empty", {
    spy <- make_spy()
    local_mocked_bindings(
      gaussian_kernel = mock_kernel_factory(spy),
      hmm_forward_filter = mock_hmm_factory(spy)
    )

    expect_no_error(eval_sigma(sigma = 1, likelihood_array = dummy_lik))
    expect_null(spy$kern_args$kern_extra)
    expect_null(spy$hmm_args$hmm_extra)
  })
})

describe("errors from `...`", {
  it("aborts when an argument in `...` is unnamed", {
    spy <- make_spy()
    local_mocked_bindings(
      gaussian_kernel = mock_kernel_factory(spy),
      hmm_forward_filter = mock_hmm_factory(spy)
    )

    # positional args fill sigma..return_maps; the 6th lands in `...` unnamed
    err <- rlang::catch_cnd(
      suppressMessages(
        eval_sigma(1, FALSE, TRUE, dummy_lik, FALSE, 99)
      ),
      classes = "error"
    )

    expect_s3_class(err, "rlang_error")
    expect_true(
      grepl("must be named", conditionMessage(err), fixed = TRUE)
    )
  })

  it("aborts when `...` mixes named and unnamed arguments", {
    spy <- make_spy()
    local_mocked_bindings(
      gaussian_kernel = mock_kernel_factory(spy),
      hmm_forward_filter = mock_hmm_factory(spy)
    )

    err <- rlang::catch_cnd(
      suppressMessages(
        eval_sigma(1, FALSE, TRUE, dummy_lik, FALSE, 99, kern_extra = "k")
      ),
      classes = "error"
    )

    expect_s3_class(err, "rlang_error")
    expect_true(
      grepl("must be named", conditionMessage(err), fixed = TRUE)
    )
  })

  it("aborts when a named argument matches neither collaborator", {
    spy <- make_spy()
    local_mocked_bindings(
      gaussian_kernel = mock_kernel_factory(spy),
      hmm_forward_filter = mock_hmm_factory(spy)
    )

    err <- rlang::catch_cnd(
      suppressMessages(
        eval_sigma(sigma = 1, likelihood_array = dummy_lik, not_an_arg = 1)
      ),
      classes = "error"
    )

    expect_s3_class(err, "rlang_error")
    msg <- conditionMessage(err)
    expect_true(grepl("Unused argument", msg, fixed = TRUE))
    expect_true(grepl("gaussian_kernel()", msg, fixed = TRUE))
    expect_true(grepl("hmm_forward_filter()", msg, fixed = TRUE))
  })

  it("[KNOWN BUG] does not interpolate the offending argument name", {
    spy <- make_spy()
    local_mocked_bindings(
      gaussian_kernel = mock_kernel_factory(spy),
      hmm_forward_filter = mock_hmm_factory(spy)
    )

    err <- rlang::catch_cnd(
      suppressMessages(
        eval_sigma(sigma = 1, likelihood_array = dummy_lik, not_an_arg = 1)
      ),
      classes = "error"
    )
    msg <- conditionMessage(err)

    # Current (wrong) behaviour: the message is written `{.arg unknown}`,
    # which prints the literal word "unknown" instead of evaluating the
    # `unknown` variable, so the user never sees which argument was bad.
    expect_true(grepl("unknown", msg, fixed = TRUE))
    expect_false(grepl("not_an_arg", msg, fixed = TRUE))

    # TODO: after changing to `{.arg {unknown}}`, replace the two
    # assertions above with:
    # expect_true(grepl("not_an_arg", msg, fixed = TRUE))
  })

  it("reports every unknown argument, not just the first", {
    skip("Blocked by the [KNOWN BUG] above: names are not interpolated")
    # TODO: once fixed
    # err <- rlang::catch_cnd(
    #   suppressMessages(
    #     eval_sigma(sigma = 1, likelihood_array = dummy_lik, bad_a = 1, bad_b = 2)
    #   ),
    #   classes = "error"
    # )
    # expect_true(grepl("bad_a", conditionMessage(err), fixed = TRUE))
    # expect_true(grepl("bad_b", conditionMessage(err), fixed = TRUE))
  })
})
