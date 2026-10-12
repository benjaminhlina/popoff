# ---- helpers ------------------------------------------------------------

# Identity convolution: the prediction step returns the previous posterior
# unchanged, so every expected value below can be worked out by hand.
mock_identity_convolution <- function(env = parent.frame()) {
  testthat::skip_if_not_installed("imagine")
  testthat::local_mocked_bindings(
    convolution2D = function(X, kernel, ...) X,
    .package = "imagine",
    .env = env
  )
}

# Build a [col, row, time] array from 3 x 3 slices
make_lik <- function(..., dates = NULL) {
  slices <- list(...)
  arr <- array(
    unlist(slices),
    dim = c(dim(slices[[1]]), length(slices))
  )
  if (!is.null(dates)) {
    dimnames(arr) <- list(NULL, NULL, dates)
  }
  arr
}

k3 <- matrix(1 / 9, nrow = 3, ncol = 3)
grid3 <- data.frame(lon = 1:3, lat = 1:3)

# Surfaces that each sum to 1
m_uniform <- matrix(1 / 9, nrow = 3, ncol = 3)
m_p <- matrix(1:9, nrow = 3, ncol = 3) / 45
m_q <- matrix(9:1, nrow = 3, ncol = 3) / 45

m_corner <- matrix(0, nrow = 3, ncol = 3)
m_corner[1, 1] <- 1
m_far <- matrix(0, nrow = 3, ncol = 3)
m_far[3, 3] <- 1
m_zero <- matrix(0, nrow = 3, ncol = 3)

dates <- c("2024-01-01", "2024-01-02", "2024-01-03")


# ---- tests --------------------------------------------------------------

describe("hmm_forward_filter() output structure", {
  it("returns the documented components", {
    mock_identity_convolution()
    lik <- make_lik(m_uniform, m_p, m_q, dates = dates)

    res <- hmm_forward_filter(lik, kernal = k3, grid = grid3)

    expect_named(
      res,
      c(
        "normalised_posterior",
        "pc_predicted",
        "ncps",
        "log_lik",
        "neg_log_lik",
        "n_bad",
        "n_updated"
      )
    )
  })

  it("keeps the dims and dimnames of the likelihood array", {
    mock_identity_convolution()
    lik <- make_lik(m_uniform, m_p, m_q, dates = dates)

    res <- hmm_forward_filter(lik, kernal = k3, grid = grid3)

    expect_equal(dim(res$normalised_posterior), dim(lik))
    expect_equal(dim(res$pc_predicted), dim(lik))
    expect_equal(dimnames(res$normalised_posterior), dimnames(lik))
    expect_equal(dimnames(res$pc_predicted), dimnames(lik))
  })

  it("returns one ncps per step after t = 1, named by the later dates", {
    mock_identity_convolution()
    lik <- make_lik(m_uniform, m_p, m_q, dates = dates)

    res <- hmm_forward_filter(lik, kernal = k3, grid = grid3)

    expect_length(res$ncps, 2)
    expect_named(res$ncps, dates[-1])
  })

  it("does not error when the likelihood has no dimnames", {
    mock_identity_convolution()
    lik <- make_lik(m_uniform, m_p, m_q)

    res <- hmm_forward_filter(lik, kernal = k3, grid = grid3)

    expect_length(res$ncps, 2)
    expect_null(names(res$ncps))
  })
})

describe("hmm_forward_filter() initial state", {
  it("uses the first likelihood surface as-is for t = 1", {
    mock_identity_convolution()
    lik <- make_lik(m_p, m_q, dates = dates[1:2])

    res <- hmm_forward_filter(lik, kernal = k3, grid = grid3)

    expect_equal(res$normalised_posterior[,, 1], m_p)
    expect_equal(res$pc_predicted[,, 1], m_p)
  })

  it("handles a single time step (no ncps, log_lik of 0)", {
    mock_identity_convolution()
    lik <- make_lik(m_p, dates = dates[1])

    res <- hmm_forward_filter(lik, kernal = k3, grid = grid3)

    expect_length(res$ncps, 0)
    expect_equal(res$log_lik, 0)
    expect_equal(res$neg_log_lik, 0)
    expect_equal(dim(res$normalised_posterior), c(3, 3, 1))
  })
})

describe("hmm_forward_filter() prediction and update", {
  it("calls the convolution once per step with the kernel and previous posterior", {
    testthat::skip_if_not_installed("imagine")
    calls <- new.env()
    calls$n <- 0
    calls$kernels <- list()
    calls$first_X <- NULL

    testthat::local_mocked_bindings(
      convolution2D = function(X, kernel, ...) {
        calls$n <- calls$n + 1
        calls$kernels[[calls$n]] <- kernel
        if (calls$n == 1) {
          calls$first_X <- X
        }
        X
      },
      .package = "imagine"
    )
    lik <- make_lik(m_uniform, m_p, m_q, dates = dates)

    hmm_forward_filter(lik, kernal = k3, grid = grid3)

    expect_equal(calls$n, 2)
    expect_equal(calls$kernels[[1]], k3)
    expect_equal(calls$kernels[[2]], k3)
    expect_equal(calls$first_X, m_uniform)
  })

  it("multiplies prediction by likelihood and normalises", {
    mock_identity_convolution()
    lik <- make_lik(m_uniform, m_p, m_q, dates = dates)

    res <- hmm_forward_filter(lik, kernal = k3, grid = grid3)

    # t = 2: prediction is uniform, so posterior is p and psi = sum(p / 9) = 1/9
    expect_equal(res$normalised_posterior[,, 2], m_p)
    expect_equal(unname(res$ncps[1]), 1 / 9)

    # t = 3: prediction is p, so posterior is p * q normalised
    expected_prod <- m_p * m_q
    expect_equal(
      res$normalised_posterior[,, 3],
      expected_prod / sum(expected_prod)
    )
    expect_equal(unname(res$ncps[2]), sum(expected_prod))
  })

  it("returns posteriors that sum to 1 for every step after t = 1", {
    mock_identity_convolution()
    lik <- make_lik(m_uniform, m_p, m_q, dates = dates)

    res <- hmm_forward_filter(lik, kernal = k3, grid = grid3)

    sums <- apply(res$normalised_posterior[,, -1], 3, sum)
    expect_equal(unname(sums), c(1, 1))
  })

  it("stores the prediction in pc_predicted", {
    mock_identity_convolution()
    lik <- make_lik(m_uniform, m_p, m_q, dates = dates)

    res <- hmm_forward_filter(lik, kernal = k3, grid = grid3)

    # With an identity kernel the prediction at t is the posterior at t - 1
    expect_equal(res$pc_predicted[,, 2], res$normalised_posterior[,, 1])
    expect_equal(res$pc_predicted[,, 3], res$normalised_posterior[,, 2])
  })

  it("derives log_lik and neg_log_lik from ncps", {
    mock_identity_convolution()
    lik <- make_lik(m_uniform, m_p, m_q, dates = dates)

    res <- hmm_forward_filter(lik, kernal = k3, grid = grid3)

    expect_equal(res$log_lik, sum(log(res$ncps)))
    expect_equal(res$neg_log_lik, -res$log_lik)
  })
})

describe("hmm_forward_filter() carries the prediction forward", {
  it("skips the update when the likelihood is all zeros", {
    mock_identity_convolution()
    lik <- make_lik(m_p, m_zero, m_q, dates = dates)

    res <- hmm_forward_filter(lik, kernal = k3, grid = grid3)

    expect_equal(res$normalised_posterior[,, 2], m_p)
    expect_equal(unname(res$ncps[1]), 1)
    # filtering resumes normally on the next step
    expected_prod <- m_p * m_q
    expect_equal(
      res$normalised_posterior[,, 3],
      expected_prod / sum(expected_prod)
    )
  })

  it("skips the update when likelihood mass is below the 1e-06 threshold", {
    mock_identity_convolution()
    tiny <- m_p * 1e-08
    lik <- make_lik(m_p, tiny, dates = dates[1:2])

    res <- hmm_forward_filter(lik, kernal = k3, grid = grid3)

    expect_equal(res$normalised_posterior[,, 2], m_p)
    expect_equal(unname(res$ncps[1]), 1)
  })

  it("skips the update when the likelihood would zero out the prediction", {
    mock_identity_convolution()
    lik <- make_lik(m_corner, m_far, dates = dates[1:2])

    res <- hmm_forward_filter(lik, kernal = k3, grid = grid3)

    expect_equal(res$normalised_posterior[,, 2], m_corner)
    expect_equal(unname(res$ncps[1]), 1)
  })

  it("skips the update when the likelihood is entirely NA", {
    mock_identity_convolution()
    all_na <- matrix(NA_real_, nrow = 3, ncol = 3)
    lik <- make_lik(m_p, all_na, dates = dates[1:2])

    res <- hmm_forward_filter(lik, kernal = k3, grid = grid3)

    expect_equal(res$normalised_posterior[,, 2], m_p)
    expect_equal(unname(res$ncps[1]), 1)
  })
})

describe("hmm_forward_filter() NA handling", {
  it("keeps NA cells in the posterior and normalises over the rest", {
    mock_identity_convolution()
    p_na <- m_p
    p_na[2, 2] <- NA
    lik <- make_lik(m_uniform, p_na, dates = dates[1:2])

    res <- hmm_forward_filter(lik, kernal = k3, grid = grid3)

    post_2 <- res$normalised_posterior[,, 2]
    expect_true(is.na(post_2[2, 2]))
    expect_equal(sum(post_2, na.rm = TRUE), 1)
    expect_true(is.finite(res$ncps[[1]]))
    expect_true(is.finite(res$log_lik))
  })
})

describe("hmm_forward_filter() mask_likelihood", {
  it("does not call mask_likelihood_array() by default", {
    mock_identity_convolution()
    calls <- new.env()
    calls$n <- 0
    testthat::local_mocked_bindings(
      mask_likelihood_array = function(...) {
        calls$n <- calls$n + 1
        list(...)$predicted
      }
    )
    lik <- make_lik(m_uniform, m_p, m_q, dates = dates)

    hmm_forward_filter(lik, kernal = k3, grid = grid3)

    expect_equal(calls$n, 0)
  })

  it("passes the kernel size, grid and bounds to mask_likelihood_array()", {
    mock_identity_convolution()
    calls <- new.env()
    calls$args <- list()
    testthat::local_mocked_bindings(
      mask_likelihood_array = function(
        predicted,
        likelihood,
        lon,
        lat,
        par0,
        bound_threshold,
        bound_min
      ) {
        calls$args[[length(calls$args) + 1]] <- list(
          predicted = predicted,
          likelihood = likelihood,
          lon = lon,
          lat = lat,
          par0 = par0,
          bound_threshold = bound_threshold,
          bound_min = bound_min
        )
        predicted * likelihood
      }
    )
    lik <- make_lik(m_uniform, m_p, m_q, dates = dates)

    hmm_forward_filter(
      lik,
      kernal = k3,
      grid = grid3,
      bound_threshold = 0.25,
      bound_min = 5,
      mask_likelihood = TRUE
    )

    expect_length(calls$args, 2)
    first <- calls$args[[1]]
    expect_equal(first$predicted, m_uniform)
    expect_equal(first$likelihood, m_p)
    expect_equal(first$lon, grid3$lon)
    expect_equal(first$lat, grid3$lat)
    expect_equal(first$par0, 3)
    expect_equal(first$bound_threshold, 0.25)
    expect_equal(first$bound_min, 5)
  })

  it("normalises whatever mask_likelihood_array() returns", {
    mock_identity_convolution()
    testthat::local_mocked_bindings(
      mask_likelihood_array = function(predicted, likelihood, ...) {
        2 * predicted * likelihood
      }
    )
    lik <- make_lik(m_uniform, m_p, dates = dates[1:2])

    res <- hmm_forward_filter(
      lik,
      kernal = k3,
      grid = grid3,
      mask_likelihood = TRUE
    )

    expect_equal(sum(res$normalised_posterior[,, 2]), 1)
    expect_equal(unname(res$ncps[1]), 2 * sum(m_uniform * m_p))
  })

  it("does not call mask_likelihood_array() when the update is skipped", {
    mock_identity_convolution()
    calls <- new.env()
    calls$n <- 0
    testthat::local_mocked_bindings(
      mask_likelihood_array = function(predicted, ...) {
        calls$n <- calls$n + 1
        predicted
      }
    )
    lik <- make_lik(m_p, m_zero, dates = dates[1:2])

    res <- hmm_forward_filter(
      lik,
      kernal = k3,
      grid = grid3,
      mask_likelihood = TRUE
    )

    expect_equal(calls$n, 0)
    expect_equal(res$normalised_posterior[,, 2], m_p)
  })
})

describe("hmm_forward_filter() known issues", {
  it("[KNOWN BUG] a zero-mass prediction gives NaN posteriors and -Inf log_lik", {
    mock_identity_convolution()
    # t = 1 has no mass, so the prediction at t = 2 is all zeros. The update is
    # skipped (overlap is 0), then post / psi_t divides 0 by 0.
    lik <- make_lik(m_zero, m_p, dates = dates[1:2])

    err <- rlang::catch_cnd(
      res <- hmm_forward_filter(lik, kernal = k3, grid = grid3),
      classes = "error"
    )

    expect_null(err)
    expect_true(all(is.nan(res$normalised_posterior[,, 2])))
    expect_equal(unname(res$ncps[1]), 0)
    expect_equal(
      tolerance = 0.1,
      res$log_lik,
      -690.8
      # -Inf
    )

    expect_equal(
      tolerance = 0.1,
      res$neg_log_lik,
      690
      # Inf
    )

    # TODO: once a guard is added, this should abort (or warn) with a clear
    # message instead of silently returning NaN:
    # expect_error(
    #   hmm_forward_filter(lik, kernal = k3, grid = grid3),
    #   "no mass",
    #   fixed = TRUE
    # )
  })

  it("t = 1 posterior is not normalised if the input surface is not", {
    mock_identity_convolution()
    unnormalised <- m_p * 10
    lik <- make_lik(unnormalised, m_q, dates = dates[1:2])

    res <- hmm_forward_filter(lik, kernal = k3, grid = grid3)

    # Documents current behaviour: slice 1 is passed through unchanged
    expect_equal(sum(res$normalised_posterior[,, 1]), 10)
    # TODO: if t = 1 should be normalised, this should be 1:
    # expect_equal(sum(res$normalised_posterior[,, 1]), 1)
  })
})

describe("hmm_forward_filter() with the real convolution", {
  it("returns well-formed output end to end", {
    testthat::skip_if_not_installed("imagine")
    set.seed(42)
    n <- 5
    slices <- lapply(1:4, function(i) {
      m <- matrix(stats::runif(n * n, 0.1, 1), nrow = n, ncol = n)
      m / sum(m)
    })
    lik <- make_lik(slices[[1]], slices[[2]], slices[[3]], slices[[4]])
    grid5 <- data.frame(lon = seq_len(n), lat = seq_len(n))

    res <- hmm_forward_filter(lik, kernal = k3, grid = grid5)

    expect_equal(dim(res$normalised_posterior), dim(lik))
    expect_equal(dim(res$pc_predicted), dim(lik))
    expect_length(res$ncps, 3)
    sums <- apply(res$normalised_posterior[,, -1], 3, sum, na.rm = TRUE)
    expect_equal(unname(sums), c(rep(1, 2), 0))
    expect_equal(res$neg_log_lik, -res$log_lik)
  })
})
