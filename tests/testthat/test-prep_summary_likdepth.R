describe("prep_summary_likdepth()", {
  # ---- Setup Mock Datasets ----
  test_raw_log <- data.frame(
    fish_id = "Fish_101",
    date = rep(as.Date(c("2026-01-01", "2026-01-02", "2026-01-03")), each = 4),
    date_time = seq(
      from = as.POSIXct("2026-01-01 00:00:00", tz = "UTC"),
      by = "6 hours",
      length.out = 12
    ),
    depth = c(10, 20, 15, 30, 25, 35, 40, 20, 15, 10, 5, 2),
    temperature = c(12, 11.5, 11, 10.5, 10, 9.5, 9, 10, 11, 12, 12.5, 13)
  )

  test_iniloc <- data.frame(
    fish_id = c("Fish_101", "Fish_101", "Fish_999"),
    event = c("tag", "pop", "tag"),
    date = as.POSIXct(c(
      "2026-01-01 00:12:00",
      "2026-01-02 01:12:01",
      "2025-05-01 00:00:00"
    ))
  )

  # ---- 1. Input Validation Checks ----
  it("throws error when required arguments are missing", {
    expect_error(prep_summary_likdepth(iniloc = test_iniloc, bin = 3600))
    expect_error(prep_summary_likdepth(raw_log = test_raw_log, bin = 3600))
    expect_error(prep_summary_likdepth(
      raw_log = test_raw_log,
      iniloc = test_iniloc
    ))
  })

  # ---- 2. Core Execution & Structure Output ----
  it("processes data correctly and returns expected list components", {
    res <- prep_summary_likdepth(
      raw_log = test_raw_log,
      iniloc = test_iniloc,
      bin = 3600,
      make_plot = FALSE
    )
    expect

    expect_type(res, "list")
    expect_named(
      res,
      c("raw_log_filter", "summary_table", "depth_profile", "time_res", "bin")
    )

    # Output structure assertions
    expect_s3_class(res$raw_log_filter, "data.frame")
    expect_s3_class(res$summary_table, "data.frame")
    expect_null(res$depth_profile)
    expect_equal(res$bin, 3600)
    expect_equal(res$time_res, 21600) # 6 hours = 21600 seconds
  })
  it("check msg", {
    expect_message(
      prep_summary_likdepth(
        raw_log = test_raw_log,
        iniloc = test_iniloc,
        bin = 360000,
        make_plot = FALSE
      ),
      regexp = "Time window, which is the `bin` (e.g., 3600 sec = 1 h) divided by the PSAT tranmission time (e.g., 360), is",
      fixed = TRUE
    )
  })
  # ---- 3. Filtering & Depth Inversion Logic ----
  it("filters date range correctly and inverts depth values", {
    res <- prep_summary_likdepth(
      raw_log = test_raw_log,
      iniloc = test_iniloc,
      bin = 3600,
      make_plot = FALSE
    )

    # Date range filtered strictly between tag (2026-01-01) and pop (2026-01-02)
    expect_equal(
      unique(res$raw_log_filter$date),
      as.Date(c("2026-01-01", "2026-01-02"))
    )
    expect_true(all(res$raw_log_filter$date <= as.Date("2026-01-02")))

    # Inverted depth transformation check (positive becomes negative)
    expect_true(all(res$raw_log_filter$depth <= 0))
    expect_equal(res$raw_log_filter$depth[1], -10)
  })

  # ---- 4. Plot Generation & File Saving ----
  it("generates a ggplot object when make_plot is TRUE", {
    res <- prep_summary_likdepth(
      raw_log = test_raw_log,
      iniloc = test_iniloc,
      bin = 3600,
      make_plot = TRUE
    )

    expect_s3_class(res$depth_profile, "ggplot")
  })

  it("saves the plot PNG image when output_dir is provided", {
    temp_dir <- tempdir()
    expected_file <- file.path(temp_dir, "depth_profile_Fish_101.png")

    # Clean up before test
    if (file.exists(expected_file)) {
      file.remove(expected_file)
    }

    prep_summary_likdepth(
      raw_log = test_raw_log,
      iniloc = test_iniloc,
      bin = 3600,
      output_dir = temp_dir,
      make_plot = TRUE
    )

    expect_true(file.exists(expected_file))

    # Clean up after test
    if (file.exists(expected_file)) file.remove(expected_file)
  })

  # ---- 5. Safe Utility Function Coverage ----
  it("handles empty and NA vectors correctly in internal helper functions", {
    expect_equal(.safe_max(c(NA, NA)), NA_real_)
    expect_equal(.safe_max(c(1, 5, NA)), 5)

    expect_equal(.safe_range_val(c(NA, NA), type = "min"), NA_real_)
    expect_equal(.safe_range_val(c(10, 20, 5, NA), type = "min"), 5)
    expect_equal(.safe_range_val(c(10, 20, 5, NA), type = "max"), 20)
  })
})
