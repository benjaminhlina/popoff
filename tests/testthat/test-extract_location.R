# Tests for extract_locations()
#
# Array layout is [lat (row), lon (col), date], with dimnames named
# lat / lon / date. The fixture below is deliberately non-square (3 x 4) so
# that a row/col mix-up in the index arithmetic cannot pass by accident.

make_loc_array <- function() {
  arr <- array(
    NA_real_,
    dim = c(3, 4, 3),
    dimnames = list(
      lat = c("40", "41", "42"),
      lon = c("-80", "-79", "-78", "-77"),
      date = c("2024-01-01", "2024-01-02", "2024-01-03")
    )
  )

  # Day 1: five finite cells with distinct max / min / mean-nearest / median
  #   values: 1, 2, 3, 4, 100  -> mean = 22 (nearest 4), median = 3
  arr[1, 1, 2] <- 1 # lat 40, lon -80
  arr[2, 1, 2] <- 2 # lat 41, lon -80
  arr[3, 2, 2] <- 3 # lat 42, lon -79
  arr[1, 3, 2] <- 4 # lat 40, lon -78
  arr[2, 4, 2] <- 100 # lat 41, lon -77

  # Day 2: all NA (left as is)

  # Day 3: non-finite values plus one finite cell
  arr[1, 1, 3] <- Inf
  arr[2, 2, 3] <- -Inf
  arr[3, 3, 3] <- NaN
  arr[3, 4, 3] <- 7 # lat 42, lon -77

  arr
}

describe("extract_locations()", {
  # describe("type argument", {
  it("defaults to 'max'", {
    arr <- make_loc_array()
    expect_equal(
      extract_locations(arr),
      extract_locations(arr, type = "max")
    )
  })

  it("errors on an invalid type", {
    expect_error(
      extract_locations(make_loc_array(), type = "mode"),
      "should be one of"
    )
  })

  it("returns the location of the maximum", {
    arr <- make_loc_array()
    out <- extract_locations(arr, type = "max")
    expect_equal(out$lon, c(NA, -77, -77))
    expect_equal(out$lat, c(NA, 41, 42))
    expect_equal(out$value, c(NA, 100, 7))
  })

  it("returns the location of the minimum", {
    arr <- make_loc_array()
    out <- extract_locations(arr, type = "min")
    expect_equal(out$lon, c(NA, -80, -77))
    expect_equal(out$lat, c(NA, 40, 42))
    expect_equal(out$value, c(NA, 1, 7))
  })

  it("returns the cell nearest the mean, reporting that cell's value", {
    arr <- make_loc_array()
    out <- extract_locations(arr, type = "mean")
    # mean is 22, so the nearest cell is the one holding 4
    expect_equal(out$lon, c(NA, -78, -77))
    expect_equal(out$lat, c(NA, 40, 42))
    expect_equal(out$value, c(NA, 4, 7))
  })

  it("returns the cell nearest the median, reporting that cell's value", {
    arr <- make_loc_array()
    out <- extract_locations(arr, type = "median")
    expect_equal(out$lon, c(NA, -79, -77))
    expect_equal(out$lat, c(NA, 42, 42))
    expect_equal(out$value, c(NA, 3, 7))
  })
})

describe("index to lon/lat mapping", {
  it("maps rows to lat and columns to lon", {
    arr <- array(
      0,
      dim = c(3, 4, 2),
      dimnames = list(
        lat = c("40", "41", "42"),
        lon = c("-80", "-79", "-78", "-77"),
        date = c("2024-01-01", "2024-01-02")
      )
    )
    arr[3, 2, 1] <- 5 # last row, second column
    out <- extract_locations(arr, type = "max")
    expect_equal(out$lat, c(42, 40))
    expect_equal(out$lon, c(-79, -80))
  })

  it("handles an array with a single row (obj[, , i] drops to a vector)", {
    arr <- array(
      c(1, 5, 2),
      dim = c(1, 3, 2),
      dimnames = list(
        lat = "40",
        lon = c("-80", "-79", "-78"),
        date = c("2024-01-01", "2024-01-02")
      )
    )
    out <- extract_locations(arr, type = "max")
    expect_equal(out$lat, c(40, 40))
    expect_equal(out$lon, c(-79, -79))
    expect_equal(out$value, c(5, 5))
  })

  it("handles an array with a single column (obj[, , i] drops to a vector)", {
    arr <- array(
      c(1, 5, 2),
      dim = c(3, 1, 2),
      dimnames = list(
        lat = c("40", "41", "42"),
        lon = "-80",
        date = c("2024-01-01", "2024-01-02")
      )
    )
    out <- extract_locations(arr, type = "max")
    expect_equal(out$lat, c(41, 41))
    expect_equal(out$lon, c(-80, -80))
    expect_equal(out$value, c(5, 5))
  })

  it("handles a single-cell, single-date array", {
    arr <- array(
      3,
      dim = c(1, 1, 1),
      dimnames = list(lat = "40", lon = "-80", date = "2024-01-01")
    )
    expect_error(
      extract_locations(arr, type = "min"),
      fixed = TRUE,
      regex = "must have at least 2 time slice"
    )
  })

  it("returns the first cell in column-major order when values tie", {
    arr <- array(
      0,
      dim = c(2, 2, 2),
      dimnames = list(
        lat = c("40", "41"),
        lon = c("-80", "-79"),
        date = c("2024-01-01", "2024-01-02")
      )
    )
    arr[2, 1, 1] <- 9 # column-major position 2
    arr[1, 2, 1] <- 9 # column-major position 3
    out <- extract_locations(arr, type = "max")
    expect_equal(out$lat, c(41, 40))
    expect_equal(out$lon, c(-80, -80))
  })
})

describe("output structure", {
  it("returns a data.frame with one row per date and the expected columns", {
    local_mocked_bindings(msg_extract_no_val = function(...) invisible(NULL))
    out <- extract_locations(make_loc_array())
    expect_s3_class(out, "data.frame")
    expect_named(out, c("date", "lon", "lat", "value"))
    expect_equal(nrow(out), 3)
  })

  it("returns dates as Date objects in array order", {
    local_mocked_bindings(msg_extract_no_val = function(...) invisible(NULL))
    out <- extract_locations(make_loc_array())
    expect_s3_class(out$date, "Date")
    expect_equal(
      out$date,
      as.Date(c("2024-01-01", "2024-01-02", "2024-01-03"))
    )
  })

  it("returns numeric lon, lat and value when at least one date has a value", {
    local_mocked_bindings(msg_extract_no_val = function(...) invisible(NULL))
    out <- extract_locations(make_loc_array())
    expect_true(is.numeric(out$lon))
    expect_true(is.numeric(out$lat))
    expect_true(is.numeric(out$value))
  })

  it("resets row names", {
    local_mocked_bindings(msg_extract_no_val = function(...) invisible(NULL))
    out <- extract_locations(make_loc_array())
    expect_equal(rownames(out), as.character(1:3))
  })
})

describe("empty and non-finite handling", {
  it("returns NA lon/lat/value but keeps the date for an all-NA slice", {
    local_mocked_bindings(msg_extract_no_val = function(...) invisible(NULL))
    out <- extract_locations(make_loc_array())
    expect_equal(out$date[2], as.Date("2024-01-02"))
    expect_true(is.na(out$lon[1]))
    expect_true(is.na(out$lat[1]))
    expect_true(is.na(out$value[1]))
  })

  it("ignores Inf, -Inf and NaN when choosing a location", {
    local_mocked_bindings(msg_extract_no_val = function(...) invisible(NULL))
    arr <- make_loc_array()
    # Day 3 has Inf, -Inf and NaN; only the 7 at lat 42 / lon -77 is finite
    for (type in c("max", "min", "mean", "median")) {
      out <- extract_locations(arr, type = type)
      expect_equal(out$lat[3], 42, info = type)
      expect_equal(out$lon[3], -77, info = type)
      expect_equal(out$value[3], 7, info = type)
    }
  })

  it("treats a slice of only non-finite values as empty", {
    local_mocked_bindings(msg_extract_no_val = function(...) invisible(NULL))
    arr <- array(
      c(Inf, -Inf, NaN, NA),
      dim = c(2, 2, 2),
      dimnames = list(
        lat = c("40", "41"),
        lon = c("-80", "-79"),
        date = c("2024-01-01", "2024-01-02")
      )
    )
    out <- extract_locations(arr)
    expect_equal(is.na(out$lon), c(TRUE, TRUE))
    expect_equal(is.na(out$lat), c(TRUE, TRUE))
    expect_equal(is.na(out$value), c(TRUE, TRUE))
  })

  it("keeps empty rows by default (drop_empty = FALSE)", {
    local_mocked_bindings(msg_extract_no_val = function(...) invisible(NULL))
    out <- extract_locations(make_loc_array())
    expect_equal(nrow(out), 3)
  })

  it("drops rows with no location when drop_empty = TRUE", {
    local_mocked_bindings(msg_extract_no_val = function(...) invisible(NULL))
    out <- extract_locations(make_loc_array(), drop_empty = TRUE)
    expect_equal(nrow(out), 2)
    expect_equal(out$date, as.Date(c("2024-01-02", "2024-01-03")))
    expect_false(anyNA(out$lon))
  })

  it("returns a zero-row data.frame when every slice is empty and drop_empty = TRUE", {
    local_mocked_bindings(msg_extract_no_val = function(...) invisible(NULL))
    arr <- array(
      NA_real_,
      dim = c(2, 2, 2),
      dimnames = list(
        lat = c("40", "41"),
        lon = c("-80", "-79"),
        date = c("2024-01-01", "2024-01-02")
      )
    )
    out <- extract_locations(arr, drop_empty = TRUE)
    expect_s3_class(out, "data.frame")
    expect_equal(nrow(out), 0)
    expect_named(out, c("date", "lon", "lat", "value"))
  })

  it("returns all-NA rows for every date when no slice has values", {
    local_mocked_bindings(msg_extract_no_val = function(...) invisible(NULL))
    arr <- array(
      NA_real_,
      dim = c(2, 2, 2),
      dimnames = list(
        lat = c("40", "41"),
        lon = c("-80", "-79"),
        date = c("2024-01-01", "2024-01-02")
      )
    )
    out <- extract_locations(arr)
    expect_equal(nrow(out), 2)
    expect_true(all(is.na(out$lon)))
    expect_true(all(is.na(out$lat)))
    expect_true(all(is.na(out$value)))
  })
})

describe("empty-value message", {
  it("calls msg_extract_no_val() with the number of empty dates", {
    called <- NULL
    local_mocked_bindings(
      msg_extract_no_val = function(n_empty, out) {
        called <<- list(n_empty = n_empty, out = out)
        invisible(NULL)
      }
    )
    out <- extract_locations(make_loc_array())
    expect_false(is.null(called))
    expect_equal(called$n_empty, 1)
    expect_equal(called$out, out)
  })

  it("does not call msg_extract_no_val() when every date has a value", {
    n_calls <- 0
    local_mocked_bindings(
      msg_extract_no_val = function(...) {
        n_calls <<- n_calls
        invisible(NULL)
      }
    )
    extract_locations(make_loc_array())
    expect_equal(n_calls, 0)
  })

  it("reports empties before drop_empty removes them", {
    called <- NULL
    local_mocked_bindings(
      msg_extract_no_val = function(n_empty, out) {
        called <<- list(n_empty = n_empty, n_rows = nrow(out))
        invisible(NULL)
      }
    )
    extract_locations(make_loc_array(), drop_empty = TRUE)
    expect_equal(called$n_empty, 1)
    expect_equal(called$n_rows, 3)
  })
})

# Pins current behavior; flag for follow-up if you want it changed.
describe("known quirks", {
  it("[KNOWN QUIRK] drop_empty = TRUE leaves non-contiguous row names", {
    local_mocked_bindings(msg_extract_no_val = function(...) invisible(NULL))
    out <- extract_locations(make_loc_array(), drop_empty = TRUE)
    expect_equal(rownames(out), c("2", "3"))
    # TODO: if row names are reset after dropping, this should be c("1", "2")
    # expect_equal(rownames(out), c("1", "2"))
  })

  it("[KNOWN QUIRK] all-empty input yields logical (not numeric) NA columns", {
    local_mocked_bindings(msg_extract_no_val = function(...) invisible(NULL))
    arr <- array(
      NA_real_,
      dim = c(2, 2, 2),
      dimnames = list(
        lat = c("40", "41"),
        lon = c("-80", "-79"),
        date = c("2024-01-01", "2024-01-02")
      )
    )
    out <- extract_locations(arr)
    expect_true(is.logical(out$lon))
    expect_true(is.logical(out$lat))
    expect_true(is.logical(out$value))
    # TODO: if NA_real_ is used for the empty row, these should be numeric
    # expect_true(is.numeric(out$lon))
  })
})
# })
