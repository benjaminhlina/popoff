#' Load functions
#'
#' These functions load different data into R
#'
#' @param bathy A bathymetric feature to be brought in as a `SpatRaster` object
#' @param crs a character string that is the coordinate reference system (CRS) of the raster object e.g., `"ESPG:4326"`.
#' Default is `NULL` resulting in not supplying a CRS. See `terra::crs()` for valid CRS format.
#' @param extent default is `NULL` and if supplied `bathy` object is created by `download_noaa_bathy()`
#' then the extent is already pulled. If not then no extent will be applied. If
#' desiring an extent for the `bathy` object you need to supply a vector containing
#' the following `c(xmin, xmax, ymin, ymax)` as numerics.
#' @param ... addtional arguments to pass `terra::rast()`.
#' @details returns a `SpatRaster` object to be used with an attribute of index. To access
#' index use `attr(obj, "idx")`.
#' @seealso [terra::rast()]
#'
#' @name load_functions
#' @export

load_raster <- function(bathy, extent = NULL, crs = NULL, ...) {
  # general message
  msg_bathy(bathy = bathy)
  error_load(bathy)

  if (inherits(bathy, "SpatRaster")) {
    attr(bathy, "idx") <- .valid_idx(bathy)
    attr(bathy, "idx_na") <- .valid_idx(bathy, na = TRUE)

    if (!is.null(crs)) {
      error_crs(crs)
      terra::crs(bathy) <- crs
    }
    return(bathy)
  }

  # use terra to return a spat rast if bahty is raster or matrix
  if (inherits(bathy, "RasterLayer") || is.matrix(bathy)) {
    if (is.matrix(bathy)) {
      if (all(c("xlim", "ylim") %in% names(attributes(bathy)))) {
        extent <- c(attr(bathy, "xlim"), attr(bathy, "ylim"))
      }
      if (!is.null(extent)) {
        error_extent(extent)
        extent <- terra::ext(extent)
      }

      bathy <- terra::rast(bathy, extent = extent, ...)
    } else {
      bathy <- terra::rast(bathy, ...)
    }

    attr(bathy, "idx") <- .valid_idx(bathy)
    attr(bathy, "idx_na") <- .valid_idx(bathy, na = TRUE)

    if (!is.null(crs)) {
      error_crs(crs)
      terra::crs(bathy) <- crs
    }
    return(bathy)
  }

  if (is.character(bathy)) {
    if (grepl("\\.RData$", bathy, ignore.case = TRUE)) {
      error_path(bathy)

      e <- new.env()
      load(bathy, envir = e)
      obj_name <- ls(e)

      if (length(obj_name) == 0) {
        error_empty_file(bathy)
      }
      obj <- get(obj_name[1], envir = e)

      msg_load_rdata(obj, obj_name)
    }
    if (grepl("\\.rds$", bathy, ignore.case = TRUE)) {
      obj <- readRDS(bathy)
    }
    if (grepl("\\.tif$", bathy, ignore.case = TRUE)) {
      msg_load_raster(bathy)
      obj <- terra::rast(bathy, ...)
    }

    return(load_raster(obj))
  }
}
