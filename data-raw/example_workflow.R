# ----- load packages ----

{
  library(dplyr)
  library(ggplot2)
  library(poppr)
  library(purrr)
  library(scales)
  library(terra)
  library(tidyterra)
}

# ----- load bathymetry for the area of interest ----

bathy_matrix <- download_noaa_bathy(
  xlim = c(-72.41, -40.0),
  ylim = c(41.52, 55.0),
  res = 1
)
bathy_matrix
# in this case Gulf of St. Lawrence

# we can check that it is SpatRaster

bathy <- load_raster(bathy_matrix, crs = "EPSG:4326")

# ----- next we can plot the spat raster using tidyterra ----

ggplot() +
  geom_spatraster(data = bathy) +
  coord_sf() +
  scale_fill_gradientn(
    colours = c(
      "#031E40",
      "#0B3C73",
      "#17A2FF",
      "#7FDBFF",
      "#EAF6FF",
      "#FFE8A3",
      "#E6F38F",
      "#9AD850"
    ),
    values = scales::rescale(c(-4000, -600, -150, -60, 0, 50, 150, 400)),
    limits = c(-4000, 400),
    oob = squish,
    name = "Elevation / Depth (m)"
  )

# aesome the bathymetr looks great

# ----- next prepare bathy for likelihood calculations -----
preped_bathy <- prep_bathy_likdepth(bathy = bathy, focal_dim = 3)

# we can view the output which is a list containing 1) bathy 2) sd of bathy
# 3) idx locations 4) idx_na locations, 4) focal dimensions, 5) if weighted
preped_bathy
# ------ prep log for summary prior to likelihood calculation ------
preped_log <- prep_summary_likdepth(raw_log = pdt, iniloc = iniloc, bin = 3600)

# pepared summary is a list containging 1) raw log, 2) summary table 3) a ggplot
# of the depth and temperature profiles, 4) time resoulation, 5) bin value

preped_log

# ---- calculaite likeilhoood depth -----
likedepth_transmit <- calc_likdepth_transmit(
  preped_log = preped_log,
  preped_bathy = preped_bathy,
  daily_sd_depth = 0.57
)
dim(likedepth_transmit)
# likedepth_transmit is a matrix with dims of bathymetry by the number
# of days observed in this case 203, each day is a new likelihood est.

# ----- we can plot the first couple of days -----

days <- seq(1, 4, by = 1)

xlim <- attr(bathy_matrix, "xlim")
ylim <- attr(bathy_matrix, "ylim")


# ---- we can convert them into spatrasters using load functions
likelihood_spr <- days |>
  map(
    ~ load_raster(
      likedepth_transmit[,, .x],
      extent = c(xlim, ylim),
      crs = "EPSG:4326"
    ),
    .progress = TRUE
  )
likelihood_spr
# ----- then polt them using ggplot and tidyterra -----
likelihood_spr |>
  map(
    ~ ggplot() +
      geom_spatraster(data = .x) +
      # coord_sf() +
      scale_fill_viridis_c(option = "H")
  )


# ---- Take idx values that equla 0.5 or greater and pull the depth of those
#idx vlaues -----

?terra::extract()


bathy_vec <- bathy |>
  as.polygons(values = TRUE)
s <- sf::st_as_sf(bathy_vec)

s
ggplot() +
  geom_sf(data = s)

locs <- iniloc |>
  dplyr::select(lon, lat)
terra::extract(likelihood_spr[[1]], locs)


g <- loaded_data$grid.peter

r_extent <- terra::ext(
  min(g$xgrid) - dflon[1] / 2,
  max(g$xgrid) + dflon[length(dflon)] / 2,
  min(g$ygrid) - dflat[1] / 2,
  max(g$ygrid) + dflat[length(dflat)] / 2
)
r_extent


terra::ext(bathy) <- r_extent
terra::crs(bathy) <- "EPSG:4326"
bathy
likelihood_spr <- likelihood_spr |>
  map(
    ~ {
      terra::ext(.x) <- r_extent
      .x
    }
  )


likelihood_spr <- likelihood_spr |>
  map(
    ~ {
      terra::crs(.x) <- "EPSG:4326"
      .x
    }
  )


pts <- cbind(iniloc$lon, iniloc$lat)

cellId <- likelihood_spr |>
  map(
    ~ terra::cellFromXY(.x, pts)
  ) # linear cell index


vals <- likelihood_spr |>
  map(~ terra::extract(.x, pts)) # likelihood value(s) at those cells
vals
