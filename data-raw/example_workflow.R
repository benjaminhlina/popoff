# ----- load packages ----

{
  library(ggplot2)
  library(poppr)
  library(purrr)
  library(scales)
  library(terra)
  library(tidyterra)
}

# ----- load bathymetry for the area of interest ----
# in this case Gulf of St. Lawrence

bathy <- load_bathy_raster(bathy_matrix, crs = "EPSG:6500")

# we can check that it is SpatRaster
bathy

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
    values = rescale(c(-4000, -600, -150, -60, 0, 50, 150, 400)),
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

# likedepth_transmit is a matrix with dims of bathymetry by the number
# of days observed in this case 203, each day is a new likelihood est.

# ----- we can plot the first couple of days -----

days <- seq(1, 4, by = 1)

# ---- we can convert them into spatrasters using load functions
likelihood_spr <- days |>
  map(
    ~ load_bathy_raster(likedepth_transmit[,, .x], crs = "EPSG:6500"),
    .progress = TRUE
  )

# ----- then polt them using ggplot and tidyterra -----
likelihood_spr |>
  map(
    ~ ggplot() +
      geom_spatraster(data = .x) +
      coord_sf() +
      scale_fill_viridis_c(option = "H")
  )
