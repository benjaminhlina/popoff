# ----- load packages ----

{
  library(ggplot2)
  library(popr)
  library(terra)
  library(tidyterra)
  library(scales)
}

# ----- load bathymetry for area ----

bathy <- load_bathy_raster(bathy_matrix, crs = "EPSG:6500")

bathy
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

# ----- prep bathy for likelihood calcualtion -----
preped_bathy <- prep_bathy_likdepth(bathy = bathy, focal_dim = 3)
preped_bathy
# ------ prep log for summary prior to likelihood calculation ------
preped_log <- prep_summary_likdepth(raw_log = pdt, iniloc = iniloc, bin = 3600)

preped_log

# ---- calculaite likeilhoood depth -----
likedepth_transmit <- calc_likdepth_transmit(
  preped_log = preped_log,
  preped_bathy = preped_bathy,
  daily_sd_depth = 0.57
)

# ----- plot the first couple ----
ldt_dims <- dim(likedepth_transmit)

days <- seq(1, 2, by = 1)

bathy_matrix |>
  str()
likedepth_transmit[,, .1] |>
  str()
# this just confirms that it its working how it is intended will be making tests

test <- days |>
  purrr::map(
    ~ load_bathy_raster(likedepth_transmit[,, .x], crs = "EPSG:6500"),
    # t(),
    .progress = TRUE
  )

# m <- likedepth_transmit[,, 1]

# test <- load_bathy_raster(m, crs = "EPSG:6500")

# test <- terra::t(test)
test |>
  purrr::map(
    ~ ggplot() +
      geom_spatraster(data = .x) +
      coord_sf() +
      scale_fill_viridis_c(option = "H")
  )
