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

days <- unique(preped_log$summary_table$date)

xlim <- attr(bathy_matrix, "xlim")
ylim <- attr(bathy_matrix, "ylim")


# ---- we can convert them into spatrasters using load functions
likelihood_spr <- days |>
  set_names() |>
  map(
    ~ load_raster(
      likedepth_transmit[,, as.character(.x)],
      extent = c(xlim, ylim),
      crs = "EPSG:4326"
    ),
    .progress = TRUE
  )
likelihood_spr
# ----- then polt them using ggplot and tidyterra -----
likelihood_spr |>
  imap(
    ~ ggplot() +
      geom_spatraster(data = .x) +
      coord_sf() +
      scale_fill_viridis_c(option = "H") +
      theme_bw(
        base_size = 15
      ) +
      theme(
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)
      ) +
      labs(title = .y)
  )


# ---- Take idx values that equla 0.5 or greater and pull the depth of those
#idx vlaues -----

bb <- sf::st_bbox(
  c(xmin = xlim[1], xmax = xlim[2], ymin = ylim[1], ymax = ylim[2]),
  crs = sf::st_crs(4326)
)

bb # an sf bbox
area_sfc <- sf::st_as_sfc(bb) # rectangle polygon (sfc)
area_sf <- sf::st_as_sf(area_sfc) # same, as an sf data frame
grid <- area_sf |>
  sf::st_transform(32620) |>
  sf::st_make_grid(cellsize = 10000, square = FALSE, what = "corners") |>
  sf::st_as_sf() |>
  rename(
    geometry = x
  )
grid

ggplot() +
  geom_spatraster(data = bathy) +
  geom_sf(data = grid, alpha = 0.05) +
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


pts <- cbind(iniloc$lon, iniloc$lat)

cellId <- likelihood_spr |>
  map(
    ~ terra::cellFromXY(.x, pts)
  ) # linear cell index


vals <- likelihood_spr |>
  map(
    \(x) {
      tibble(
        likelihood = terra::extract(x, pts)[[1]],
        event = iniloc$event
      )
    }
  ) |>
  dplyr::bind_rows(.id = "day") |>
  mutate(
    day = as.numeric(day)
  ) # likelihood value(s) at those cells
vals


ggplot(data = vals, aes(x = day, y = likelihood, fill = event)) +
  geom_point(
    shape = 21,
    size = 3,
    alpha = 0.5,
    position = position_jitter(height = 0.03, seed = 4)
  ) +
  scale_x_continuous(breaks = seq(0, 205, 10)) +
  ggview::canvas(14, 11)
