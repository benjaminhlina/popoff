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

bathy_matrix_download <- download_noaa_bathy(
  xlim = c(-72.41, -40.0),
  ylim = c(41.52, 55.0),
  res = 1
)
bathy_matrix_download
# in this case Gulf of St. Lawrence

# we can check that it is SpatRaster

bathy <- load_raster(bathy_matrix_download, crs = "EPSG:4326")

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
prepped_bathy <- prep_bathy_likdepth(bathy = bathy, focal_dim = 3)

# we can view the output which is a list containing 1) bathy 2) sd of bathy
# 3) idx locations 4) idx_na locations, 4) focal dimensions, 5) if weighted
prepped_bathy
# ------ prep log for summary prior to likelihood calculation ------
prepped_log <- prep_summary_likdepth(
  raw_log = pdt,
  event_log = event_log,
  bin = 3600
)

# pepared summary is a list containging 1) raw log, 2) summary table 3) a ggplot
# of the depth and temperature profiles, 4) time resoulation, 5) bin value

prepped_log

# ---- calculaite likeilhoood depth -----
likedepth_transmit <- calc_likdepth_transmit(
  prepped_log = prepped_log,
  prepped_bathy = prepped_bathy,
  event_log = event_log,
  daily_sd_depth = 0.57
)
dim(likedepth_transmit)
# likedepth_transmit is a matrix with dims of bathymetry by the number
# of days observed in this case 203, each day is a new likelihood est.

# ----- we can plot the first couple of days -----

days <- unique(prepped_log$summary_table$date)

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


# ----- then polt them using ggplot and tidyterra -----
likelihood_spr[204] |>
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
likelihood_spr[204]
r <- likelihood_spr$`2022-02-28`

terra::where.max(r) # cell number + value
cell <- terra::where.max(r)[, "cell"]
terra::xyFromCell(r, cell) # lon/lat of that cell

event_log

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


pts <- cbind(event_log$lon, event_log$lat)

cell_id <- likelihood_spr |>
  map(
    ~ terra::cellFromXY(.x, pts)
  )


vals <- likelihood_spr |>
  map(
    \(x) {
      tibble(
        likelihood = terra::extract(x, pts)[[1]],
        event = event_log$event
      )
    }
  ) |>
  dplyr::bind_rows(.id = "date") |>
  mutate(
    date = as.Date(date),
    event = factor(event, levels = c("tag", "mrpat", "pop"))
  )
vals


p <- ggplot(data = vals, aes(x = date, y = likelihood, fill = event)) +
  geom_point(
    shape = 21,
    size = 3,
    alpha = 0.5,
    position = position_jitter(height = 0.03, seed = 4)
  ) +
  scale_fill_brewer(
    name = "Event",
    type = "div",
    palette = 8,
    labels = c("Tag", "Mrpat", "Pop")
  ) +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%m-%d") +
  theme_bw(
    base_size = 15
  ) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.95, 0.5),
    legend.background = element_rect(colour = "black", linewidth = 0.5)
  ) +
  labs(
    x = "Date",
    y = "Likelihood"
  )

p

p + ggview::canvas(14, 8.5)
