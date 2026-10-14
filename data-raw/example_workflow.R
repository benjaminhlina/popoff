# ----- load packages ----

{
  library(dplyr)
  library(ggplot2)
  library(poppr)
  library(purrr)
  library(scales)
  library(sf)
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

bathy_utm <- project(bathy, "EPSG:32620", method = "bilinear")


# ----- next we can plot the spat raster using tidyterra ----
bathy_utm

ggplot() +
  geom_spatraster(data = bathy_utm) +
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
prepped_bathy <- prep_bathy_likdepth(bathy = bathy_utm, focal_dim = 3)
ggplot() +
  geom_spatraster(data = prepped_bathy$sdi) +
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
# we can view the output which is a list containing 1) bathy 2) sd of bathy
# 3) idx locations 4) idx_na locations, 4) focal dimensions, 5) if weighted
prepped_bathy
event_log_sf <- event_log |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

event_log_utms <- event_log_sf |>
  st_transform(32620) |>
  mutate(
    lon = st_coordinates(geometry)[, "X"],
    lat = st_coordinates(geometry)[, "Y"]
  ) |>
  st_drop_geometry()
event_log_utms
# ------ prep log for summary prior to likelihood calculation ------

prepped_log <- prep_summary_likdepth(
  raw_log = pdt,
  event_log = event_log_utms,
  bin = 3600
)

# pepared summary is a list containging 1) raw log, 2) summary table 3) a ggplot
# of the depth and temperature profiles, 4) time resoulation, 5) bin value

prepped_log

# ---- calculaite likeilhoood depth -----
likedepth_transmit <- calc_likdepth_transmit(
  prepped_log = prepped_log,
  prepped_bathy = prepped_bathy,
  event_log = event_log_utms,
  daily_sd_depth = 0.57
)
dim(likedepth_transmit)
str(likedepth_transmit)
# likedepth_transmit is a matrix with dims of bathymetry by the number
# of days observed in this case 203, each day is a new likelihood est.

# ----- we can plot the first couple of days -----

days <- unique(prepped_log$summary_table$date)

xlim <- prepped_bathy$xlim
ylim <- prepped_bathy$ylim


# ---- we can convert them into spatrasters using load functions
likelihood_spr <- days |>
  set_names() |>
  map(
    ~ load_raster(
      likedepth_transmit[,, as.character(.x)],
      # crs = "EPSG:4326",
      extent = c(xlim, ylim),
      crs = "EPSG:32620"
    ),
    # .progress = TRUE
  )

# ----- then polt them using ggplot and tidyterra -----
pl <- likelihood_spr[1:10] |>
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
pl
# # ---- Take idx values that equla 0.5 or greater and pull the depth of those
# #idx vlaues -----
# likelihood_spr[204]
# r <- likelihood_spr$`2022-02-28`

# terra::where.max(r) # cell number + value
# cell <- terra::where.max(r)[, "cell"]
# terra::xyFromCell(r, cell) # lon/lat of that cell

# Get the resolution of your bathy raster in degrees
# res_deg <- terra::res(bathy)[1] # e.g., 1/60 of a degree if res = 1

# # Get the approximate latitude of your study area center (e.g., Gulf of St. Lawrence ~ 48°N)
# # Or calculate it dynamically from your raster extent:
# mean_lat <- mean(terra::ext(bathy)[3:4])

# # Conversion factors:
# # 1 degree of latitude is ~111.32 km everywhere
# # 1 degree of longitude is ~111.32 * cos(latitude in radians) km
# km_per_deg_lat <- 111.32
# km_per_deg_lon <- 111.32 * cos(mean_lat * pi / 180)

# # Real-world size of 1 pixel in kilometers
# pixel_width_km <- res_deg * km_per_deg_lon
# pixel_height_km <- res_deg * km_per_deg_lat

# If 1 pixel ≈ 1.85 km:
# A 5 km movement scale = 5 / 1.85 ≈ 2.7 pixels
# A 50 km movement scale = 50 / 1.85 ≈ 27 pixels

sigma_pixels <- seq(2, 27, by = 1)

# # You can easily map them for your own reference:
# sigma_df <- data.frame(
#   sigma_pixels = sigma_pixels,
#   approx_km = sigma_pixels * pixel_width_km
# )
# sigma_df
# print(sigma_df)
# ggplot() +
#   geom_spatraster(data = bathy_utm) +
#   geom_sf(data = grid_sf, alpha = 0.05) +
#   coord_sf() +
#   scale_fill_gradientn(
#     colours = c(
#       "#031E40",
#       "#0B3C73",
#       "#17A2FF",
#       "#7FDBFF",
#       "#EAF6FF",
#       "#FFE8A3",
#       "#E6F38F",
#       "#9AD850"
#     ),
#     values = scales::rescale(c(-4000, -600, -150, -60, 0, 50, 150, 400)),
#     limits = c(-4000, 400),
#     oob = squish,
#     name = "Elevation / Depth (m)"
#   )

# pts <- cbind(event_log_$lon, event_log$lat)

# cell_id <- likelihood_spr |>
#   map(
#     ~ terra::cellFromXY(.x, pts)
#   )

# vals <- likelihood_spr |>
#   map(
#     \(x) {
#       tibble(
#         likelihood = terra::extract(x, pts)[[1]],
#         event = event_log$event
#       )
#     }
#   ) |>
#   dplyr::bind_rows(.id = "date") |>
#   mutate(
#     date = as.Date(date),
#     event = factor(event, levels = c("tag", "mrpat", "pop"))
#   )
# vals

# p <- ggplot(data = vals, aes(x = date, y = likelihood, fill = event)) +
#   geom_point(
#     shape = 21,
#     size = 3,
#     alpha = 0.5,
#     position = position_jitter(height = 0.03, seed = 4)
#   ) +
#   scale_fill_brewer(
#     name = "Event",
#     type = "div",
#     palette = 8,
#     labels = c("Tag", "Mrpat", "Pop")
#   ) +
#   scale_x_date(date_breaks = "2 weeks", date_labels = "%m-%d") +
#   theme_bw(
#     base_size = 15
#   ) +
#   theme(
#     legend.position = "inside",
#     legend.position.inside = c(0.95, 0.5),
#     legend.background = element_rect(colour = "black", linewidth = 0.5)
#   ) +
#   labs(
#     x = "Date",
#     y = "Likelihood"
#   )

# p

# p + ggview::canvas(14, 8.5)

sigma_val <- seq(0.25, 0.35, 0.015)
# sigma_val <- seq(0.1, 3, 0.1)

length(sigma_val)
sigma_val
gc()

# test <- eval_sigma(
#   sigma = sigma_val[1],
#   normalize = TRUE,
#   likelihood_array = likedepth_transmit,
#   mask_likelihood = TRUE,
#   return_maps = TRUE
# )

# sigma_val
# ?eval_sigma_iterate
out <- eval_sigma_iterate(
  workers = 6,
  normalize = TRUE,
  sigma_vec = sigma_val,
  mask_likelihood = TRUE,
  likelihood_array = likedepth_transmit
)

out

ggplot(data = out, aes(x = sigma, y = neg_log_lik)) +
  scale_y_log10() +
  geom_line(linewidth = 1.1) +
  geom_point()


optimiaized_sigma <- optimize_sigma(
  sigma_iterate_obj = out,
  workers = 6,
  likelihood_array = likedepth_transmit,
  normalize = TRUE,
  mask_likelihood = TRUE,
  .progress = TRUE
)

optimiaized_sigma
final <- eval_sigma(
  sigma = optimiaized_sigma$optimal_sigma,
  normalize = TRUE,
  likelihood_array = likedepth_transmit,
  mask_likelihood = TRUE,
  return_maps = TRUE
)

kerns <- gaussian_kernel(
  sigma = optimiaized_sigma$optimal_sigma,
  normalize = TRUE,
)


smoothed <- hmm_backward_smooth_old(
  forward_filter_obj = final,
  likelihood_array = likedepth_transmit,
  kernel = kerns
)
smoothed


smoothed_spr <- days |>
  set_names() |>
  map(
    ~ load_raster(
      smoothed[,, as.character(.x)],
      # crs = "EPSG:4326",
      extent = c(xlim, ylim),
      crs = "EPSG:32620"
    ),
    # .progress = TRUE
  )

# ----- then polt them using ggplot and tidyterra -----
pl <- smoothed_spr[1:10] |>
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
pl

max_per_day_sf <- max_per_day |>
  st_as_sf(coords = c("lon", "lat"), crs = 32620)


max_per_day_sf |>
  mapgl::maplibre_view()


p <- ggplot() +
  geom_sf(
    data = lennox::eastern_canada |>
      st_transform(4326)
  ) +
  geom_sf(
    data = max_per_day_sf |>
      st_transform(4326),
    aes(fill = max_value),
    shape = 21
  ) +
  scale_fill_viridis_c(option = "H", name = "Maximum\nlikelihood") +
  coord_sf(
    xlim = c(-62, -57),
    ylim = c(44, 48)
  ) +
  ggspatial::annotation_north_arrow(location = "tr") +
  ggspatial::annotation_scale() +
  theme_bw(
    base_size = 15
  ) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.92, 0.5),
    legend.background = element_rect(colour = "black")
  )

p + ggview::canvas(8.5, 10)
