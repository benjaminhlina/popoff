# ----- load packages ----

{
  library(dplyr)
  library(ggplot2)
  library(popoff)
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

# bathy_matrix_download
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
# ggplot() +
#   geom_spatraster(data = prepped_bathy$sdi) +
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
# dim(likedepth_transmit)
# str(likedepth_transmit)
# # likedepth_transmit is a matrix with dims of bathymetry by the number
# # of days observed in this case 203, each day is a new likelihood est.

# # ----- we can plot the first couple of days -----

days <- unique(prepped_log$summary_table$date)

xlim <- prepped_bathy$xlim
ylim <- prepped_bathy$ylim

# # ---- we can convert them into spatrasters using load functions
# likelihood_spr <- days |>
#   set_names() |>
#   map(
#     ~ load_raster(
#       likedepth_transmit[,, as.character(.x)],
#       # crs = "EPSG:4326",
#       extent = c(xlim, ylim),
#       crs = "EPSG:32620"
#     ),
#     # .progress = TRUE
#   )

# bbox_latlon <- st_polygon(list(matrix(
#   c(
#     -59.5,
#     44.1,
#     -58.5,
#     44.1,
#     -58.5,
#     44.9,
#     -59.5,
#     44.9,
#     -59.5,
#     44.1
#   ),
#   ncol = 2,
#   byrow = TRUE
# )))

# bbox_sf <- st_sfc(bbox_latlon, crs = 4326)

# # 2. Transform to UTM Zone 20N (EPSG:32620)
# bbox_utm <- st_transform(bbox_sf, crs = 32620)

# # 3. Extract the new bounding box coordinates (in meters)
# range <- st_bbox(bbox_utm)
# range[1]
# # ----- then polt them using ggplot and tidyterra -----
# pl <- likelihood_spr[1:10] |>
#   imap(
#     ~ ggplot() +
#       geom_spatraster(data = .x) +
#       geom_sf(
#         data = event_log_sf |>
#           mutate(date = as.Date(date)),
#         aes(colour = event),
#       ) +
#       # coord_sf() +
#       scale_fill_viridis_c(option = "H") +
#       theme_bw(
#         base_size = 15
#       ) +
#       coord_sf(
#         xlim = c(range[1], range[3]),
#         ylim = c(range[2], range[4])
#       ) +
#       theme(
#         panel.grid = element_blank(),
#         plot.title = element_text(hjust = 0.5)
#       ) +
#       labs(title = .y)
#   )
# # pl

# sigma_val <- seq(0.25, 0.35, 0.015)
sigma_val <- seq(0.2, 0.3, 0.01)

length(sigma_val)
sigma_val
gc()


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
optimiaized_sigma$optimal_sigma

kernel <- gaussian_kernel(
  sigma = 0.3,
  normalize = TRUE,
)

final <- hmm_forward_filter(
  kernel = kernel,
  likelihood_array = likedepth_transmit,
  mask_likelihood = TRUE,
  return_maps = TRUE
)
# mass <- apply(final$normalised_posterior, 3, function(x) sum(x, na.rm = TRUE))
# which(!is.finite(mass) | mass <= 0) # dead filtered posteriors
# final$n_bad_days # should be > 0 if this is the cause
# final$n_updated

# final$n_bad_days # should now be small, not 79
# which(final$ncps <= 0 | !is.finite(final$ncps)) + 1 # steps that genuinely collapsed
# mass <- apply(final$normalised_posterior, 3, function(x) sum(x, na.rm = TRUE))
# all(is.finite(mass) & mass > 0) # must be TRUE before smoothing

# # 1. Does convolution2D spread NA?
# m <- matrix(1, 15, 15)
# m[8, 8] <- NA
# sum(is.na(imagine::convolution2D(m, kerns))) # 1 = fine, >1 = NA dilates

# # 2. Where does pred die? (needs return_maps = TRUE)
# pm <- apply(final$pc_predicted, 3, function(x) sum(x, na.rm = TRUE))
# nna <- apply(final$pc_predicted, 3, function(x) sum(is.na(x)))
# data.frame(step = 120:130, pred_mass = pm[120:130], n_na = nna[120:130])
# sum(is.na(likedepth_transmit[,, 1])) # land cells, for comparison

smoothed <- hmm_backward_smooth(
  forward_filter_obj = final
)
# smoothed

# ------ exxtract ------
# max_per_day <- extract_locations(obj = final$normalised_posterior)

# max_per_day_sf <- max_per_day |>
#   st_as_sf(coords = c("lon", "lat"), crs = 32620)

smoothed_per_day <- extract_locations(smoothed)

smoothed_per_day
# max_per_day_sf |>
#   mapgl::maplibre_view()
smoothed_per_day_sf <- smoothed_per_day |>
  # filter(!is.na(lon)) |>
  st_as_sf(coords = c("lon", "lat"), crs = 32620)


track <- smoothed_per_day_sf |>
  mutate(date = as.Date(date)) |> # Ensure proper date sorting
  arrange(date) |> # Sort chronologically (June 18 -> 19 -> 20...)
  st_combine() |>
  st_cast("LINESTRING")


st_length(track)
51206.16 / 1000

smoothed_per_day
# max_per_day_sf |>
#   mapgl::maplibre_view()
smoothed_per_day_new_sf <- smoothed_per_day_new |>
  filter(!is.na(lon)) |>
  st_as_sf(coords = c("lon", "lat"), crs = 32620)


track_new <- smoothed_per_day_new_sf |>
  mutate(date = as.Date(date)) |> # Ensure proper date sorting
  arrange(date) |> # Sort chronologically (June 18 -> 19 -> 20...)
  st_combine() |>
  st_cast("LINESTRING")


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


# st_length(track)
# 51206.16 / 1000

p <- ggplot() +
  geom_sf(
    data = lennox::eastern_canada |>
      st_transform(4326)
  ) +
  # geom_spatraster(data = bathy_utm) +
  # scale_fill_gradientn(
  #   colours = c(
  #     "#031E40",
  #     "#0B3C73",
  #     "#17A2FF",
  #     "#7FDBFF",
  #     "#EAF6FF",
  #     "#FFE8A3",
  #     "#E6F38F",
  #     "#9AD850"
  #   ),
  #   values = scales::rescale(c(-4000, -600, -150, -60, 0, 50, 150, 400)),
  #   limits = c(-4000, 400),
  #   oob = squish,
  #   name = "Elevation / Depth (m)"
  # ) +
  geom_sf(
    data = track |>
      st_transform(4326),
    # filter(date %in% as.Date(event_log$date)),
  ) +
  geom_sf(
    data = smoothed_per_day_sf |>
      st_transform(4326),
    # filter(date %in% as.Date(event_log$date)),
    aes(fill = value),
    shape = 21
  ) +

  # geom_sf(
  #   data = m_sf |>
  #     st_transform(4326) |>
  #     filter(date %in% as.Date(event_log$date)),
  #   aes(fill = value)
  # ) +
  geom_sf(
    data = event_log_sf |>
      mutate(date = as.Date(date)),
    aes(colour = event),
  ) +

  # facet_wrap(~date) +
  scale_fill_viridis_c(option = "H", name = "Maximum\nlikelihood") +
  # coord_sf(
  #   xlim = c(-62, -57),
  #   ylim = c(44, 48)
  # ) +
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


# p + ggview::canvas(8.5, 10)

p1 <- p +
  coord_sf(
    xlim = c(-59.5, -58.5),
    ylim = c(44.1, 44.9)
  )

p
p1
p1 + ggview::canvas(8.5, 10)
p + ggview::canvas(8.5, 10)

load(
  "/Users/benhlina/Library/CloudStorage/Dropbox/Dal-Post Doc/data/Geolocation for Jena/trajectories/backfil.transmit.depth_TRUE.199637.RData"
)

load(
  "/Users/benhlina/Library/CloudStorage/Dropbox/Dal-Post Doc/data/Geolocation for Jena/trajectories/backfil.transmit.AR1.depth_TRUE.199637.RData"
)


# Loop through each of the 256 slices along the third dimension
max_results <- do.call(
  rbind,
  lapply(1:dim(s_right)[3], function(k) {
    slice <- s_right[,, k]

    # Check if the entire slice is NA to avoid errors
    if (all(is.na(slice))) {
      return(data.frame(slice = k, max_val = NA, row = NA, col = NA))
    }

    # Find the maximum value (ignoring NAs)
    max_val <- max(slice, na.rm = TRUE)

    # Find the row and column coordinates of the maximum value
    # arr.ind = TRUE ensures it returns row and column indices instead of a flat vector index
    indices <- which(slice == max_val, arr.ind = TRUE)

    # If multiple cells share the exact same maximum value, this takes the first occurrence
    data.frame(
      slice = k,
      max_val = max_val,
      row = indices[1, 1],
      col = indices[1, 2]
    )
  })
)

# View the first few slices of the resulting table
head(max_results)

max_results |>
  print()
summary(max_results)
