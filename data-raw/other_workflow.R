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
  res = 4
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

likelihood_array <- likedepth_transmit

likelihood_array <- aperm(likelihood_array)
mask_likelihood <- FALSE
progress <- TRUE
# ----- get dims and dimnames --------
dims <- dim(likelihood_array)
dims

dn <- dimnames(likelihood_array)
dn

kernel <- gaussian_kernel(
  sigma = 0.3,
  normalize = TRUE,
)
# ---- create timesteps and rows and columns -----
time_step <- dims[3]

row <- dims[2]
col <- dims[1]
# m <- 2 # 2 behaviour
# pred <- array(0, dim = c(m, T, col, row))
# phi <- array(0, dim = c(m, T, col, row))
# phi[2, 1, , ] <- L[1, , ]
# pred[2, 1, , ] <- L[1, , ]

# ----- create arrays and starting location ------
pred <- array(0, dim = c(col, row, time_step), dimnames = dn)
phi <- pred

phi[,, 1] <- likelihood_array[,, 1]
pred[,, 1] <- likelihood_array[,, 1]

psi <- rep(0, time_step - 1)
#K <- imager::as.cimg(K)
# K1 <- imager::as.cimg(K1)
# K2 <- imager::as.cimg(K2)
popoff:::msg_progress(
  progress = TRUE,
  n_time = time_step,
  name = "HMM forward filtering",
  type = "Filtering"
)

for (t in 2:time_step) {
  # p1 = imager::as.cimg(t(phi[1, t - 1, , ]))
  # p2 = imager::as.cimg(t(phi[2, t - 1, , ]))
  # q1 = imager::convolve(p1, K1)
  # q2 = imager::convolve(p2, K2)
  # q1 = t(as.matrix(q1))
  # q2 = t(as.matrix(q2))
  # pred[1, t, , ] <- P[1, 1] * q1 + P[2, 1] * q2
  # pred[2, t, , ] <- P[1, 2] * q1 + P[2, 2] * q2

  #p = imager::as.cimg(t(phi[ t - 1, , ]))
  #q = imager::convolve(p, K)
  #q = t(as.matrix(q))

  pred[,, t] <- imagine::convolution2D(phi[,, t - 1], kernel) #q

  # in case of data update led to an empty matrix
  check.data.update <- sum(pred[,, t] * likelihood_array[,, t], na.rm = TRUE)
  # sumL = sum(L[t, , ])
  # na.rm=T has been added as there is still some na in the data # 170718
  sumL <- sum(likelihood_array[,, t], na.rm = TRUE)
  if (sumL > 1e-06 & !check.data.update %in% c(0, NA, Inf, -Inf)) {
    if (mask_likelihood) {
      print("code for mask.L not available yet")
      stop()
      # post1 <- mask.L(pred.t = pred[1, t, , ], L.t = L[t,
      #                                                  , ], lon = g$lon[1, ], lat = g$lat[, 1], par0 = dim(K1)[1],
      #                 bound.thr = bound.thr, minBounds = minBounds)
      # post2 <- mask.L(pred.t = pred[2, t, , ], L.t = L[t,
      #                                                  , ], lon = g$lon[1, ], lat = g$lat[, 1], par0 = dim(K1)[1],
      #                 bound.thr = bound.thr, minBounds = minBounds)

      post <- mask_likelihood_array(
        pred.t = pred[t, , ],
        L.t = L[t, , ],
        lon = g$lon[1, ],
        lat = g$lat[, 1],
        par0 = dim(K)[1],
        bound.thr = bound.thr,
        minBounds = minBounds
      )
    } else {
      # post1 <- pred[1, t, , ] * L[t, , ]
      # post2 <- pred[2, t, , ] * L[t, , ]
      post <- pred[,, t] * likelihood_array[,, t]
    }
  } else {
    # post1 <- pred[1, t, , ]
    # post2 <- pred[2, t, , ]
    post <- pred[,, t]
  }
  # psi[t - 1] <- sum(as.vector(post1), na.rm = T) + sum(as.vector(post2), na.rm = T)
  psi[t - 1] <- sum(as.vector(post), na.rm = TRUE)

  # phi[1, t, , ] <- post1/(psi[t - 1] + 1e-15)
  #  phi[2, t, , ] <- post2/(psi[t - 1] + 1e-15)
  phi[,, t] <- post / (psi[t - 1]) # removed +1e-15 post/(psi[t - 1] + 1e-15)

  if (isTRUE(progress)) {
    cli::cli_progress_update()
  }
}

# enforce recapture position
# ajout pg 31.05.18
# phi[ T, , ] <- L[T, , ]
#pred[ T, , ] <- L[T, , ]

list(phi = phi, pred = pred, psi = psi)

str(phi)

nllf. <- -sum(log(psi))

nllf.

nllf_2 <- -sum(log(psi[psi > 0]))
nllf_2

# library(terra)
# library(dplyr)
# library(purrr)
# library(tibble)

foward_obj <- days |>
  set_names() |>
  map(
    ~ load_raster(
      phi[,, as.character(.x)],
      # crs = "EPSG:4326",
      extent = c(xlim, ylim),
      crs = "EPSG:32620"
    ),
    # .progress = TRUE
  )


max_locs <- imap_dfr(foward_obj, function(r, date) {
  # cell index of max value (first one if ties)
  wm <- where.max(r, value = TRUE) # cols: layer, cell, value
  xy <- xyFromCell(r, wm[, "cell"])

  # reproject UTM 20N (EPSG:32620) -> WGS84 lon/lat
  pts <- vect(xy, type = "points", crs = crs(r))
  ll <- crds(project(pts, "EPSG:4326"))

  tibble(
    date = as.Date(date),
    lon = ll[, 1],
    lat = ll[, 2],
    value = wm[, "value"]
  )
})


max_locs
smoothed_per_day_new_sf <- max_locs |>
  # filter(!is.na(lon)) |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326)


track_new <- smoothed_per_day_new_sf |>
  mutate(date = as.Date(date)) |> # Ensure proper date sorting
  arrange(date) |> # Sort chronologically (June 18 -> 19 -> 20...)
  st_combine() |>
  st_cast("LINESTRING")


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
    data = track_new |>
      st_transform(4326),
    # filter(date %in% as.Date(event_log$date)),
  ) +
  geom_sf(
    data = smoothed_per_day_new_sf |>
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
