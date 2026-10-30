  
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
usethis::use_data(bathy_matrix_download)
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
  


T <- dim(L)[1]
  row <- dim(g$lon)[1]
  col <- dim(g$lon)[2]
  # m <- 2 # 2 behaviour
  # pred <- array(0, dim = c(m, T, col, row))
  # phi <- array(0, dim = c(m, T, col, row))
  # phi[2, 1, , ] <- L[1, , ]
  # pred[2, 1, , ] <- L[1, , ]

  pred <- array(0, dim = c(T, col, row))
  phi <- array(0, dim = c(T, col, row))

  phi[1, , ] <- L[1, , ]
  pred[1, , ] <- L[1, , ]

  psi <- rep(0, T - 1)
  #K <- imager::as.cimg(K)
  # K1 <- imager::as.cimg(K1)
  # K2 <- imager::as.cimg(K2)
  for (t in 2:T) {
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

    pred[t, , ] <- imagine::convolution2D(phi[t - 1, , ], K) #q

    # in case of data update led to an empty matrix
    check.data.update <- sum(pred[t, , ] * L[t, , ], na.rm = T)
    # sumL = sum(L[t, , ])
    # na.rm=T has been added as there is still some na in the data # 170718
    sumL <- sum(L[t, , ], na.rm = T)
    if (sumL > 1e-06 & !check.data.update %in% c(0, NA, Inf, -Inf)) {
      if (maskL) {
        print("code for mask.L not available yet")
        stop()
        # post1 <- mask.L(pred.t = pred[1, t, , ], L.t = L[t,
        #                                                  , ], lon = g$lon[1, ], lat = g$lat[, 1], par0 = dim(K1)[1],
        #                 bound.thr = bound.thr, minBounds = minBounds)
        # post2 <- mask.L(pred.t = pred[2, t, , ], L.t = L[t,
        #                                                  , ], lon = g$lon[1, ], lat = g$lat[, 1], par0 = dim(K1)[1],
        #                 bound.thr = bound.thr, minBounds = minBounds)

        post <- mask.L(
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
        post <- pred[t, , ] * L[t, , ]
      }
    } else {
      # post1 <- pred[1, t, , ]
      # post2 <- pred[2, t, , ]
      post <- pred[t, , ]
    }
    # psi[t - 1] <- sum(as.vector(post1), na.rm = T) + sum(as.vector(post2), na.rm = T)
    psi[t - 1] <- sum(as.vector(post), na.rm = T)

    # phi[1, t, , ] <- post1/(psi[t - 1] + 1e-15)
    #  phi[2, t, , ] <- post2/(psi[t - 1] + 1e-15)
    phi[t, , ] <- post / (psi[t - 1]) # removed +1e-15 post/(psi[t - 1] + 1e-15)
  }

  # enforce recapture position
  # ajout pg 31.05.18
  # phi[ T, , ] <- L[T, , ]
  #pred[ T, , ] <- L[T, , ]

  list(phi = phi, pred = pred, psi = psi)
}