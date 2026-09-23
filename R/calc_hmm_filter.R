# filter deal NA
# maskL = F; bound.thr = 0.1 ; minBounds = 10
hmm_filter_1B_deal_NA <- function(
  grid,
  likelihood_array,
  kernal,
  bound_threshold = 0.1,
  bound_min = 10,
  mask_likelihood = FALSE
) {
  n_time <- dim(likelihood_array)[1]
  n_row <- dim(grid$lon)[1]
  n_col <- dim(grid$lon)[2]
  # m <- 2 # 2 behaviour
  # pred <- array(0, dim = c(m, T, col, row))
  # phi <- array(0, dim = c(m, T, col, row))
  # phi[2, 1, , ] <- L[1, , ]
  # pred[2, 1, , ] <- L[1, , ]

  predicted <- array(0, dim = c(n_time, n_col, n_row))
  phi <- array(0, dim = c(n_time, n_col, n_row))

  phi[1, , ] <- likelihood_array[1, , ]
  pred[1, , ] <- likelihood_array[1, , ]

  psi <- rep(0, n_time - 1)
  #K <- imager::as.cimg(K)
  # K1 <- imager::as.cimg(K1)
  # K2 <- imager::as.cimg(K2)
  for (t in 2:n_time) {
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

    predicted[t, , ] <- imagine::convolution2D(phi[t - 1, , ], kernal) #q

    # in case of data update led to an empty matrix
    check_data_update <- sum(
      predicted[t, , ] * likelihood_array[t, , ],
      na.rm = TRUE
    )
    # sumL = sum(L[t, , ])
    # na.rm=T has been added as there is still some na in the data # 170718
    sum_likelihood <- sum(likelihood_array[t, , ], na.rm = TRUE)
    if (sum_likelihood > 1e-06 & !check_data_update %in% c(0, NA, Inf, -Inf)) {
      if (isTRUE(mask_likelihood)) {
        # post1 <- mask.L(pred.t = pred[1, t, , ], L.t = L[t,
        #                                                  , ], lon = grid$lon[1, ], lat = grid$lat[, 1], par0 = dim(K1)[1],
        #                 bound.thr = bound.thr, minBounds = minBounds)
        # post2 <- mask.L(pred.t = pred[2, t, , ], L.t = L[t,
        #                                                  , ], lon = grid$lon[1, ], lat = grid$lat[, 1], par0 = dim(K1)[1],
        #                 bound.thr = bound.thr, minBounds = minBounds)

        post <- mask_likelihood_array(
          pred.t = pred[t, , ],
          L.t = L[t, , ],
          lon = grid$lon[1, ],
          lat = grid$lat[, 1],
          par0 = dim(K)[1],
          bound.thr = bound.thr,
          minBounds = minBounds
        )
      } else {
        # post1 <- pred[1, t, , ] * L[t, , ]
        # post2 <- pred[2, t, , ] * L[t, , ]
        post <- predicted[t, , ] * likelihood_array[t, , ]
      }
    } else {
      # post1 <- pred[1, t, , ]
      # post2 <- pred[2, t, , ]
      post <- predicted[t, , ]
    }
    # psi[t - 1] <- sum(as.vector(post1), na.rm = T) + sum(as.vector(post2), na.rm = T)
    psi[t - 1] <- sum(as.vector(post), na.rm = TRUE)

    # phi[1, t, , ] <- post1/(psi[t - 1] + 1e-15)
    #  phi[2, t, , ] <- post2/(psi[t - 1] + 1e-15)
    phi[t, , ] <- post / (psi[t - 1]) # removed +1e-15 post/(psi[t - 1] + 1e-15)
  }

  # enforce recapture position
  # ajout pg 31.05.18
  # phi[ T, , ] <- L[T, , ]
  #pred[ T, , ] <- L[T, , ]

  out <- list(phi = phi, predicted = predicted, psi = psi)
  return(out)
}
