library(dplyr)
library(lubridate)
library(terra)
library(poppr)
library(ggplot2)
library(tidyterra)

load(
  "/Users/benhlina/Library/CloudStorage/Dropbox/Dal-Post Doc/data/Geolocation for Jena/iniloc/iniloc.now.RData"
)

iniloc <- iniloc.all


iniloc <- iniloc |>
  dplyr::rename(
    fish_id = fish
  ) |>
  mutate(
    fish_id = as.character(fish_id),
    event = as.character(event),
    date = ymd_hms(date)
  ) |>
  as_tibble()
usethis::use_data(iniloc, overwrite = TRUE)

# ----- load psat data -----
# this is built for 1 tag but we can batch load using purrr
# same can be applied for iniloc data.

pdt <- readr::read_csv(
  "/Users/benhlina/Library/CloudStorage/Dropbox/Dal-Post Doc/data/Geolocation for Jena/TagData/Transmitted/199637/199637-Series.csv"
) |>
  janitor::clean_names()

glimpse(pdt)

pdt <- pdt |>
  mutate(
    day = dmy(day),
    date_time = ymd_hms(paste(day, time))
  ) |>
  relocate(date_time, .after = time) |>
  rename(
    fish_id = deploy_id,
    date = day,
  ) |>
  mutate(
    fish_id = as.character(fish_id)
  )


glimpse(pdt)

# ---- before summaries we need to trim the data based on when tagged and when popped offf ----
# we will use baser but this could be done in dplyr or data.table

# pdt <- pdt[
#   pdt$date <= as.Date(iniloc_all[iniloc_all$event == 'pop', 'date']) &
#     pdt$date >= as.Date(iniloc_all[iniloc_all$event == 'tag', 'date']),
# ]
# pdt

usethis::use_data(pdt, overwrite = TRUE)
load(
  "/Users/benhlina/Library/CloudStorage/Dropbox/Dal-Post Doc/data/Geolocation for Jena/bathy_grid/bathy.mat.gsl.ssgb.lon72.41.lat41.52.RData"
)
db

bathy <- load_raster(
  "/Users/benhlina/Library/CloudStorage/Dropbox/Dal-Post Doc/data/Geolocation for Jena/bathy_grid/bathy.mat.gsl.ssgb.lon72.41.lat41.52.rds"
)

load(
  "/Users/benhlina/Library/CloudStorage/Dropbox/Dal-Post Doc/data/Geolocation for Jena/bathy_grid/bathy.mat.gsl.ssgb.lon72.41.lat41.52.rds"
)

readr::write_rds(
  db,
  "/Users/benhlina/Library/CloudStorage/Dropbox/Dal-Post Doc/data/Geolocation for Jena/bathy_grid/bathy.mat.gsl.ssgb.lon72.41.lat41.52.rds"
)

pdt

bath <- load_raster(
  "/Users/benhlina/Library/CloudStorage/Dropbox/Dal-Post Doc/data/Geolocation for Jena/bathy_grid/bathy.mat.gsl.ssgb.lon72.41.lat41.52.rds"
)


bathy_matrix <- readRDS(
  "/Users/benhlina/Library/CloudStorage/Dropbox/Dal-Post Doc/data/Geolocation for Jena/bathy_grid/bathy.mat.gsl.ssgb.lon72.41.lat41.52.rds"
)


bathy_matrix |>
  str()

plot(bathy_matrix)


bathy_matrix_t <- t(bathy_matrix)[nrow(t(bathy_matrix)):1, ]

bathy <- load_raster(bathy_matrix_t, crs = "EPSG:4326")


# terra::ext(bathy) <- c(xmin(lon), xmax(lon), ymin(lat), ymax(lat))

bathy


terra::values(bathy) *
  -1 |>
    sort()

terra::crs(bathy) <- "EPSG:4326"

r <- terra::rast(bathy_matrix, crs = "EPSG:6500")
# r <- terra::t(r) # swap rows/cols
r <- terra::flip(r, direction = "vertical")
r <- terra::flip(r, direction = "horizontal")
r <- terra::flip(r, direction = "vertical")
r <- terra::flip(r, direction = "horizontal")
# r <- terra::rev(r)

plot(r)
bathy_matrix <- terra::as.matrix(r, wide = TRUE)

bathy_matrix
fields::image.plot(bathy_matrix)
ggplot() +
  geom_spatraster(data = r) +
  coord_sf()


usethis::use_data(bathy_matrix, overwrite = TRUE)


"/Users/benhlina/Library/CloudStorage/Dropbox/Dal-Post Doc/
  data/Geolocation for Jena/bathy_grid"

list_files <- list.files(
  here::here(
    "..",
    "..",
    "Dal-Post Doc",
    "data",
    "Geolocation for Jena",
    "bathy_grid"
  ),
  pattern = "*.RData",
  full.names = TRUE
)

list_files

temp_env <- new.env()

load_files <- list_files |>
  purrr::set_names() |>
  purrr::map(
    ~ load(.x, envir = temp_env)
  )

loaded_data <- as.list(temp_env)

loaded_data$grid.srtm |>
  dplyr::glimpse()

loaded_data$grid.gsl.ssgb |>
  dplyr::glimpse()

loaded_data$grid.peter |>
  dplyr::glimpse()

library(dplyr)
library(tidyr)


lon <- loaded_data$grid.peter$lon |>
  dplyr::as_tibble() |>
  mutate(
    id = row_number()
  )
lat <- loaded_data$grid.peter$lat |>
  dplyr::as_tibble() |>
  mutate(
    id = row_number()
  )


lon_longer <- lon |>
  pivot_longer(col = -id, names_to = "grid_no", values_to = "lon")
lon_long

lat_longer <- lat |>
  pivot_longer(col = -id, names_to = "grid_no", values_to = "lat")
lat_long

grid <- lon_longer |>
  left_join(lat_longer)

library(sf)
grid_sf <- grid |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326)


mapgl::maplibre_view(grid_sf)
