library(dplyr)
library(lubridate)


load(
  "/Users/benhlina/Library/CloudStorage/Dropbox/Dal-Post Doc/data/Geolocation for Jena/iniloc/iniloc.now.RData"
)

iniloc <- iniloc.all


iniloc <- iniloc |>
  dplyr::rename(
    fish_id = fish
  ) |>
  mutate(
    fish_id = as.character(fish_id)
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

bathy <- load_bathy_raster(
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

bath <- load_bathy_raster(
  "/Users/benhlina/Library/CloudStorage/Dropbox/Dal-Post Doc/data/Geolocation for Jena/bathy_grid/bathy.mat.gsl.ssgb.lon72.41.lat41.52.rds"
)


bathy_matrix <- readRDS(
  "/Users/benhlina/Library/CloudStorage/Dropbox/Dal-Post Doc/data/Geolocation for Jena/bathy_grid/bathy.mat.gsl.ssgb.lon72.41.lat41.52.rds"
)


usethis::use_data(bathy_matrix, overwrite = TRUE)
