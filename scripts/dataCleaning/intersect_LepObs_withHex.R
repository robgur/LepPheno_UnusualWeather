library(data.table)
library(tidyverse)
library(sf)

## read in data
leps <- fread("data/munged/allDistinctLepObs.csv")

leps <- leps %>% 
  mutate(doy = yday(lubridate::as_date(eventDate)))

leps <- filter(leps, eventDate > 1)

#mapping north america
north_america_map <- rnaturalearth::ne_countries(country = c("United States of America", "Mexico", "Canada"),returnclass = "sf")

#morphed map to be equal area
north_america_map <-st_transform(north_america_map
                                 , sp::CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 
                     +units=m +no_defs"))
#making the grid
grids = st_make_grid(north_america_map
                     , cellsize = c(250000, 250000))

grids = st_join(st_as_sf(grids), north_america_map)

grids <- grids %>% 
  filter(!is.na(type)) %>% 
  distinct(geometry)

#added grid cells onto dataset
grids = mutate(st_sf(geometry = grids), id_cells = 1:n())

leps_sf <- leps %>%
  filter(!is.na(decimalLongitude),
         !is.na(decimalLatitude)) %>%
  st_as_sf(coords = c("decimalLongitude", "decimalLatitude"),
           crs = "+proj=longlat +datum=WGS84 +no_defs")

leps_sf <- st_transform(leps_sf,crs =st_crs(grids))
leps_sf <- st_join(leps_sf, grids)

leps_df <- st_drop_geometry(leps_sf)

write.csv(x = leps_df, file = "data/munged/LepsByGrid.csv", row.names = F)
