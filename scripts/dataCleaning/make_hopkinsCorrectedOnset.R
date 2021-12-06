library(tidyverse)
library(raster)
library(sf)

# read in mdf
mdf <- read.csv('data/LMM_Data/mdf.csv') %>% 
  dplyr::filter(validName != "Pyrgus albescens")	

# link to lon and lat of grids
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
grids_centroid <- st_centroid(grids) 
coords <- st_coordinates(grids_centroid) %>% 
  as.data.frame() %>% 
  dplyr::rename(lon = 1, lat = 2)

grids_coords <- cbind(st_drop_geometry(grids_centroid), coords)

mdf <- left_join(mdf, grids_coords)

# read in raster data
r <- raster("data/elevation_raster/elevation_NorthAmerica.tif")
rproj <- projectRaster(r, crs =  "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs ")
rdf <- raster::as.data.frame(rproj, xy = T) %>% 
  dplyr::rename("lon" = "x", "lat" = "y")

ggplot() + 
  geom_tile(rdf, mapping = aes(x = lon, y = lat, fill = elevation_NorthAmerica)) +
  geom_point(mdf, mapping = aes(x = lon, y = lat), color = "yellow") +
  scale_fill_viridis_c(na.value = "transparent") +
  theme_classic()

# join model df to elevation df
mdf2 <- left_join(mdf, rdf)

coordinates(mdf) <- ~ lon + lat
elev <- extract(rproj, mdf) # extract elevation data

mdf2 <- as.data.frame(mdf) %>% 
  dplyr::mutate(elevation = elev)

mdf2 <- mdf2 %>% 
  dplyr::mutate(elevation = if_else(
    condition = id_cells == 464, 
    true = 367.93410,
    false = elevation

  ))

ggplot() + 
  geom_tile(mdf2, mapping = aes(x = lon, y = lat, fill = elevation)) +
  scale_fill_viridis_c(trans = "log", na.value = "transparent") +
  theme_classic()

# find cell with lowest latitude
lowest_cell <- mdf2 %>% 
  filter(lat == min(mdf2$lat))

lowest_cell_lon <- lowest_cell$lon[1]
lowest_cell_lat <- lowest_cell$lat[1]
lowest_cell_elev <- lowest_cell$elevation[1]

mdf3 <- mdf2 %>% 
  mutate(lon_cor = lon - lowest_cell_lon,
         lat_cor = lat - lowest_cell_lat,
         elev_cor = round((elevation * 3.28) - (lowest_cell_elev * 3.28)))

# for every m North, phenology delays (4/100000) days
# for every m east, phenology delays (4/500000) days
# for every foot higher in elevation, phenology delays (1/100) days

mdf3 <- mdf3 %>% 
  mutate(hopkins_cor =  (lat_cor * -(4/100000)) + 
           (lon_cor * -(4/500000)) +
           (elev_cor * -(1/100)))

ggplot() + 
  geom_tile(mdf3, mapping = aes(x = lon, y = lat, fill = hopkins_cor)) +
  scale_fill_gradient2() +
  theme_classic()

# add hopkins correction to original onset estimates
mdf4 <- mdf3 %>% 
  dplyr::mutate(hopkins_onset = q5 + hopkins_cor,
                hopkins_offset = q95 + hopkins_cor,
                hopkins_peak = q50 + hopkins_cor)

write.csv(x = mdf4, file = "data/LMM_Data/mdf_hopkins_onset.csv", row.names = F)

spp_onset <- mdf4 %>% 
  group_by(validName) %>% 
  summarise(mean_onset = mean(q5), sd_onset = sd(q5))

ggplot() + 
  geom_histogram(spp_onset, mapping = aes(x = mean_onset), fill = "turquoise") +
  geom_vline(xintercept = c(79, 172, 265)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() 

## Can't use equinox to group species into traits. let's try quantiles
quantile(spp_onset$mean_onset, c(0.2,0.5, 0.8))

ggplot() + 
  geom_histogram(spp_onset, mapping = aes(x = mean_onset), fill = "turquoise") +
  geom_vline(xintercept = c(129, 183)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() 
