library(tidyverse)

## read in data
mdf_sd_pruned <- read.csv("data/LMM_Data/mdf_removeOutliersResiduals.csv") %>% 
  dplyr::filter(overwinteringStage != "None")

## remove species and cells without at least three estimates
spp_sum <- mdf_sd_pruned %>% 
  group_by(validName) %>% 
  summarise(count = n())

cell_sum <- mdf_sd_pruned %>% 
  group_by(id_cells) %>% 
  summarise(count = n())

enough_cells <- dplyr::filter(cell_sum, count >= 2)

mdf_sd_pruned2 <- mdf_sd_pruned %>% 
  dplyr::filter(id_cells %in% enough_cells$id_cells)

spp_sum <- mdf_sd_pruned2 %>% 
  group_by(validName) %>% 
  summarise(count = n())

enough_spp <- dplyr::filter(spp_sum, count >= 2)

mdf_sd_pruned3 <- mdf_sd_pruned2 %>% 
  dplyr::filter(validName %in% enough_spp$validName)

## now make sure we're good on cells still
cell_sum <- mdf_sd_pruned3 %>% 
  group_by(id_cells) %>% 
  summarise(count = n())

## one needs to go -- cell number 75
mdf_sd_pruned3 <- mdf_sd_pruned3 %>% 
  filter(id_cells != 75)

#save the mdf after making seasonality trait
# lets' read in FFP and see if controlling for FFP changes seasonality trait
library(terra)
library(sf)
ffp <- rast("data/FFP/Normal_1991_2020_bFFP.tif")
ffp <- terra::aggregate(ffp, 25)
ffp <- project(ffp, y = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
ffp_df <- raster::as.data.frame(ffp, xy = T) %>% 
  filter(!is.na(bFFP))
ffp_sf <- st_as_sf(ffp_df, coords = c("x","y"), 
    crs = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
# get this at a cell basis
library(sf)
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

#fishnet begin of ffp to id_cells
stj <- st_join(ffp_sf, grids)

gb_sf <- stj %>% 
  group_by(id_cells) %>% 
  summarise(mean_bFFP = base::mean(bFFP, na.rm = T))

ggplot() +
  geom_sf(gb_sf, mapping = aes(color = mean_bFFP))

ffp_df <- st_drop_geometry(gb_sf)

# calculate number of days butterflies are flying after beginning of frost free period
mdf_sd_pruned4 <- left_join(mdf_sd_pruned3, gb_sf, by = "id_cells") %>% 
  mutate(butt_delay = q5 - mean_bFFP)

spp_ave_onset <- mdf_sd_pruned4 %>% 
  group_by(validName) %>% 
  summarize(meanDelay = mean(butt_delay, na.rm = T))

## Can't use equinox to group species into traits. let's try quantiles
quantile(spp_ave_onset$meanDelay, c(0.15,0.5, 0.85))

ggplot() + 
  geom_histogram(spp_ave_onset, mapping = aes(x = meanDelay), fill = "turquoise") +
  geom_vline(xintercept = c(-1.12, 68.4)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() 



mdf_sd_pruned4 <- left_join(mdf_sd_pruned4, spp_ave_onset) %>% 
  mutate(Seas = case_when(
    meanDelay < -1.12 ~ "Spring",
    meanDelay >= -1.12 & meanDelay <= 68.4 ~ "Summer",
    meanDelay > 68.4 ~ "Fall"
  ))

# REMOVE rows for mdf
mdf_sd_pruned4 <- mdf_sd_pruned4 %>% 
  dplyr::select(-residual, -mean_bFFP, -geometry, -butt_delay, -meanDelay)

gdf <- data.frame(id_cells = grids$id_cells,
                  lon = st_coordinates(st_centroid(grids))[,1],
                  lat = st_coordinates(st_centroid(grids))[,2])

mdf_sd_pruned5 <- left_join(mdf_sd_pruned4, gdf, by = "id_cells")

write.csv(x = mdf_sd_pruned5,
          file = "data/LMM_Data/mdf_removeOutliersResiduals_wSeasonalityTrait.csv",
          row.names = F)
