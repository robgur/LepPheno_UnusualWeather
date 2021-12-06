library(tidyverse)
library(lme4)
library(lmerTest)
library(car)
library(sjPlot)
library(MuMIn)
library(sf)


## read in estimated phenometrics
pheno <- read.csv("outputs/lep_phenometrics_noNAs.csv") %>% 
  rename(id_cells = "HEXcell") %>% 
  mutate(gcode = paste(code, year, id_cells, sep = "."))


## read in spp_list
spp_list <- read.csv("data/traits/spp_list_traits_mike_withValidNames.csv")

names_df <- spp_list %>% 
  dplyr::select(scientific_name, Syn, validName)

# join to names by scientific_name
pheno_nameJoin <- left_join(pheno, names_df, 
                           by = c("code" = "scientific_name")) %>% 
  dplyr::distinct(gcode, .keep_all = T)

pheno_nameJoin2 <- left_join(pheno_nameJoin, names_df, 
                            by = c("code" = "Syn")) %>% 
  dplyr::distinct(gcode, .keep_all = T)

pheno_validNames <- pheno_nameJoin2 %>% 
  mutate(validName = if_else(condition = is.na(validName.x),
                             true = validName.y,
                             false = validName.x))

pheno_validNames <- pheno_validNames %>% 
  filter(!is.na(validName)) %>% 
  dplyr::select(-validName.x, -validName.y)

rm(pheno_nameJoin, pheno_nameJoin)
## remove tropical cells
tc <- read.csv("data/CellsWithKoeppenTropics.csv")

pheno_validNames_noTropics <- pheno_validNames %>% 
  filter(!id_cells %in% tc$id_cells)

## read in climate data
clim <- read.csv("data/climate/chelsa_climateData.csv")
pheno_clim <- left_join(pheno_validNames_noTropics, clim)

## group to traits
traits <- read.csv("data/traits/spp_list_traits_mike_withValidNames.csv") %>% 
  distinct(validName, .keep_all = T)

pheno_clim_traits <- left_join(pheno_clim, traits, by = "validName") %>% 
  mutate(decade = year - year %% 10) 

## remove spp with multiple overwintering stages
pheno_clim_traits <- pheno_clim_traits %>% 
  filter(overwinteringStage != "L or P") %>% 
  filter(overwinteringStage != "P or L")

# remove NA is climate 
pheno_clim_traits <- pheno_clim_traits %>% 
  filter(!is.na(annualTemp),
         !is.na(tempSeas),
         !is.na(annualPrec),
         !is.na(precSeas)) %>% 
  dplyr::select(-code, -Syn.x, -scientific_name.x, - scientific_name.y,
                -notes.x, -onsetProblems, -offsetProblems, 
                -overwintering.stage.y, - voltinism.y, - notes.y, -Syn.y)

## save mdf
write.csv(x = pheno_clim_traits, "data/LMM_Data/mdf.csv", row.names = F)

## where is the modelable data?
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

## create some decade plots
dd_sum <- pheno_clim_traits %>% 
  group_by(id_cells, decade) %>% 
  summarise(nPhenometrics = n(),
            nSpp = length(unique(validName)))

dd_sum_grid <- left_join(grids, dd_sum) %>% 
  filter(!is.na(nPhenometrics)) %>% 
  st_sf()

ggplot() +
  geom_sf(north_america_map, mapping = aes()) +
  geom_sf(dd_sum_grid, mapping = aes(fill = nSpp)) +
  scale_fill_viridis_c() +
  facet_wrap(~decade) +
  theme_classic()
