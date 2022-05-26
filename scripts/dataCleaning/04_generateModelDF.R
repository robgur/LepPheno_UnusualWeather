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