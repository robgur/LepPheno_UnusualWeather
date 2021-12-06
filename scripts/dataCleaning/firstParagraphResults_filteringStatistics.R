library(data.table)
library(tidyverse)
library(splitstackshape)
library(pbapply)

## read in data
leps <- fread("data/munged/LepsByGrid.csv") %>% 
  mutate(binomial = paste(genus, specificEpithet, sep = " "))

## only adult data -- to the best of our knowledge
notAdults <- c("Nymph", "Immature", "Pupa", "Larva", "Egg",
               "Caterpillar", "Juvenile", "juvenile; larva", 
               "1st instar", "Larval", "larval", "pupal case",
               "Larvae", "pupa", "Pupae", "egg", "juvenile;larva",
               "Inmaduro", "juvenile", "Exuviae", "exuviae", "Egg sac",
               "Cocoon", "cocoon", "larva", "1 exuvia w/egg", "Pupal case + larva exuvia",
               "Case", "Larval?", "larval case", "larvae", "larval in sac", "Pupa case")
leps <- leps %>% 
  filter(!lifeStage %in% notAdults)

ggplot(leps, mapping = aes(x = doy)) +
  geom_histogram(bins = 365) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw()

# lines were used to visualize what the data looked like for first of years
#doy1 <- filter(leps, doy == 1)
#doy32 <- filter(leps, doy == 32) 

spp_list <- read.csv("data/traits/spp_list_traits_mike_withValidNames.csv")
#spp_list <- spp_list$scientific_name

names_df <- spp_list %>% 
  dplyr::select(scientific_name, Syn, validName)

## filter lep obs to spp_list
leps <- leps %>% 
  #  filter(binomial %in% spp_list) %>% 
  mutate(gcode = paste(binomial, year, id_cells, sep = ".")) %>% 
  filter(!is.na(doy)) %>% 
  filter(day > 1) 


# join to names by scientific_name
leps_nameJoin <- left_join(leps, names_df, 
                           by = c("binomial" = "scientific_name")) %>% 
  dplyr::distinct(id, .keep_all = T)

leps_nameJoin2 <- left_join(leps_nameJoin, names_df, 
                            by = c("binomial" = "Syn")) %>% 
  dplyr::distinct(id, .keep_all = T)

leps_validNames <- leps_nameJoin2 %>% 
  mutate(validName = if_else(condition = is.na(validName.x),
                             true = validName.y,
                             false = validName.x))

leps_validNames <- leps_validNames %>% 
  filter(!is.na(validName)) %>% 
  dplyr::select(-validName.x, -validName.y)

### read in mdf ###
mdf <- read.csv("data/LMM_Data/mdf_pruned_hopkins.csv")

leps_validNames2 <- leps_validNames %>% 
  filter(validName %in% mdf$validName)

unique(leps_validNames2$validName) %>% length()
unique(leps_validNames2$Syn) %>% length()
unique(c(leps_validNames$Syn, leps_validNames2$scientific_name)) %>% length()


## how many spp in original trait database
t <- read.csv("data/traits/spp_list_traits_mike_withValidNames.csv")
unique(t$validName) %>% length()

tt <- t %>% 
  filter(!is.na(overwinteringStage),
         !is.na(voltinism),
         !is.na(diurnality)) %>% 
  filter(overwinteringStage != "",
         voltinism != "",
        diurnality != "")

unique(tt$validName) %>% length()
