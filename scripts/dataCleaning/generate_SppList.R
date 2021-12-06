library(data.table)
library(tidyverse)
library(splitstackshape)
library(pbapply)

## read in data
leps <- fread("data/munged/LepsByGrid.csv") %>% 
  mutate(binomial = paste(genus, specificEpithet, sep = " "))

tropics <- read.csv("data/CellsWithKoeppenTropics.csv")

leps <- leps %>% 
  filter(!id_cells %in% tropics$id_cells)

ggplot(leps, mapping = aes(x = doy)) +
  geom_histogram(bins = 365) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw()

# lines were used to visualize what the data looked like for first of years
#doy1 <- filter(leps, doy == 1)
#doy32 <- filter(leps, doy == 32) 

#spp_list <- read.csv("spp_list.csv")
#spp_list <- spp_list$scientific_name

## filter lep obs to spp_list
leps <- leps %>% 
  #  filter(binomial %in% spp_list) %>% 
  mutate(gcode = paste(binomial, year, id_cells, sep = ".")) %>% 
  filter(!is.na(doy)) %>% 
  filter(doy > 1) # maybe this is not needed

leps_sum <- leps %>% 
  group_by(year, binomial, id_cells) %>% 
  summarize(ndstcol = n_distinct(recordedBy), 
            dstdoy=n_distinct(doy), 
            nObs=n())

enough_data <- leps_sum %>% 
  filter(ndstcol >= 3, dstdoy >= 4, nObs >= 5) %>% 
  mutate(gcode = paste(binomial, year, id_cells, sep = ".")) %>% 
  mutate(spaces = str_count(binomial, pattern = " ")) %>% 
  filter(!is.na(id_cells))

spp_list <- data.frame(spp = unique(enough_data$binomial))
spp_list <- spp_list %>% 
  mutate(spaces = str_count(string = spp, pattern ="\\S+")) %>% 
  filter(spaces >= 2)
length(spp_list$spp)

enough_data <- enough_data %>% 
  filter(binomial %in% spp_list$spp)

# remove spp without at least 3 things
passing_spp <- enough_data %>% 
  group_by(binomial) %>% 
  summarise(count = n()) %>% 
  filter(count >= 5)

enough_data <- enough_data %>% 
  filter(binomial %in% passing_spp$binomial)

passing_cells <- enough_data %>% 
  group_by(id_cells) %>% 
  summarize(count = n()) %>% 
  filter(count >= 5)

enough_data <- enough_data %>% 
  filter(id_cells %in% passing_cells$id_cells)

length(unique(enough_data$binomial))

## need to be in 5 cells and 5 years
spp_list_df <- data.frame(scientific_name = unique(enough_data$binomial))

trait_df <- read.csv("data/traits/LepSpeciesTraits_Sheet2_Synonyms.csv")
traits_cross <- read.csv("TraitsButterflyCrossley.csv")

spp_list_traits <- left_join(spp_list_df, trait_df)
spp_list_traits <- left_join(spp_list_traits, traits_cross, 
                             by = c("scientific_name" = "Species"))

write.csv(spp_list_traits, file = "data/traits/spp_list_traits_mike.csv", row.names = F)

## given spp list, generate synonyms
library(taxotools)
spp_list <- read.csv('data/traits/spp_list_traits_mike.csv')
spp_list <- unique(spp_list$scientific_name)
wiki_syn <- list_wiki_syn(spp_list)

wiki_syn_filt <- filter(wiki_syn, Syn != "")

## join to Mike's List
mike_traits <- read.csv("data/traits/spp_list_traits_mike.csv")
mike_traits$scientific_name <- as.character(mike_traits$scientific_name)
wiki_syn_filt$Name <- as.character(wiki_syn_filt$Name)
wiki_syn_filt$WikiName <- as.character(wiki_syn_filt$WikiName)
wiki_syn_filt$Syn <- as.character(wiki_syn_filt$Syn)
wiki_syn_filt$OrigSyn <- as.character(wiki_syn_filt$OrigSyn)

mike_traits <- left_join(mike_traits, wiki_syn_filt,
        by = c("scientific_name" = "Name")) %>% 
  mutate(validName = WikiName)

mike_traits2 <- mike_traits %>% 
    mutate(validName = if_else(condition = is.na(validName), 
                                true = scientific_name, 
                                false = validName))

write.csv(mike_traits2, file = "data/traits/spp_list_traits_mike_withValidNames.csv", row.names = F)
