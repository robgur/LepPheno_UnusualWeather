library(dplyr)
library(sf)

#' script for more data cleaning that occurred post generating phenometrics
#' but after gathering the trait information. Year x species x cell combinations
#' that don't have enough data because they are not obligate univoltine will 
#' be removed. We will also flag cell x species (x year) combinations where 
#' species likely fly year round

## read in estimated phenometrics
pheno <- read.csv("outputs/lep_phenometrics_noNAs.csv") %>% 
  rename(id_cells = HEXcell) %>% 
  rename(species = code)

## read in traits
traits <-  read.csv("data/traits/LepSpeciesTraits_Sheet2_Synonyms.csv") %>% 
  distinct(validName, .keep_all = T)

## join traits with phenometris
pheno_traits <- left_join(pheno, traits, 
                          by = c("species" = "scientific_name"))

# minimum number of ndstcol, dstdoy, and nObs was 3,4,5 for obligate Univoltine
# for everything else it will be 3,8,10

pheno_traits_filt <- pheno_traits %>% 
  filter(voltinism == "U" |
           voltinism != "U" & dstdoy >= 8 & nObs >= 10)

ggplot(pheno_traits_filt) +
  geom_histogram(mapping = aes(x = dstdoy)) +
  scale_x_continuous(limits = c(0,20)) +
  facet_wrap(~voltinism)

ggplot(pheno_traits_filt) +
  geom_histogram(mapping = aes(x = nObs)) +
  scale_x_continuous(limits = c(0,20)) +
  facet_wrap(~voltinism)

## let's look for outliers in phenometrics to remove

outliers_onset <- ggplot(pheno_traits_filt) +
  geom_histogram(mapping = aes(x = q5), alpha = 0.5) +
  geom_point(mapping = aes(x = q5, y = 0, color = id_cells)) +
  scale_color_viridis_c() +
  theme_classic() +
  facet_wrap(~validName, ncol = 5, scales = "free")

ggsave(filename = "outputs/onset_outliers.pdf", plot = outliers_onset, 
       width = 10, height = 60,
       limitsize = FALSE)


## there is now a flagged column for species with questionable onset spreads
traits <-  read.csv("data/traits/LepSpeciesTraits_Sheet2_Synonyms.csv") %>% 
  distinct(validName, .keep_all = T)

## filter to species with okay onsets
traits_okaySpp <- traits %>% 
  filter(onsetProblems != "Y")

pheno_traits_filt2 <- pheno_traits_filt %>% 
  filter(validName %in% traits_okaySpp$validName) %>% 
  mutate(gcode = paste(validName, year, id_cells, sep = "."))

outliers_onset2 <- ggplot(pheno_traits_filt2) +
  geom_histogram(mapping = aes(x = q5), alpha = 0.5) +
  geom_point(mapping = aes(x = q5, y = 0, color = id_cells)) +
  scale_color_viridis_c() +
  theme_classic() +
  facet_wrap(~validName, ncol = 5, scales = "free")

ggsave(filename = "outputs/onset_outliers2.pdf", plot = outliers_onset2, 
       width = 10, height = 20,
       limitsize = FALSE)

## custom filter of rows to remove
tr1 <- pheno_traits_filt2 %>% 
  filter(validName == "Battus polydamas" & q5 > 200)

tr2 <- pheno_traits_filt2 %>% 
  filter(validName == "Aphrissa statira" & q5 > 175)

tr3 <- pheno_traits_filt2 %>% 
  filter(validName == "Callophrys niphon" & q5 < 100)

tr4 <- pheno_traits_filt2 %>% 
  filter(validName == "Campaea perlata" & q5 < 105)

tr5 <- pheno_traits_filt2 %>% 
  filter(validName == "Catonephele numilia" & q5 > 175)

tr6 <- pheno_traits_filt2 %>% 
  filter(validName == "Colias eurytheme" & q5 < 25)

tr7 <- pheno_traits_filt2 %>% 
  filter(validName == "Elkalyce comyntas" & id_cells < 100)

tr8 <- pheno_traits_filt2 %>% 
  filter(validName == "Hypercompe permaculata" & q5 < 110)

tr9 <- pheno_traits_filt2 %>% 
  filter(validName == "Nymphalis antiopa" & q5 > 150 & id_cells < 100)

tr10 <- pheno_traits_filt2 %>% 
  filter(validName == "Papilio glaucus" & q5 > 150)

tr11 <- pheno_traits_filt2 %>% 
  filter(validName == "Pholisora catullus" & q5 > 210)

tr12 <- pheno_traits_filt2 %>% 
  filter(validName == "Poanes zabulon" & q5 > 150 & id_cells < 100)

tr13 <- pheno_traits_filt2 %>% 
  filter(validName == "Vanessa atalanta" & q5 > 150 & id_cells < 100)

toRemove <- rbind(tr1, tr2, tr3, tr4, tr5,
                  tr6, tr7, tr8, tr9, tr10,
                  tr11, tr12, tr13) %>% 
  mutate(gcode = paste(validName, year, id_cells, sep = "."))

## remove the custom outliers
pheno_traits_filt3 <- pheno_traits_filt2 %>% 
  filter(!gcode %in% toRemove$gcode)

## now let's examine offset outliers
outliers_offset <- ggplot(pheno_traits_filt3) +
  geom_histogram(mapping = aes(x = q95), alpha = 0.5) +
  geom_point(mapping = aes(x = q95, y = 0, color = id_cells)) +
  scale_color_viridis_c() +
  theme_classic() +
  facet_wrap(~validName, ncol = 5, scales = "free")

ggsave(filename = "outputs/offset_outliers.pdf", plot = outliers_offset, 
       width = 10, height = 20,
       limitsize = FALSE)

## custom filter of rows to remove
fr1 <- pheno_traits_filt3 %>% 
  filter(validName == "Boloria chariclea" & q95 > 250)

fr2 <- pheno_traits_filt3 %>% 
  filter(validName == "Campaea perlata" & q95 > 350)

fr3 <- pheno_traits_filt3 %>% 
  filter(validName == "Chlosyne rosita")

fr4 <- pheno_traits_filt3 %>% 
  filter(validName == "Colias eurytheme" & q95 < 200)

fr5 <- pheno_traits_filt3 %>% 
  filter(validName == "Cyclogramma bacchis" & q95 < 280)

fr6 <- pheno_traits_filt3 %>% 
  filter(validName == "Euptoieta claudia" & q95 < 200)

fr7 <- pheno_traits_filt3 %>% 
  filter(validName == "Hamadryas atlantis" & q95 < 200)

fr8 <- pheno_traits_filt3 %>% 
  filter(validName == "Hesperocharis costaricensis" & q95 < 280)

fr9 <- pheno_traits_filt3 %>% 
  filter(validName == "Junonia coenia" & q95 < 200 & id_cells < 100)

fr10 <- pheno_traits_filt3 %>% 
  filter(validName == "Limenitis archippus" & q95 <150)

fr11 <- pheno_traits_filt3 %>% 
  filter(validName == "Phyciodes pallescens" & q95 <150)

fr <- rbind(fr1, fr2, fr3, fr4, fr5,
            fr6, fr7, fr8, fr9, fr10,
            fr11)

pheno_traits_filt4 <- pheno_traits_filt3 %>% 
  filter(!gcode %in% fr$gcode)

# are there species with less than 3 phenometrics? if so remove
pheno_sum <- pheno_traits_filt4 %>% 
  group_by(validName) %>% 
  summarise(count = n()) # we all good three is the lowest

# are there cells with less than 3 phenometrics?
pheno_sum <-  pheno_traits_filt4 %>% 
  group_by(id_cells) %>% 
  summarise(count = n()) # yes many

tooFewData <- filter(pheno_sum, count < 3)

pheno_traits_filt5 <- pheno_traits_filt4 %>% 
  filter(!id_cells %in% tooFewData$id_cells)

write.csv(x = pheno_traits_filt5, 
          file = "outputs/filteredPhenometricsTraits.csv",
          row.names = F)
