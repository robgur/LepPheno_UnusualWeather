library(tidyverse)
library(data.table)

# read in scan data
scan1 <- fread("data/SCAN1/occurrences.csv", 
               select = c("id", "institutionCode", "collectionCode", "basisOfRecord", 
                          "occurrenceID", "order", "family", "genus", "specificEpithet",
                          "scientificName", "identifiedBy", "eventDate", "year", "month", 
                          "day", "verbatimEventDate",
                          "recordedBy", "lifeStage", "decimalLongitude", "decimalLatitude", 
                          "locality", "county", "coordinateUncertaintyInMeters"))
scan2 <- fread("data/SCAN2/occurrences.csv",
               select = c("id", "institutionCode", "collectionCode", "basisOfRecord", 
                          "occurrenceID","order", "family", "genus", "specificEpithet",
                          "scientificName", "identifiedBy", "eventDate", "year", "month",
                          "day",  "verbatimEventDate",
                          "recordedBy", "lifeStage", "decimalLongitude", "decimalLatitude", 
                          "locality", "county", "coordinateUncertaintyInMeters"))
scan3 <- fread("data/SCAN3/occurrences.csv",
               select = c("id", "institutionCode", "collectionCode", "basisOfRecord", 
                          "occurrenceID","order", "family", "genus", "specificEpithet",
                          "scientificName", "identifiedBy", "eventDate", "year", "month",
                          "day",  "verbatimEventDate",
                          "recordedBy", "lifeStage", "decimalLongitude", "decimalLatitude", 
                          "locality", "county", "coordinateUncertaintyInMeters"))

#combine together and select columns
scan <- rbind(scan1, scan2, scan3) %>% 
  distinct(id, .keep_all = T)  %>% 
  mutate(source = "scan") %>% 
  mutate(id = id)

# read in gbif data
gbif <- fread("data/GBIF/occurrence.txt",
              select = c("gbifID", "institutionCode", "collectionCode", "basisOfRecord", 
                         "occurrenceID","order","family", "genus", "specificEpithet",
                         "scientificName", "identifiedBy", "eventDate", "year", "month", 
                         "day",  "verbatimEventDate",
                         "recordedBy", "lifeStage", "decimalLongitude", "decimalLatitude", 
                         "locality", "county", "coordinateUncertaintyInMeters"))

gbif <- gbif %>% 
  distinct(gbifID, .keep_all = T) %>% 
  mutate(source = "gbif") %>% 
  rename(id = gbifID)

# read in iDigBio data
idig <- fread("data/iDigBio/occurrence_raw.csv",
              select = c("coreid", "dwc:institutionCode", "dwc:collectionCode", "dwc:basisOfRecord", 
                         "dwc:occurrenceID", "dwc:order", "dwc:family", "dwc:genus", "dwc:specificEpithet",
                         "dwc:scientificName", "dwc:identifiedBy", "dwc:eventDate", "dwc:year", "dwc:month",
                         "dwc:day", "dwc:verbatimEventDate",
                         "dwc:recordedBy", "dwc:lifeStage", "dwc:decimalLongitude", "dwc:decimalLatitude", 
                         "dwc:locality", "dwc:county", "dwc:coordinateUncertaintyInMeters"))


idig <- idig %>% 
  distinct(coreid, .keep_all = T) %>% 
  mutate(source = "idigbio") 

idig <- idig %>% 
  rename_with(~gsub(pattern = "dwc:",replacement =  "", .x)) %>% 
  rename(id = coreid)

gbif$id <- as.character(gbif$id)
scan$id <- as.character(scan$id)

### combine these all together
tdf <- rbind(gbif, scan, idig)

# distinct ID's only
d  <- duplicated(tdf$id)
tdf <- mutate(tdf, dup_id = d)

#distinct occurence ID's
d <- duplicated(tdf$occurrenceID)
tdf <- mutate(tdf, dup_occID = d)

tdf$year <- as.integer(tdf$year)

year_sum <- tdf %>% 
  group_by(year, source) %>% 
  summarise(count = n()) %>% 
  filter(!is.na(year), year != "", year > 1859, year < 2022) 

tdf_noDuplicated <- tdf %>% 
  filter(dup_occID == FALSE & dup_id == FALSE)

year_sum_noDup <- tdf_noDuplicated %>% 
  group_by(year, source) %>% 
  summarise(count = n()) %>% 
  filter(!is.na(year), year != "", year > 1859, year < 2022) 


a <- ggplot(year_sum, mapping = aes(x = year, y = count)) + 
  geom_area(mapping = aes(fill = source)) +
  scale_fill_viridis_d() +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(limits = c(1850, 2020)) +
  labs(x = "Year", y = "Records") +
  ggtitle("All data") +
  theme_classic() +
  theme(legend.position="bottom") 

b <- ggplot(year_sum_noDup, mapping = aes(x = year, y = count)) + 
  geom_area(mapping = aes(fill = source)) +
  scale_fill_viridis_d() +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(limits = c(1850, 2020)) +
  labs(x = "Year", y = "Records") +
  ggtitle("Unique records only") +
  theme_classic() +
  theme(legend.position="bottom") 

## now let's make species per source for unique records only
## maybe fix this after phenometrics (a & b seem useful for now though)

spp_sum <- tdf_noDuplicated %>% 
  group_by(source) %>% 
  summarise(count = length(unique(scientificName))) 

ggplot(spp_sum, mapping = aes(x = source, y = count)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic()

cp <- cowplot::plot_grid(a,b, labels = c("A", "B"))
cp

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

panelC <- ggplot(leps, mapping = aes(x = doy)) +
  geom_histogram(bins = 365) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "Day of year", y = "Records") +
  ggtitle("Unfiltered") +
  theme_classic()

# lines were used to visualize what the data looked like for first of years
#doy1 <- filter(leps, doy == 1)
#doy32 <- filter(leps, doy == 32) 

## filter lep obs to spp_list
leps2 <- leps %>% 
  #  filter(binomial %in% spp_list) %>% 
  mutate(gcode = paste(binomial, year, id_cells, sep = ".")) %>% 
  filter(!is.na(doy)) %>% 
  filter(day > 1) 

panelD <- ggplot(leps2, mapping = aes(x = doy)) +
  geom_histogram(bins = 365) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "Day of year", y = "Records") +
  ggtitle("Filtered") +
  theme_classic()

cp2 <- cowplot::plot_grid(panelC,panelD, labels = c("C", "D"))
cp2

total_cp <- egg::ggarrange(a, b, panelC, panelD, nrow = 2, 
                           labels = c("A", "B", "C", "D"))
total_cp

ggsave(plot = total_cp, filename = "figures/data_biases.png", 
       dpi = 600, width = 10, height = 6)
