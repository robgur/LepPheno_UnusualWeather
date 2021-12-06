library(data.table)
library(tidyverse)
library(splitstackshape)
library(pbapply)

## read in data
library(tidyverse)
library(data.table)

# read in scan data
scan1 <- fread("data/SCAN1/occurrences.csv", 
               select = c("id", "institutionCode", "collectionCode", "basisOfRecord", 
                          "occurrenceID", "order", "family", "genus", "specificEpithet",
                          "scientificName", "identifiedBy", "eventDate", "year", "month", "day",
                          "recordedBy", "lifeStage", "decimalLongitude", "decimalLatitude", 
                          "locality", "county", "coordinateUncertaintyInMeters", "references"))
scan2 <- fread("data/SCAN2/occurrences.csv",
               select = c("id", "institutionCode", "collectionCode", "basisOfRecord", 
                          "occurrenceID","order", "family", "genus", "specificEpithet",
                          "scientificName", "identifiedBy", "eventDate", "year", "month", 
                          "recordedBy", "lifeStage", "decimalLongitude", "decimalLatitude", 
                          "locality", "county", "coordinateUncertaintyInMeters", "references"))
scan3 <- fread("data/SCAN3/occurrences.csv",
               select = c("id", "institutionCode", "collectionCode", "basisOfRecord", 
                          "occurrenceID","order", "family", "genus", "specificEpithet",
                          "scientificName", "identifiedBy", "eventDate", "year", "month", 
                          "recordedBy", "lifeStage", "decimalLongitude", "decimalLatitude", 
                          "locality", "county", "coordinateUncertaintyInMeters", "references"))

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
                         "recordedBy", "lifeStage", "decimalLongitude", "decimalLatitude", 
                         "locality", "county", "coordinateUncertaintyInMeters", "references"))

gbif <- gbif %>% 
  distinct(gbifID, .keep_all = T) %>% 
  mutate(source = "gbif") %>% 
  rename(id = gbifID)

# read in iDigBio data
idig <- fread("data/iDigBio/occurrence_raw.csv",
              select = c("coreid", "dwc:institutionCode", "dwc:collectionCode", "dwc:basisOfRecord", 
                         "dwc:occurrenceID", "dwc:order", "dwc:family", "dwc:genus", "dwc:specificEpithet",
                         "dwc:scientificName", "dwc:identifiedBy", "dwc:eventDate", "dwc:year", "dwc:month", 
                         "dwc:recordedBy", "dwc:lifeStage", "dwc:decimalLongitude", "dwc:decimalLatitude", 
                         "dwc:locality", "dwc:county", "dwc:coordinateUncertaintyInMeters", "dwc:associatedMedia"))


idig <- idig %>% 
  distinct(coreid, .keep_all = T) %>% 
  mutate(source = "idigbio") %>% 
  rename("references" = "dwc:associatedMedia")

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
  ggtitle("All data") +
  theme_classic()

b <- ggplot(year_sum_noDup, mapping = aes(x = year, y = count)) + 
  geom_area(mapping = aes(fill = source)) +
  scale_fill_viridis_d() +
  ggtitle("Unique records only") +
  theme_classic()

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
## and then make a bifurcated map of where the unique scan and gbif records are found

tdf <- rbind(gbif, scan, idig) 

tdf2 <- dplyr::distinct(tdf, occurrenceID, .keep_all = T)
tdf3 <- dplyr::distinct(tdf2, id, .keep_all = T)
tdf3 <- filter(tdf3, order == "Lepidoptera")

first_of_months <- tdf3 %>% 
  mutate(doy = yday(lubridate::as_date(eventDate))) %>% 
  filter(doy %in% c(1,32,60, 91,121,152,182,213,244,174305,335)) %>% 
  filter(!is.na(references)) %>% 
  filter(references != "")

samp <- slice_sample(first_of_months, n = 750)

write.csv(samp, file = "outputs/sample_of_firstOfMonths.csv", row.names = F)
