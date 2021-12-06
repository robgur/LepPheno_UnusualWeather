library(tidyverse)
library(sf)
library(taxotools)
library(dbplyr)
library(RSQLite)

## read in data
mdf <- read.csv("data/LMM_Data/mdf_pruned_hopkins.csv") %>% 
  dplyr::filter(overwinteringStage != "None")


## get family names
#fams <- taxotools::list_higher_taxo(indf = mdf, canonical = "validName") # done once
myConn <- dbConnect(drv = SQLite(), dbname= "taxo.db")
dbListTables(myConn)
fams <- dbReadTable(myConn, "taxo") %>% 
  collect()


fams_sum <- fams %>% 
  group_by(Family) %>% 
  summarise(Species = n()) %>% 
  na.omit()

fams_sum_10 <- dplyr::top_n(fams_sum, 10, wt = Species)
fams_sum_10 <- fams_sum_10 %>% 
  mutate(Family = fct_reorder(Family, Species, .desc = F))

ggplot() +
  geom_bar(fams_sum_10, mapping = aes(y = Family, x = Species), stat = "identity") +
  scale_x_continuous(expand = c(0,0)) + 
  labs(x = "Number of Species", y = "") +
  theme_classic()

ggsave("figures/speciesPerFamily.png", dpi = 600, width = 4, height = 3)
