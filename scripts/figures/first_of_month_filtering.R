library(data.table)
library(tidyverse)
library(splitstackshape)
library(pbapply)

## read in data
leps <- fread("data/munged/LepsByGrid.csv") %>% 
  mutate(binomial = paste(genus, specificEpithet, sep = " "))

uncorrected <- ggplot(leps, mapping = aes(x = doy)) +
  geom_histogram(bins = 365) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "Day of year", y = "Number of records") +
  theme_bw()

## remove first of months
leps_noFirstOfMonths <- leps %>% 
  filter(day > 1)

corrected <- ggplot(leps_noFirstOfMonths, mapping = aes(x = doy)) +
  geom_histogram(bins = 365) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "Day of year", y = "Number of records") +
  theme_bw()



cp <- cowplot::plot_grid(uncorrected, corrected, labels = c("A", "B"),
                         nrow = 2)
cp
