library(tidyverse)

## read in data
mdf <- read.csv("data/LMM_Data/mdf_hopkins_onset.csv") %>% 
  dplyr::filter(overwinteringStage != "None") %>% 
  mutate(voltinism = if_else(condition = voltinism == "F",
                             true = "M",
                             false = voltinism))
mdf <- mdf %>% 
  mutate(volt2 = if_else(voltinism == "M",
                         true = "Multivoltine",
                         false = "Univoltine"))


spp_sd <- mdf %>% 
  group_by(validName) %>% 
  summarise(onset_low_sd = mean(hopkins_onset, na.rm = T) - 2*sd(hopkins_onset, na.rm = T),
            offset_low_sd = mean(hopkins_offset, na.rm = T) - 2*sd(hopkins_offset, na.rm = T),
            onset_high_sd = mean(hopkins_onset, na.rm = T) + 2*sd(hopkins_onset, na.rm = T),
            offset_high_sd = mean(hopkins_offset, na.rm = T) + 2*sd(hopkins_offset, na.rm = T))

mdf_sd <- left_join(mdf, spp_sd) 

## remove outliers
mdf_sd_outlier <- mdf_sd %>% 
  dplyr::mutate(Outlier = if_else(condition = 
              onset_low_sd < hopkins_onset & hopkins_onset < onset_high_sd &
              offset_low_sd < hopkins_offset & hopkins_offset < offset_high_sd,
              true = "No",
              false = "Yes")) %>% 
  dplyr::filter(!is.na(Outlier))


ggplot(mdf_sd_outlier) + 
  geom_hline(yintercept = 8, linetype = "dashed") +
  geom_jitter(aes(x=year, y=dstdoy, 
                  size = Outlier, color = Outlier), alpha = 0.5) + 
  scale_size_manual(values = c(1, 4)) +
  scale_color_manual(values = c("dodgerblue3", "gold")) +
  labs(x = "Year", y = "Distinct observation days") +
  facet_wrap(~volt2) +
  theme_classic()

ggsave("figures/Outlier_Phenometrics.png", dpi = 600,
       width = 8, height = 4)
