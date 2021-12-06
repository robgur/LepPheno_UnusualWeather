library(tidyverse)

## read in data
mdf <- read.csv("data/LMM_Data/mdf_hopkins_onset.csv") %>% 
  dplyr::filter(overwinteringStage != "None")

spp_sd <- mdf %>% 
  group_by(validName) %>% 
  summarise(onset_low_sd = mean(hopkins_onset, na.rm = T) - 2*sd(hopkins_onset, na.rm = T),
         offset_low_sd = mean(hopkins_offset, na.rm = T) - 2*sd(hopkins_offset, na.rm = T),
         onset_high_sd = mean(hopkins_onset, na.rm = T) + 2*sd(hopkins_onset, na.rm = T),
         offset_high_sd = mean(hopkins_offset, na.rm = T) + 2*sd(hopkins_offset, na.rm = T))

mdf_sd <- left_join(mdf, spp_sd)

outliers_onset <- ggplot(mdf_sd) +
  geom_histogram(mapping = aes(x = hopkins_onset), alpha = 0.5) +
  geom_point(mapping = aes(x = hopkins_onset, y = 0, color = id_cells)) +
  geom_vline(mapping = aes(xintercept = onset_low_sd)) +
  geom_vline(mapping = aes(xintercept = onset_high_sd)) +
  scale_color_viridis_c() +
  theme_classic() +
  facet_wrap(~validName, ncol = 5, scales = "free")

ggsave(filename = "outputs/onset_outliers.pdf", plot = outliers_onset, 
       width = 10, height = 60,
       limitsize = FALSE)

outliers_offset <- ggplot(mdf_sd) +
  geom_histogram(mapping = aes(x = hopkins_offset), alpha = 0.5) +
  geom_point(mapping = aes(x = hopkins_offset, y = 0, color = id_cells)) +
  geom_vline(mapping = aes(xintercept = offset_low_sd)) +
  geom_vline(mapping = aes(xintercept = offset_high_sd)) +
  scale_color_viridis_c() +
  theme_classic() +
  facet_wrap(~validName, ncol = 5, scales = "free")

ggsave(filename = "outputs/offset_outliers.pdf", plot = outliers_offset, 
       width = 10, height = 60,
       limitsize = FALSE)

## remove outliers
mdf_sd_pruned <- mdf_sd %>% 
  dplyr::filter(onset_low_sd < hopkins_onset & hopkins_onset < onset_high_sd &
                offset_low_sd < hopkins_offset & hopkins_offset < offset_high_sd)

## remove species and cells without at least three estimates
spp_sum <- mdf_sd_pruned %>% 
  group_by(validName) %>% 
  summarise(count = n())

cell_sum <- mdf_sd_pruned %>% 
  group_by(id_cells) %>% 
  summarise(count = n())

enough_cells <- dplyr::filter(cell_sum, count >= 2)

mdf_sd_pruned2 <- mdf_sd_pruned %>% 
  dplyr::filter(id_cells %in% enough_cells$id_cells)

spp_sum <- mdf_sd_pruned2 %>% 
  group_by(validName) %>% 
  summarise(count = n())

enough_spp <- dplyr::filter(spp_sum, count >= 2)

mdf_sd_pruned3 <- mdf_sd_pruned2 %>% 
  dplyr::filter(validName %in% enough_spp$validName)

## now make sure we're good on cells still
cell_sum <- mdf_sd_pruned3 %>% 
  group_by(id_cells) %>% 
  summarise(count = n())

# we good. save the mdf after making seasonality trait
spp_ave_onset <- mdf_sd_pruned3 %>% 
  group_by(validName) %>% 
  summarize(meanOnset = mean(hopkins_onset))

## Can't use equinox to group species into traits. let's try quantiles
quantile(spp_ave_onset$meanOnset, c(0.15,0.5, 0.85))

ggplot() + 
  geom_histogram(spp_ave_onset, mapping = aes(x = meanOnset), fill = "turquoise") +
  geom_vline(xintercept = c(50, 123.5)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() 



mdf_sd_pruned3 <- left_join(mdf_sd_pruned3, spp_ave_onset) %>% 
  mutate(Seas = case_when(
    meanOnset < 49.96 ~ "Spring",
    meanOnset >= 49.96 & meanOnset <= 123.45 ~ "Summer",
    meanOnset > 123.45 ~ "Fall"
  ))

write.csv(x = mdf_sd_pruned3,
          file = "data/LMM_Data/mdf_pruned_hopkins.csv",
          row.names = F)
