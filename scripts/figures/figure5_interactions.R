library(ggplot2)
library(gridExtra)
library(magrittr)
library(egg)
library(dplyr)

#midpoint temp seas 
load("singleInteractionFigures/midpoint_Seasonality_temp.Rdata")

mstd <- mst$data %>% 
  rename(lowci = 3, highci = 5)

p3 <- ggplot(mst$data, mapping = aes(x = V1, y = mean)) +
  geom_ribbon(mapping = aes(ymin = `0.025quant`, ymax = `0.975quant`,
                            fill = V2), 
              alpha = 0.15) +
  geom_line(mapping = aes(color = V2, linetype = V2), size = 1.05) +
  labs(x = "Temperature", y = "Midpoint", 
       fill = "Season flying", 
       color = "Season flying",
       linetype = "Season flying") + 
  scale_fill_viridis_d(option = "inferno") +
  scale_color_viridis_d(option = "inferno") +
  theme_classic() +
  theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))

# duration temp voltinism
load("singleInteractionFigures/termination_voltinism_temp.Rdata")
p4 <- ggplot(tvt$data, mapping = aes(x = V1, y = mean)) +
  geom_ribbon(mapping = aes(ymin = `0.025quant`, ymax = `0.975quant`,
                            fill = V2), 
              alpha = 0.15) +
  geom_line(mapping = aes(color = V2, linetype = V2), size = 1.05) +
  labs(x = "Temperature", y = "Termination", 
       fill = "Voltinism", color = "Voltinism", linetype = "Voltinism") + 
  scale_fill_viridis_d(option = "rocket", begin = 0.3, end = 0.7) +
  scale_color_viridis_d(option = "rocket", begin = 0.3, end = 0.7) +
  ggtitle("") +
  theme_classic() +
  theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))

#onset year ows
load("singleInteractionFigures/onset_overwintering_year.Rdata")
p1 <-  ggplot(ooy$data, mapping = aes(x = V1, y = mean)) +
  geom_ribbon(mapping = aes(ymin = `0.025quant`, ymax = `0.975quant`,
                            fill = V2), 
              alpha = 0.15) +
  geom_line(mapping = aes(color = V2, linetype = V2), size = 1.05) +
  labs(x = "Year", y = "Onset",
       fill = "Diapause stage", color = "Diapause stage", linetype = "Diapause stage") + 
  scale_fill_viridis_d() +
  scale_color_viridis_d() +
  theme_classic() +
  theme(legend.position = "bottom")
#offset year ows
load("singleInteractionFigures/termination_overwintering_year.Rdata")
p2 <- ggplot(toy$data, mapping = aes(x = V1, y = mean)) +
  geom_ribbon(mapping = aes(ymin = `0.025quant`, ymax = `0.975quant`,
                            fill = V2), 
              alpha = 0.15) +
  geom_line(mapping = aes(color = V2, linetype = V2), size = 1.05) +
  labs(x = "Year", y = "Termination", 
       fill = "Diapause stage", color = "Diapause stage", linetype = "Diapause stage") + 
  scale_fill_viridis_d() +
  scale_color_viridis_d() +
  ggtitle(label = "") + 
  theme_classic() +
  theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))


ep <- egg::ggarrange( p3, p4, p1, p2,
                      nrow = 2, ncol = 2, 
                      labels = c("A", "B", "C", "D"))


ggsave(plot = ep, filename = "figures/Fig5_interactions.png",
       dpi = 600, width = 11, height = 7)
