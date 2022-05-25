library(ggplot2)
library(gridExtra)
library(magrittr)
library(egg)

# onset precxyear first
load("singleInteractionFigures/onset_prec_year.Rdata")

a <- ggplot(opy$data, mapping = aes(x = V1, y = mean)) +
  geom_ribbon(mapping = aes(ymin = `0.025quant`, ymax = `0.975quant`,
                            fill = V2), 
              alpha = 0.15) +
  geom_line(mapping = aes(color = V2, linetype = V2), size = 1.05) +
  labs(x = "Year", y = "Onset", 
       fill = "Precipitation", color = "Precipitation", linetype = "Precipitation") + 
  scale_color_manual(values = c("brown", "cyan", "Blue")) +
  scale_fill_manual(values = c("brown", "cyan", "Blue")) +
  theme_classic() +
  theme(legend.title = element_text(size = 12),
        legend.text = element_text(size  = 12), 
        axis.title = element_text(size = 16),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5))


# now midpoint tempSeas x yewar
load("singleInteractionFigures/mid_tempSeas_year.Rdata")

b <- ggplot(mpSy$data, mapping = aes(x = V1, y = mean)) +
  geom_ribbon(mapping = aes(ymin = `0.025quant`, ymax = `0.975quant`,
                            fill = V2), 
              alpha = 0.15) +
  geom_line(mapping = aes(color = V2, linetype = V2), size = 1.05) +
  labs(x = "Year", y = "Midpoint", 
       fill = "Temperature seasonality", color = "Temperature seasonality", 
       linetype = "Temperature seasonality") + 
  scale_color_manual(values = c("brown", "cyan", "Blue")) +
  scale_fill_manual(values = c("brown", "cyan", "Blue")) +
  theme_classic() +
  theme(legend.title = element_text(size = 12),
        legend.text = element_text(size  = 12), 
        axis.title = element_text(size = 16),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5)) 

## now get the termination plots -- tmep, tempSesa, precSeas
#tempxyear
load("singleInteractionFigures/termination_temp_year.Rdata")

d <- ggplot(tyt$data, mapping = aes(x = V1, y = mean)) +
  geom_ribbon(mapping = aes(ymin = `0.025quant`, ymax = `0.975quant`,
                            fill = V2), 
              alpha = 0.15) +
  geom_line(mapping = aes(color = V2, linetype = V2), size = 1.05) +
  labs(x = "Year", y = "Termination", 
       fill = "Temperature", color = "Temperature", 
       linetype = "Temperature") + 
  scale_color_manual(values = rev(c("firebrick", "slategray", "turquoise"))) +
  scale_fill_manual(values = rev(c("firebrick", "slategray", "turquoise"))) +
  theme_classic() +
  theme(legend.title = element_text(size = 12),
        legend.text = element_text(size  = 12), 
        axis.title = element_text(size = 16),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5)) 


# now offset tempSeas x yewar
load("singleInteractionFigures/termination_tempSeas_year.Rdata")

e <- ggplot(tytS$data, mapping = aes(x = V1, y = mean)) +
  geom_ribbon(mapping = aes(ymin = `0.025quant`, ymax = `0.975quant`,
                            fill = V2), 
              alpha = 0.15) +
  geom_line(mapping = aes(color = V2, linetype = V2), size = 1.05) +
  labs(x = "Year", y = "Termination", 
       fill = "Temperature seasonality", color = "Temperature seasonality", 
       linetype = "Temperature seasonality") + 
  scale_color_viridis_d(option = "turbo") +
  scale_fill_viridis_d(option = "turbo") +
  theme_classic() +
  theme(legend.title = element_text(size = 12),
        legend.text = element_text(size  = 12), 
        axis.title = element_text(size = 16),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5)) 


# last one is precSeas x year

load("singleInteractionFigures/termination_precSeas_year.Rdata")

f <- ggplot(mpSy$data, mapping = aes(x = V1, y = mean)) +
  geom_ribbon(mapping = aes(ymin = `0.025quant`, ymax = `0.975quant`,
                            fill = V2), 
              alpha = 0.15) +
  geom_line(mapping = aes(color = V2, linetype = V2), size = 1.05) +
  labs(x = "Year", y = "Termination", 
       fill = "Precipitation seasonality", color = "Precipitation seasonality", 
       linetype = "Precipitation seasonality") + 
  scale_color_viridis_d(option = "mako") +
  scale_fill_viridis_d(option = "mako") +
  theme_classic() +
  theme(legend.title = element_text(size = 12),
        legend.text = element_text(size  = 12), 
        axis.title = element_text(size = 16),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5)) 




p1 <- ggplot() + theme_void()


ep <- egg::ggarrange(a,b,p1,d,e,f,
                     labels = c("A", "B",
                                "", "C",
                                "D", "E"), nrow = 2, ncol = 3)
ep

ggsave(ep, 
       filename = "figures/Supplemental_YearXClimate.png",
       width = 15,  height = 8)


