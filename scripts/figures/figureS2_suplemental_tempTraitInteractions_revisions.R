library(ggplot2)
library(gridExtra)
library(magrittr)
library(egg)

#OWS 
load("singleInteractionFigures/onset_overwintering_temp.Rdata")

a <- ggplot(oot$data, mapping = aes(x = V1, y = mean)) +
  geom_ribbon(mapping = aes(ymin = `0.025quant`, ymax = `0.975quant`,
                            fill = V2), 
              alpha = 0.15) +
  geom_line(mapping = aes(color = V2, linetype = V2), size = 1.05) +
  labs(x = "Temperature", y = "Onset", 
       fill = "Diapause stage", color = "Diapause stage",
       linetype = "Diapause stage") + 
  scale_fill_viridis_d() +
  scale_color_viridis_d() +
  theme_classic() +
  theme(legend.title = element_text(size = 12),
        legend.text = element_text(size  = 12), 
        axis.title = element_text(size = 16),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5)) 

#midpoint ows
load("singleInteractionFigures/midpoint_overwintering_temp.Rdata")

b <- ggplot(mot$data, mapping = aes(x = V1, y = mean)) +
  geom_ribbon(mapping = aes(ymin = `0.025quant`, ymax = `0.975quant`,
                            fill = V2), 
              alpha = 0.15) +
  geom_line(mapping = aes(color = V2, linetype = V2), size = 1.05) +
  labs(x = "Temperature", y = "Midpoint", 
       fill = "Diapause stage", color = "Diapause stage", linetype = "Diapause stage") + 
  scale_fill_viridis_d() +
  scale_color_viridis_d() +
  theme_classic() + 
  theme(legend.title = element_text(size = 12),
                          legend.text = element_text(size  = 12), 
                          axis.title = element_text(size = 16),
                          legend.position = "bottom",
                          plot.title = element_text(hjust = 0.5)) 


#offset ows
load("singleInteractionFigures/termination_overwintering_temp.Rdata")

c <- ggplot(tot$data, mapping = aes(x = V1, y = mean)) +
  geom_ribbon(mapping = aes(ymin = `0.025quant`, ymax = `0.975quant`,
                            fill = V2), 
              alpha = 0.15) +
  geom_line(mapping = aes(color = V2, linetype = V2), size = 1.05) +
  labs(x = "Temperature", y = "Termination", 
       fill = "Diapause stage", color = "Diapause stage", linetype = "Diapause stage") + 
  scale_fill_viridis_d() +
  scale_color_viridis_d() +
  theme_classic() +
  theme(legend.title = element_text(size = 12),
        legend.text = element_text(size  = 12), 
        axis.title = element_text(size = 16),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5)) 


# now seas, seas and then volt
# offset
load("singleInteractionFigures/termination_Seasonality_temp.Rdata")

d <- ggplot(tst$data, mapping = aes(x = V1, y = mean)) +
  geom_ribbon(mapping = aes(ymin = `0.025quant`, ymax = `0.975quant`,
                            fill = V2), 
              alpha = 0.15) +
  geom_line(mapping = aes(color = V2, linetype = V2), size = 1.05) +
  labs(x = "Temperature", y = "Termination", 
       fill = "Season flying", color = "Season flying", linetype = "Season flying") + 
  scale_fill_viridis_d(option = "inferno") +
  scale_color_viridis_d(option = "inferno") +
  theme_classic() +
  theme(legend.title = element_text(size = 12),
        legend.text = element_text(size  = 12), 
        axis.title = element_text(size = 16),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5)) 


# now do seas and duration
load("singleInteractionFigures/dur_Seasonality_temp.Rdata")

e <- ggplot(dst$data, mapping = aes(x = V1, y = mean)) +
  geom_ribbon(mapping = aes(ymin = `0.025quant`, ymax = `0.975quant`,
                            fill = V2), 
              alpha = 0.15) +
  geom_line(mapping = aes(color = V2, linetype = V2), size = 1.05) +
  labs(x = "Temperature", y = "Duration", 
       fill = "Season flying", color = "Season flying", linetype = "Season flying") + 
  scale_fill_viridis_d(option = "inferno") +
  scale_color_viridis_d(option = "inferno") +
  theme_classic() +
  theme(legend.title = element_text(size = 12),
        legend.text = element_text(size  = 12), 
        axis.title = element_text(size = 16),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5)) 


# last one offset temp*voltinism
load("singleInteractionFigures/dur_voltinism_temp.Rdata")

f <- ggplot(dvt$data, mapping = aes(x = V1, y = mean)) +
  geom_ribbon(mapping = aes(ymin = `0.025quant`, ymax = `0.975quant`,
                            fill = V2), 
              alpha = 0.15) +
  geom_line(mapping = aes(color = V2, linetype = V2), size = 1.05) +
  labs(x = "Temperature", y = "Duration", fill = "Voltinism", 
       color = "Voltinism", linetype = "Voltinism") + 
  scale_fill_viridis_d(option = "rocket", begin = 0.3, end = 0.7) +
  scale_color_viridis_d(option = "rocket", begin = 0.3, end = 0.7) +
  theme_classic() +
  theme(legend.title = element_text(size = 12),
        legend.text = element_text(size  = 12), 
        axis.title = element_text(size = 16),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5)) 



ep <- egg::ggarrange(a,b,c,d,e,f,
                     labels = c("A", "B",
                                "C", "D",
                                "E", "F"), nrow = 2, ncol = 3)
ep

ggsave(ep, 
       filename = "figures/Supplemental_TempXTraits.png",
       width = 19,  height = 8)
