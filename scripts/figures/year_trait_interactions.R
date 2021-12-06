library(ggplot2)
library(gridExtra)
library(magrittr)
library(egg)

#onset 
load("singleInteractionFigures/onset_overwintering_year.Rdata")
load("singleInteractionFigures/onset_prec_year.Rdata")
#midpoint
load("singleInteractionFigures/mid_precSeas_year.Rdata")
#duration
load("singleInteractionFigures/dur_Seasonality_year.Rdata")


p1 <- {ooy + 
    theme_classic() +
    theme(legend.position = "bottom")}

p2 <- {opy + 
    theme_classic() +
    theme(legend.position = "bottom")}


p3 <- {mpSy + 
    scale_color_viridis_d(option = "mako", direction = -1) +
    scale_fill_viridis_d(option = "mako", direction = -1) + 
    theme_classic() +
    theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))}

p4 <- {dsy + 
    theme_classic() +
    theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))}


ep <- egg::ggarrange(p1, p2, p3, p4,
                     nrow = 2, ncol = 2, 
                     labels = c("A", "B", "C", "D"))


ggsave(plot = ep, filename = "figures/year_trait_interactions.png",
       dpi = 600, width = 8, height = 6)

