library(ggplot2)
library(gridExtra)
library(magrittr)
library(egg)

#midpoint temp seas
load("singleInteractionFigures/midpoint_Seasonality_temp.Rdata")
p1 <- {mst +
    labs(color = "Season flying", fill = "Season flying") +
    theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))}
# offset temp volt
load("singleInteractionFigures/termination_Seasonality_temp.Rdata")
p2 <- {tst +
    theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))}

#onset year ows
load("singleInteractionFigures/onset_overwintering_year.Rdata")
p3 <- {ooy + 
    theme_classic() +
    theme(legend.position = "bottom")}
#duration year seas
load("singleInteractionFigures/dur_Seasonality_year.Rdata")
p4 <- {dsy + 
    theme_classic() +
    labs(color = "Season flying", fill = "Season flying") +
    theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))}


ep <- egg::ggarrange(p1, p2, p3, p4,
                     nrow = 2, ncol = 2, 
                     labels = c("A", "B", "C", "D"))


ggsave(plot = ep, filename = "figures/Fig5_interactions.png",
       dpi = 600, width = 9, height = 6)
