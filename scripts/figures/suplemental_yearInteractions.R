library(ggplot2)
library(gridExtra)
library(magrittr)
library(egg)

#onset 
load("singleInteractionFigures/onset_prec_year.Rdata")
#midpoint
load("singleInteractionFigures/mid_precSeas_year.Rdata")
#midpoint
load("singleInteractionFigures/midpoint_overwintering_year.Rdata")
load("singleInteractionFigures/mid_prec_year.Rdata")

p1 <-  {opy + 
    theme_classic() +
    theme(legend.position = "bottom")}
p2 <- {moy + 
    theme_classic() +
    theme(legend.position = "bottom")}
p3 <-  {mpy + 
    theme_classic() +
    theme(legend.position = "bottom")}
p4 <- {mpSy + 
    scale_color_viridis_d(option = "mako", direction = -1) +
    scale_fill_viridis_d(option = "mako", direction = -1) + 
    theme_classic() +
    theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))}




ep <- egg::ggarrange(p1, p2, 
                     p3, p4,
                     labels = c("A", "B",
                                "C", "D"))
ep

ggsave(ep, 
       filename = "figures/SupplementalYearEffects.png",
       width = 8, height = 6)
