library(ggplot2)
library(egg)
library(gridExtra)
library(magrittr)

#midpoint
load("singleInteractionFigures/midpoint_overwintering_year.Rdata")
load("singleInteractionFigures/mid_prec_year.Rdata")


ep <- egg::ggarrange(moy, mpy, labels = c("A", "B"))
ep

ggsave(ep, 
       filename = "figures/SupplementalMidPoint_YearEffects.png",
       width = 6, height = 6)
)