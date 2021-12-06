library(ggplot2)
library(gridExtra)
library(magrittr)

#onset 
load("singleInteractionFigures/onset_prec_temp.Rdata")
#midpoint
load("singleInteractionFigures/mid_prec_temp.Rdata")
#duration
load("singleInteractionFigures/dur_prec_temp.Rdata")

####
ep <- egg::ggarrange(opt, mpt, dpt, labels = c("A", "B", "C"),
                     ncol = 2)

ggsave(filename = "figures/supplemental_temp_prec_interaction.png",
       plot = ep, dpi = 600,
       width = 8, height = 4)
