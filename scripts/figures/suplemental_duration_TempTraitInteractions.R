library(ggplot2)
library(egg)

#duration
load("singleInteractionFigures/dur_overwintering_temp.Rdata")
load("singleInteractionFigures/dur_Seasonality_temp.Rdata")
load("singleInteractionFigures/dur_voltinism_temp.Rdata")

dot <- dot + labs(title = "") + theme(legend.position = "right")
dst <- dst + theme(legend.position = "right")
dvt <- dvt + theme(legend.position = "right")

egg_plot <- egg::ggarrange(dot, dst, dvt, nrow = 2, ncol = 2, 
                           labels = c("A", "B", "C"))


ggsave(egg_plot, 
       filename = "figures/SupplementalDuration_TempTraitInteractions.png",
       width = 8, height = 4)

