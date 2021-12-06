library(ggplot2)
library(gridExtra)
library(magrittr)

#onset 
load("singleInteractionFigures/onset_overwintering_temp.Rdata")
#midpoint
load("singleInteractionFigures/midpoint_overwintering_temp.Rdata")
load("singleInteractionFigures/midpoint_Seasonality_temp.Rdata")
#offset
load("singleInteractionFigures/termination_overwintering_temp.Rdata")
load("singleInteractionFigures/termination_Seasonality_temp.Rdata")
load("singleInteractionFigures/termination_voltinism_temp.Rdata")


v <- {tvt + 
    theme(legend.position = "left", legend.title = element_blank(),
          text = element_text(size=14)) }%>% 
  ggpubr::get_legend() %>% arrangeGrob(left = "Voltinism") 

s <- {mst + 
    theme(legend.position = "left", legend.title = element_blank(),
          text = element_text(size=14)) }%>% 
  ggpubr::get_legend() %>% arrangeGrob(left = "Seasonality") 

d <- {oot + 
    theme(legend.position = "left", legend.title = element_blank(),
          text = element_text(size=14)) }%>% 
  ggpubr::get_legend() %>% arrangeGrob(left = "Diapause stage") 


oot <- oot + labs(title = expression(paste(underline("Onset")))) + 
  theme(legend.position = "none", text = element_text(size=20)) 
oot <- ggpubr::ggarrange(oot, labels = "A",
                         font.label = list(size = 24, face = "bold.italic"))
mot <- mot + labs(title = expression(paste(underline("Midpoint")))) +
  theme(text = element_text(size=20))
mot <- ggpubr::ggarrange(mot, labels = "B",
                         font.label = list(size = 24, face = "bold.italic"))
tot <- tot + labs(title = expression(paste(underline("Termination")))) +
  theme(text = element_text(size=20))
tot <- ggpubr::ggarrange(tot, labels = "C",
                         font.label = list(size = 24, face = "bold.italic"))

mst <- mst + theme(legend.position = "none", text = element_text(size=20))
mst <- ggpubr::ggarrange(mst, labels = "D",
                         font.label = list(size = 24, face = "bold.italic"))
tst <- tst + theme(text = element_text(size=20))
tst <- ggpubr::ggarrange(tst, labels = "E",
                         font.label = list(size = 24, face = "bold.italic"))

tvt <- tvt + theme(legend.position = "none", text = element_text(size=20))
tvt <- ggpubr::ggarrange(tvt, labels = "F", 
                         font.label = list(size = 24, face = "bold.italic"))

ga <- grid.arrange(d, oot, mot, tot,
                         s,mst ,tst,
                         v,tvt,
                         layout_matrix = 
                           matrix(data = c(
                             NA,1,rep(2,3),rep(3,3),rep(4,3),
                             rep(NA, 4),6,rep(7,3),rep(8,3),
                             rep(NA,7),10,rep(11,3)), 
                             nrow = 3, ncol = 11, byrow = T))

ggsave(plot = ga, filename = "figures/temperature_trait_interactions.png",
       dpi = 600, width = 15, height = 9)
