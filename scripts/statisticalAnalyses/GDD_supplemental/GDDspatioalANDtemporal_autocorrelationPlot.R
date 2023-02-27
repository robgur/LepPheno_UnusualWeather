library(cowplot)
library(ggplot2)

# temporal autocorrelation 

load("Figures/GDDonset_acf.rds")
load("Figures/GDDoffset_acf.rds")
load("Figures/GDDdur_acf.rds")

cowplot::plot_grid(on_acf, off_acf, dur_acf, ncol = 1, 
                   labels = c("Onset", "Termination", "Duration"))

ggsave("ManuscriptFigures/SI3_Fig4.png", height = 6, width = 4)

#spatial autocorrelation

load("Figures/GDDonset_correlogram.rds")
load("Figures/GDDoffset_correlogram.rds")
load("Figures/GDDdur_correlogram.rds")

cowplot::plot_grid(onset_correlogram_plot, 
                   offset_correlogram_plot, 
                   duration_correlogram_plot, ncol = 1, 
                   labels = c("Onset", "Termination", "Duration"),
                   label_x = 0.15)

ggsave("ManuscriptFigures/SI3_Fig5.png", height = 8, width = 5)
