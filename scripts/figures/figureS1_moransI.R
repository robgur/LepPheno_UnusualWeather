library(ggplot2)

load("figures/onset_correlogram_plot_revision.Rdata")
load("figures/peak_correlogram_plot_revision.Rdata")
load("figures/offset_correlogram_plot_revision.Rdata")
load("figures/dur_correlogram_plot_revision.Rdata")

p1 <- {onset_correlogram_plot + 
    scale_y_continuous(limits = c(-.005,.005)) +
    theme_classic() }

p2 <- {peak_correlogram_plot + 
    scale_y_continuous(limits = c(-.005,.005)) +
    theme_classic() }

p3 <- {offset_correlogram_plot + 
    scale_y_continuous(limits = c(-.005,.005)) +
    theme_classic() }

p4 <- {dur_correlogram_plot + 
    scale_y_continuous(limits = c(-.005,.005)) +
    theme_classic() }


cp <- cowplot::plot_grid(p1,p2,p3,p4,
                         ncol = 2, labels = c("Emergence",
                                              "Midpoint",
                                              "Termination",
                                              "Duration"),
                         label_x = 0.15)
cp


ggsave("figures/MoransI.png", width = 9, height = 4)

