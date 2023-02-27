library(phyr)
library(car)
library(ggplot2)
library(ncf)
library(fpp2)

# get data ready for pglmm models
source('scripts/pglmm_dataSetUp_GDDModels.R')

#model_onset_with_anomaly_reduced_vif_ok
onset_m <- pglmm(onset ~ accum_GDD_sc + UnusualColdSum + UnusualWarmSum + 
                          voltinism + overwinteringStage
                        + (1 | unique_name__) + (1 | id_cells) + 
                          UnusualColdSum:UnusualWarmSum,
                        data = mdf_phylo,
                        cov_ranef = list(unique_name = tt),
                        bayes = T)


# WAIC(onset_m_reduce)
# 4621
#WAIC(onset_m_noAnom_reduced)
#4673

plot_bayes(onset_m)
#ggsave(filename = "Figures/onset_posteriorDistributions_GDDModels.pdf",
#       width = 8, height = 9)

summary(onset_m)
rr2::R2(onset_m) #0.7985

#save model results as table
im <- onset_m$inla.model
im_fix <- im$summary.fixed %>% 
  tibble::rownames_to_column("Effects") %>% 
  dplyr::select(Effects, mean, `0.025quant`, `0.975quant`)
sjPlot::tab_df(im_fix, title = "PGLMM Onset", file = "Tables/onset_GDDPGLMM.doc")
write.csv(im_fix, file = "Tables/onsetResults_GDDPGLMM.csv", row.names = F)

## Model Assumption Checks
# residuals
resids <- residuals(onset_m)
hist(resids)

resid_onset <- mdf_phylo %>% 
  mutate(resid = resids)

## look into temporal autocorrelation
resid_onset_sort <- resid_onset[order(resid_onset$year),]

acfPlot <- acf(resid_onset_sort$resid) 

on_acf <- ggAcf(resid_onset_sort$resid) +
  theme_bw() +
  ggtitle("")

save(on_acf, file = "Figures/GDDonset_acf.rds")

ggplot(resid_onset, mapping = aes(x = year, y = resid)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(x = "Year", y = "Residual")+
  theme_bw()

## Examine spatial autocorrelation a second way
rdf <- mutate(mdf_phylo, residuals = resids)

fit = correlog(x = rdf$lon, y = rdf$lat, z = rdf$residuals,
               increment = 250000, latlon = F, resamp = 100)
plot(fit)

fit_df <- data.frame(distance = fit$mean.of.class, 
                     correlation = fit$correlation,
                     p.value = if_else(fit$p < 0.025, 
                                       true = "Sig",
                                       false = "Not Sig")) %>% 
  filter(distance < 5000000 & distance > 1)


onset_correlogram_plot <- ggplot() +
  geom_line(fit_df,
            mapping = aes(x = distance, y = correlation)) +
  geom_point(fit_df,
             mapping = aes(x = distance, y = correlation, fill = p.value),
             size = 3, shape = 21, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("white", "black")) +
  labs(x = "Distance", y = "Moran's I", fill = "P-Value") +
  theme_classic()

onset_correlogram_plot # no evidence of significant spatial autocorrelation

save(onset_correlogram_plot, file = "Figures/GDDonset_correlogram.rds")

################################################################################
################################################################################
#############Interaction Figures ###############################################
################################################################################
################################################################################
# make interaction figures rob included in slack
inla_mdf <- mdf_phylo %>% 
  dplyr::select(onset, 
                accum_GDD_sc, 
                UnusualColdSum, 
                UnusualWarmSum, 
                voltinism, 
                overwinteringStage, 
                id_cells,
                unique_name)

pred_vals <- -3:3

#UnusualColdDays is x axis (so pred_vals),
#UnusualWarmDays is legend, so start with -1
term1 <- data.frame(onset = rep(NA, length(pred_vals)),
                      accum_GDD_sc = rep(NA, length(pred_vals)), 
                      UnusualColdSum = pred_vals, 
                      UnusualWarmSum = -1,
                      voltinism = rep(NA, length(pred_vals)), 
                      overwinteringStage = rep(NA, length(pred_vals)), 
                      id_cells = rep(NA, length(pred_vals)), 
                      unique_name = rep(NA, length(pred_vals)))

#UnusualColdDays is x axis (so pred_vals),
#UnusualWarmDays is legend, so start with 0
term2 <- data.frame(onset = rep(NA, length(pred_vals)),
                    accum_GDD_sc = rep(NA, length(pred_vals)), 
                    UnusualColdSum = pred_vals, 
                    UnusualWarmSum = 0,
                    voltinism = rep(NA, length(pred_vals)), 
                    overwinteringStage = rep(NA, length(pred_vals)), 
                    id_cells = rep(NA, length(pred_vals)), 
                    unique_name = rep(NA, length(pred_vals)))

#UnusualColdDays is x axis (so pred_vals),
#UnusualWarmDays is legend, so start with 0
term3 <- data.frame(onset = rep(NA, length(pred_vals)),
                    accum_GDD_sc = rep(NA, length(pred_vals)), 
                    UnusualColdSum = pred_vals, 
                    UnusualWarmSum = 1,
                    voltinism = rep(NA, length(pred_vals)), 
                    overwinteringStage = rep(NA, length(pred_vals)), 
                    id_cells = rep(NA, length(pred_vals)), 
                    unique_name = rep(NA, length(pred_vals)))

pred.df.1 <- rbind(inla_mdf, term1)
pred.df.2 <- rbind(inla_mdf, term2)
pred.df.3 <- rbind(inla_mdf, term3)

# fit models
pglmm1 <-  pglmm(onset ~ accum_GDD_sc + UnusualColdSum + UnusualWarmSum + 
                   voltinism + overwinteringStage
                 + (1 | unique_name__) + (1 | id_cells) + 
                   UnusualColdSum:UnusualWarmSum,
                 data = pred.df.1,
                 cov_ranef = list(unique_name = tt),
                 bayes = T)

pglmm2 <-  pglmm(onset ~ accum_GDD_sc + UnusualColdSum + UnusualWarmSum + 
                   voltinism + overwinteringStage
                 + (1 | unique_name__) + (1 | id_cells) + 
                   UnusualColdSum:UnusualWarmSum,
                 data = pred.df.2,
                 cov_ranef = list(unique_name = tt),
                 bayes = T)

pglmm3 <-  pglmm(onset ~ accum_GDD_sc + UnusualColdSum + UnusualWarmSum + 
                     voltinism + overwinteringStage
                   + (1 | unique_name__) + (1 | id_cells) + 
                     UnusualColdSum:UnusualWarmSum,
                   data = pred.df.3,
                   cov_ranef = list(unique_name = tt),
                   bayes = T)


## make df of results
rdf1.1 <- pglmm1$inla.model$summary.linear.predictor[(nrow(mdf_phylo)+1):nrow(pred.df.1),] %>% 
  mutate(V1 = pred_vals,
         V2 = "-1 - Low")

rdf2.1 <- pglmm2$inla.model$summary.linear.predictor[(nrow(mdf_phylo)+1):nrow(pred.df.1),] %>% 
  mutate(V1 = pred_vals,
         V2 = "0 - Average")

rdf3.1 <- pglmm3$inla.model$summary.linear.predictor[(nrow(mdf_phylo)+1):nrow(pred.df.1),] %>% 
  mutate(V1 = pred_vals,
         V2 = "1 - High")

rdf_total <- rbind(rdf1.1, rdf2.1, rdf3.1)

rdf_total$V2 <- factor(rdf_total$V2, levels = c("-1 - Low", "0 - Average", "1 - High"))

onset_plot1 <- ggplot(rdf_total, mapping = aes(x = V1, y = mean)) +
  geom_ribbon(mapping = aes(ymin = `0.025quant`, ymax = `0.975quant`,
                            fill = V2), alpha = 0.05 ) +
  geom_path(mapping = aes(y = `0.025quant`, color = V2), size = 0.25, linetype = 2) +
  geom_path(mapping = aes(y = `0.975quant`, color = V2), size = 0.25, linetype = 2) +
  geom_line(mapping = aes(color = V2), size = 1.05) +
  labs(x = "Unusual cold days", y = "Onset (Day of year)",
       fill = "Unusual warm days", color = "Unusual warm days") +  
  scale_color_manual(values = c("#0DFFCA", "#47B4C7", "#2A4B61")) +
  scale_fill_manual(values = c("#0DFFCA", "#47B4C7", "#2A4B61")) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    panel.border = element_blank()
  )

onset_plot1

ggsave(plot = onset_plot1, filename = "ManuscriptFigures/SI3_Fig1.png",
       height = 3, width = 6)
