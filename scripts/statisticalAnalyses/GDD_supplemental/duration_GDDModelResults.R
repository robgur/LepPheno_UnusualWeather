library(phyr)
library(car)
library(ggplot2)
library(ncf)
library(fpp2)

# get data ready for pglmm models
source('scripts/pglmm_dataSetUp_GDDModels.R')

#model_dur_with_anomaly_reduced_vif_ok
dur_m <- pglmm(duration ~ accum_GDD_sc + annualPrec + 
                        UnusualColdSum + UnusualWarmSum + voltinism +
                        dstdoy + overwinteringStage + 
                        (1 | unique_name__) + (1 | id_cells) + 
                        accum_GDD_sc:voltinism 
                      + annualPrec:UnusualWarmSum + UnusualColdSum:voltinism, 
                      data = mdf_phylo,
                      cov_ranef = list(unique_name = tt),
                      bayes = T)


# WAIC(dur_m_reduce)
# 4621
#WAIC(dur_m_noAnom_reduced)
#4673

plot_bayes(dur_m)
#ggsave(filename = "Figures/dur_posteriorDistributions_GDDModels.pdf",
#       width = 8, height = 9)

summary(dur_m)
rr2::R2(dur_m) #0.836

#save model results as table
im <- dur_m$inla.model
im_fix <- im$summary.fixed %>% 
  tibble::rownames_to_column("Effects") %>% 
  dplyr::select(Effects, mean, `0.025quant`, `0.975quant`)
sjPlot::tab_df(im_fix, title = "PGLMM dur", file = "Tables/dur_GDDPGLMM.doc")
write.csv(im_fix, file = "Tables/durResults_GDDPGLMM.csv", row.names = F)

## Model Assumption Checks
# residuals
resids <- residuals(dur_m)
hist(resids)

resid_dur <- mdf_phylo %>% 
  mutate(resid = resids)

## look into temporal autocorrelation
resid_dur_sort <- resid_dur[order(resid_dur$year),]

acfPlot <- acf(resid_dur_sort$resid) 

dur_acf <- ggAcf(resid_dur_sort$resid) +
  theme_bw() +
  ggtitle("")

save(dur_acf, file = "Figures/GDDdur_acf.rds")

ggplot(resid_dur, mapping = aes(x = year, y = resid)) +
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


dur_correlogram_plot <- ggplot() +
  geom_line(fit_df,
            mapping = aes(x = distance, y = correlation)) +
  geom_point(fit_df,
             mapping = aes(x = distance, y = correlation, fill = p.value),
             size = 3, shape = 21, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("white", "black")) +
  labs(x = "Distance", y = "Moran's I", fill = "P-Value") +
  theme_classic()

dur_correlogram_plot # no evidence of significant spatial autocorrelation

save(dur_correlogram_plot, file = "Figures/GDDdur_correlogram.rds")

################################################################################
################################################################################
#############Interaction Figures ###############################################
################################################################################
################################################################################
# make interaction figures
inla_mdf <- mdf_phylo %>% 
  dplyr::select(duration, 
                accum_GDD_sc, 
                annualPrec, 
                UnusualColdSum, 
                UnusualWarmSum, 
                voltinism, 
                dstdoy,
                overwinteringStage, 
                id_cells,
                unique_name)

pred_vals <- -3:3

# unusualcoldsum is x axis (so pred_vals),
#voltinism is legend, start with U
term1.U <- data.frame(duration = rep(NA, length(pred_vals)),
                      accum_GDD_sc = rep(NA, length(pred_vals)), 
                      annualPrec = rep(NA, length(pred_vals)),
                      UnusualColdSum =  pred_vals,
                      UnusualWarmSum = rep(NA, length(pred_vals)),
                      voltinism = rep("U", length(pred_vals)), 
                      dstdoy = rep(NA, length(pred_vals)),
                      overwinteringStage = rep(NA, length(pred_vals)),
                      id_cells = rep(NA, length(pred_vals)),
                      unique_name = rep(NA, length(pred_vals)))

# unusualcoldsum is x axis (so pred_vals),
#voltinism is legend, move to M
term1.M <- data.frame(duration = rep(NA, length(pred_vals)),
                      accum_GDD_sc = rep(NA, length(pred_vals)), 
                      annualPrec = rep(NA, length(pred_vals)),
                      UnusualColdSum =  pred_vals,
                      UnusualWarmSum = rep(NA, length(pred_vals)),
                      voltinism = rep("M", length(pred_vals)), 
                      dstdoy = rep(NA, length(pred_vals)),
                      overwinteringStage = rep(NA, length(pred_vals)),
                      id_cells = rep(NA, length(pred_vals)),
                      unique_name = rep(NA, length(pred_vals)))


pred.df.1.U <- rbind(inla_mdf, term1.U)
pred.df.1.M <- rbind(inla_mdf, term1.M)

# fit models
pglmm1.U <- pglmm(duration ~ accum_GDD_sc + annualPrec + 
                    UnusualColdSum + UnusualWarmSum + voltinism +
                    dstdoy + overwinteringStage + 
                    (1 | unique_name__) + (1 | id_cells) + 
                    accum_GDD_sc:voltinism 
                  + annualPrec:UnusualWarmSum + UnusualColdSum:voltinism, 
                  data = pred.df.1.U,
                  cov_ranef = list(unique_name = tt),
                  bayes = T)


pglmm1.M <- pglmm(duration ~ accum_GDD_sc + annualPrec + 
                    UnusualColdSum + UnusualWarmSum + voltinism +
                    dstdoy + overwinteringStage + 
                    (1 | unique_name__) + (1 | id_cells) + 
                    accum_GDD_sc:voltinism 
                  + annualPrec:UnusualWarmSum + UnusualColdSum:voltinism, 
                  data = pred.df.1.M,
                  cov_ranef = list(unique_name = tt),
                  bayes = T)


## make df of results
rdf1.1 <- pglmm1.U$inla.model$summary.linear.predictor[(nrow(mdf_phylo)+1):nrow(pred.df.1.U),] %>% 
  mutate(V1 = pred_vals,
         V3 = "Univoltine")

rdf1.2 <- pglmm1.M$inla.model$summary.linear.predictor[(nrow(mdf_phylo)+1):nrow(pred.df.1.M),] %>% 
  mutate(V1 = pred_vals,
         V3 = "Multivoltine")

rdf_total <- rbind(rdf1.1,
                   rdf1.2)

rdf_total$V3 <- factor(rdf_total$V3, levels = c("Univoltine", "Multivoltine"))

ColdDay_Volt <- ggplot(rdf_total, mapping = aes(x = V1, y = mean)) +
  geom_ribbon(mapping = aes(ymin = `0.025quant`, ymax = `0.975quant`,
                            fill = V3), alpha = 0.05 ) +
  geom_path(mapping = aes(y = `0.025quant`, color = V3), size = 0.25, linetype = 2) +
  geom_path(mapping = aes(y = `0.975quant`, color = V3), size = 0.25, linetype = 2) +
  geom_line(mapping = aes(color = V3), size = 1.05) +
  labs(x = "Unusual cold days", y = "Duration (Days)", fill = "Voltinism", color = "Voltinism") + 
  scale_color_manual(values = c("#CBA328","#3A445D")) +
  scale_fill_manual(values = c("#CBA328","#3A445D")) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    panel.border = element_blank()
  )

ColdDay_Volt

ggsave(plot = ColdDay_Volt, filename = "Figures/durationGDD_coldDay_voltinism.png",
       height = 3, width = 7)


################################################################################
################################################################################
#############Interaction Figure 2 ###############################################
################################################################################
################################################################################
# make interaction figures Prec : UnusualWarmSum

#unusualWarmSum is x axis (so pred_vals),
# Prec is legend, so start with -1
term1 <- data.frame(duration = rep(NA, length(pred_vals)),
                      accum_GDD_sc = rep(NA, length(pred_vals)), 
                      annualPrec = rep(-1, length(pred_vals)),
                      UnusualColdSum = rep(NA, length(pred_vals)),
                      UnusualWarmSum = pred_vals,
                      voltinism = rep(NA, length(pred_vals)), 
                      dstdoy = rep(NA, length(pred_vals)),
                      overwinteringStage = rep(NA, length(pred_vals)),
                      id_cells = rep(NA, length(pred_vals)),
                      unique_name = rep(NA, length(pred_vals)))

# Prec is legend, so move to 0
term2 <- data.frame(duration = rep(NA, length(pred_vals)),
                    accum_GDD_sc = rep(NA, length(pred_vals)), 
                    annualPrec = rep(0,length(pred_vals)),
                    UnusualColdSum = rep(NA, length(pred_vals)),
                    UnusualWarmSum = pred_vals,
                    voltinism = rep(NA, length(pred_vals)), 
                    dstdoy = rep(NA, length(pred_vals)),
                    overwinteringStage = rep(NA, length(pred_vals)),
                    id_cells = rep(NA, length(pred_vals)),
                    unique_name = rep(NA, length(pred_vals)))


# Prec is legend, so move to 1
term3 <- data.frame(duration = rep(NA, length(pred_vals)),
                    accum_GDD_sc = rep(NA, length(pred_vals)), 
                    annualPrec = rep(1,length(pred_vals)),
                    UnusualColdSum = rep(NA, length(pred_vals)),
                    UnusualWarmSum = pred_vals,
                    voltinism = rep(NA, length(pred_vals)), 
                    dstdoy = rep(NA, length(pred_vals)),
                    overwinteringStage = rep(NA, length(pred_vals)),
                    id_cells = rep(NA, length(pred_vals)),
                    unique_name = rep(NA, length(pred_vals)))

pred.df.1 <- rbind(inla_mdf, term1)
pred.df.2 <- rbind(inla_mdf, term2)
pred.df.3 <- rbind(inla_mdf, term3)

# fit models
pglmm1 <- pglmm(duration ~ accum_GDD_sc + annualPrec + UnusualColdSum + 
                    UnusualWarmSum + voltinism + dstdoy +
                    overwinteringStage + 
                    (1 | unique_name__) + (1 | id_cells)
                  + accum_GDD_sc:UnusualColdSum + annualPrec:UnusualWarmSum
                  + UnusualColdSum:voltinism + UnusualWarmSum:voltinism,
                  data = pred.df.1,
                  cov_ranef = list(unique_name = tt),
                  bayes = T)

pglmm2 <- pglmm(duration ~ accum_GDD_sc + annualPrec + UnusualColdSum + 
                    UnusualWarmSum + voltinism + dstdoy +
                    overwinteringStage + 
                    (1 | unique_name__) + (1 | id_cells)
                  + accum_GDD_sc:UnusualColdSum + annualPrec:UnusualWarmSum
                  + UnusualColdSum:voltinism + UnusualWarmSum:voltinism,
                  data = pred.df.2,
                  cov_ranef = list(unique_name = tt),
                  bayes = T)


pglmm3 <- pglmm(duration ~ accum_GDD_sc + annualPrec + UnusualColdSum + 
                    UnusualWarmSum + voltinism + dstdoy +
                    overwinteringStage + 
                    (1 | unique_name__) + (1 | id_cells)
                  + accum_GDD_sc:UnusualColdSum + annualPrec:UnusualWarmSum
                  + UnusualColdSum:voltinism + UnusualWarmSum:voltinism,
                  data = pred.df.3,
                  cov_ranef = list(unique_name = tt),
                  bayes = T)

## make df of results
rdf1.1 <- pglmm1$inla.model$summary.linear.predictor[(nrow(mdf_phylo)+1):nrow(pred.df.1),] %>% 
  mutate(V1 = pred_vals,
         V2 = "-1 - Low")

rdf2.1 <- pglmm2$inla.model$summary.linear.predictor[(nrow(mdf_phylo)+1):nrow(pred.df.2),] %>% 
  mutate(V1 = pred_vals,
         V2 = "0 - Average")

rdf3.1 <- pglmm3$inla.model$summary.linear.predictor[(nrow(mdf_phylo)+1):nrow(pred.df.3),] %>% 
  mutate(V1 = pred_vals,
         V2 = "1 - High")

rdf_total <- rbind(rdf1.1, rdf2.1, rdf3.1)

rdf_total$V2 <- factor(rdf_total$V2, levels = c("-1 - Low", "0 - Average", "1 - High"))

WarmDay_Prec <- ggplot(rdf_total, mapping = aes(x = V1, y = mean)) +
  geom_ribbon(mapping = aes(ymin = `0.025quant`, ymax = `0.975quant`,
                            fill = V2), alpha = 0.05 ) +
  geom_path(mapping = aes(y = `0.025quant`, color = V2), size = 0.25, linetype = 2) +
  geom_path(mapping = aes(y = `0.975quant`, color = V2), size = 0.25, linetype = 2) +
  geom_line(mapping = aes(color = V2), size = 1.05) +
  labs(x = "Unusual warm days", y = "Duration (Days)", 
       fill = "Annual prec", color = "Annual prec") + 
  scale_color_manual(values = c("#003249", "#8390FA", "#3BF4FB")) +
  scale_fill_manual(values = c("#003249", "#8390FA",  "#3BF4FB")) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    panel.border = element_blank()
  ) 

WarmDay_Prec

ggsave(plot = WarmDay_Prec, filename = "Figures/Duration_GDDModel_WarmDay_prec.png",
       height = 3, width = 7)

cowplot::plot_grid(ColdDay_Volt, WarmDay_Prec,
                   nrow = 2, ncol = 1,
                   labels = c("A", "B"))

ggsave(filename = "ManuscriptFigures/SI3_Fig3.png",
       height = 6, width = 7)
