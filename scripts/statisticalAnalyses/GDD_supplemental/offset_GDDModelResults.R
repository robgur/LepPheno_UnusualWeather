library(phyr)
library(car)
library(ggplot2)
library(ncf)
library(fpp2)

# get data ready for pglmm models
source('scripts/pglmm_dataSetUp_GDDModels.R')

#model_offset_with_anomaly_reduced_vif_ok
offset_m <- pglmm(offset ~ accum_GDD_sc + annualPrec + UnusualColdSum + 
                    UnusualWarmSum + voltinism + dstdoy +
                    overwinteringStage + 
                    (1 | unique_name__) + (1 | id_cells)
                  + accum_GDD_sc:UnusualColdSum + annualPrec:UnusualWarmSum
                  + UnusualColdSum:voltinism + UnusualWarmSum:voltinism,
                  data = mdf_phylo,
                  cov_ranef = list(unique_name = tt),
                  bayes = T)


# WAIC(offset_m_reduce)
# 4621
#WAIC(offset_m_noAnom_reduced)
#4673

plot_bayes(offset_m)
#ggsave(filename = "Figures/offset_posteriorDistributions_GDDModels.pdf",
#       width = 8, height = 9)

summary(offset_m)
rr2::R2(offset_m) #0.873

#save model results as table
im <- offset_m$inla.model
im_fix <- im$summary.fixed %>% 
  tibble::rownames_to_column("Effects") %>% 
  dplyr::select(Effects, mean, `0.025quant`, `0.975quant`)
sjPlot::tab_df(im_fix, title = "PGLMM offset", file = "Tables/offset_GDDPGLMM.doc")
write.csv(im_fix, file = "Tables/offsetResults_GDDPGLMM.csv", row.names = F)

## Model Assumption Checks
# residuals
resids <- residuals(offset_m)
hist(resids)

resid_offset <- mdf_phylo %>% 
  mutate(resid = resids)

## look into temporal autocorrelation
resid_offset_sort <- resid_offset[order(resid_offset$year),]

acfPlot <- acf(resid_offset_sort$resid) 

off_acf <- ggAcf(resid_offset_sort$resid) +
  theme_bw() +
  ggtitle("")

save(off_acf, file = "Figures/GDDoffset_acf.rds")

ggplot(resid_offset, mapping = aes(x = year, y = resid)) +
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


offset_correlogram_plot <- ggplot() +
  geom_line(fit_df,
            mapping = aes(x = distance, y = correlation)) +
  geom_point(fit_df,
             mapping = aes(x = distance, y = correlation, fill = p.value),
             size = 3, shape = 21, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("white", "black")) +
  labs(x = "Distance", y = "Moran's I", fill = "P-Value") +
  theme_classic()

offset_correlogram_plot # no evidence of significant spatial autocorrelation

save(offset_correlogram_plot, file = "Figures/GDDoffset_correlogram.rds")

################################################################################
################################################################################
#############Interaction Figures ###############################################
################################################################################
################################################################################
# make interaction figures
inla_mdf <- mdf_phylo %>% 
  dplyr::select(offset, 
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

#unusualcoldsum is x axis (so pred_vals),
# GDD is legend, so start with -1
#voltinism is facet, start with U
term1.U <- data.frame(offset = rep(NA, length(pred_vals)),
                      accum_GDD_sc = rep(-1, length(pred_vals)), 
                      annualPrec = rep(NA, length(pred_vals)),
                      UnusualColdSum =  pred_vals,
                      UnusualWarmSum = rep(NA, length(pred_vals)),
                      voltinism = rep("U", length(pred_vals)), 
                      dstdoy = rep(NA, length(pred_vals)),
                      overwinteringStage = rep(NA, length(pred_vals)),
                      id_cells = rep(NA, length(pred_vals)),
                      unique_name = rep(NA, length(pred_vals)))

#unusualcoldsum is x axis (so pred_vals),
#unusualwarmsum is legend, now we on to 0
#voltinism is facet so keep at "U
term2.U <- data.frame(offset = rep(NA, length(pred_vals)),
                      accum_GDD_sc = rep(0, length(pred_vals)), 
                      annualPrec = rep(NA, length(pred_vals)),
                      UnusualColdSum =  pred_vals,
                      UnusualWarmSum = rep(NA, length(pred_vals)),
                      voltinism = rep("U", length(pred_vals)), 
                      dstdoy = rep(NA, length(pred_vals)),
                      overwinteringStage = rep(NA, length(pred_vals)),
                      id_cells = rep(NA, length(pred_vals)),
                      unique_name = rep(NA, length(pred_vals)))


#unusualcoldsum is x axis (so pred_vals),
#unusualwarmsum is legend, now we on to 1
#voltinism is facet so keep at "U
term3.U <- data.frame(offset = rep(NA, length(pred_vals)),
                      accum_GDD_sc = rep(1, length(pred_vals)), 
                      annualPrec = rep(NA, length(pred_vals)),
                      UnusualColdSum =  pred_vals,
                      UnusualWarmSum = rep(NA, length(pred_vals)),
                      voltinism = rep("U", length(pred_vals)), 
                      dstdoy = rep(NA, length(pred_vals)),
                      overwinteringStage = rep(NA, length(pred_vals)),
                      id_cells = rep(NA, length(pred_vals)),
                      unique_name = rep(NA, length(pred_vals)))


## now change facet to "M"
term1.M <- data.frame(offset = rep(NA, length(pred_vals)),
                      accum_GDD_sc = rep(-1, length(pred_vals)), 
                      annualPrec = rep(NA, length(pred_vals)),
                      UnusualColdSum =  pred_vals,
                      UnusualWarmSum = rep(NA, length(pred_vals)),
                      voltinism = rep("M", length(pred_vals)), 
                      dstdoy = rep(NA, length(pred_vals)),
                      overwinteringStage = rep(NA, length(pred_vals)),
                      id_cells = rep(NA, length(pred_vals)),
                      unique_name = rep(NA, length(pred_vals)))


#M
term2.M <- data.frame(offset = rep(NA, length(pred_vals)),
                      accum_GDD_sc = rep(0, length(pred_vals)), 
                      annualPrec = rep(NA, length(pred_vals)),
                      UnusualColdSum =  pred_vals,
                      UnusualWarmSum = rep(NA, length(pred_vals)),
                      voltinism = rep("M", length(pred_vals)), 
                      dstdoy = rep(NA, length(pred_vals)),
                      overwinteringStage = rep(NA, length(pred_vals)),
                      id_cells = rep(NA, length(pred_vals)),
                      unique_name = rep(NA, length(pred_vals)))

# M
term3.M <- data.frame(offset = rep(NA, length(pred_vals)),
                      accum_GDD_sc = rep(1, length(pred_vals)), 
                      annualPrec = rep(NA, length(pred_vals)),
                      UnusualColdSum =  pred_vals,
                      UnusualWarmSum = rep(NA, length(pred_vals)),
                      voltinism = rep("M", length(pred_vals)), 
                      dstdoy = rep(NA, length(pred_vals)),
                      overwinteringStage = rep(NA, length(pred_vals)),
                      id_cells = rep(NA, length(pred_vals)),
                      unique_name = rep(NA, length(pred_vals)))

pred.df.1.U <- rbind(inla_mdf, term1.U)
pred.df.2.U <- rbind(inla_mdf, term2.U)
pred.df.3.U <- rbind(inla_mdf, term3.U)

pred.df.1.M <- rbind(inla_mdf, term1.M)
pred.df.2.M <- rbind(inla_mdf, term2.M)
pred.df.3.M <- rbind(inla_mdf, term3.M)

# fit models
pglmm1.U <- pglmm(offset ~ accum_GDD_sc + annualPrec + UnusualColdSum + 
                    UnusualWarmSum + voltinism + dstdoy +
                    overwinteringStage + 
                    (1 | unique_name__) + (1 | id_cells)
                  + accum_GDD_sc:UnusualColdSum + annualPrec:UnusualWarmSum
                  + UnusualColdSum:voltinism + UnusualWarmSum:voltinism,
                  data = pred.df.1.U,
                  cov_ranef = list(unique_name = tt),
                  bayes = T)

pglmm2.U <- pglmm(offset ~ accum_GDD_sc + annualPrec + UnusualColdSum + 
                    UnusualWarmSum + voltinism + dstdoy +
                    overwinteringStage + 
                    (1 | unique_name__) + (1 | id_cells)
                  + accum_GDD_sc:UnusualColdSum + annualPrec:UnusualWarmSum
                  + UnusualColdSum:voltinism + UnusualWarmSum:voltinism,
                  data = pred.df.2.U,
                  cov_ranef = list(unique_name = tt),
                  bayes = T)


pglmm3.U <- pglmm(offset ~ accum_GDD_sc + annualPrec + UnusualColdSum + 
                    UnusualWarmSum + voltinism + dstdoy +
                    overwinteringStage + 
                    (1 | unique_name__) + (1 | id_cells)
                  + accum_GDD_sc:UnusualColdSum + annualPrec:UnusualWarmSum
                  + UnusualColdSum:voltinism + UnusualWarmSum:voltinism,
                  data = pred.df.3.U,
                  cov_ranef = list(unique_name = tt),
                  bayes = T)


pglmm1.M <- pglmm(offset ~ accum_GDD_sc + annualPrec + UnusualColdSum + 
                   UnusualWarmSum + voltinism + dstdoy +
                   overwinteringStage + 
                   (1 | unique_name__) + (1 | id_cells)
                 + accum_GDD_sc:UnusualColdSum + annualPrec:UnusualWarmSum
                 + UnusualColdSum:voltinism + UnusualWarmSum:voltinism,
                 data = pred.df.1.M,
                 cov_ranef = list(unique_name = tt),
                 bayes = T)


pglmm2.M <- pglmm(offset ~ accum_GDD_sc + annualPrec + UnusualColdSum + 
                    UnusualWarmSum + voltinism + dstdoy +
                    overwinteringStage + 
                    (1 | unique_name__) + (1 | id_cells)
                  + accum_GDD_sc:UnusualColdSum + annualPrec:UnusualWarmSum
                  + UnusualColdSum:voltinism + UnusualWarmSum:voltinism,
                  data = pred.df.2.M,
                  cov_ranef = list(unique_name = tt),
                  bayes = T)

pglmm3.M <- pglmm(offset ~ accum_GDD_sc + annualPrec + UnusualColdSum + 
                    UnusualWarmSum + voltinism + dstdoy +
                    overwinteringStage + 
                    (1 | unique_name__) + (1 | id_cells)
                  + accum_GDD_sc:UnusualColdSum + annualPrec:UnusualWarmSum
                  + UnusualColdSum:voltinism + UnusualWarmSum:voltinism,
                  data = pred.df.3.M,
                  cov_ranef = list(unique_name = tt),
                  bayes = T)


## make df of results
rdf1.1 <- pglmm1.U$inla.model$summary.linear.predictor[(nrow(mdf_phylo)+1):nrow(pred.df.1.U),] %>% 
  mutate(V1 = pred_vals,
         V2 = "-1 - Low",
         V3 = "Univoltine")

rdf2.1 <- pglmm2.U$inla.model$summary.linear.predictor[(nrow(mdf_phylo)+1):nrow(pred.df.1.U),] %>% 
  mutate(V1 = pred_vals,
         V2 = "0 - Average",
         V3 = "Univoltine")

rdf3.1 <- pglmm3.U$inla.model$summary.linear.predictor[(nrow(mdf_phylo)+1):nrow(pred.df.1.U),] %>% 
  mutate(V1 = pred_vals,
         V2 = "1 - High",
         V3 = "Univoltine")

rdf1.2 <- pglmm1.M$inla.model$summary.linear.predictor[(nrow(mdf_phylo)+1):nrow(pred.df.1.M),] %>% 
  mutate(V1 = pred_vals,
         V2 = "-1 - Low",
         V3 = "Multivoltine")

rdf2.2 <- pglmm2.M$inla.model$summary.linear.predictor[(nrow(mdf_phylo)+1):nrow(pred.df.1.M),] %>% 
  mutate(V1 = pred_vals,
         V2 = "0 - Average",
         V3 = "Multivoltine")

rdf3.2 <- pglmm3.M$inla.model$summary.linear.predictor[(nrow(mdf_phylo)+1):nrow(pred.df.1.M),] %>% 
  mutate(V1 = pred_vals,
         V2 = "1 - High",
         V3 = "Multivoltine")


rdf_total <- rbind(rdf1.1, rdf2.1, rdf3.1,
                   rdf1.2, rdf2.2, rdf3.2)

rdf_total$V2 <- factor(rdf_total$V2, levels = c("-1 - Low", "0 - Average", "1 - High"))
rdf_total$V3 <- factor(rdf_total$V3, levels = c("Univoltine", "Multivoltine"))

cold_GDD_volt_plot <- ggplot(rdf_total, mapping = aes(x = V1, y = mean)) +
  geom_ribbon(mapping = aes(ymin = `0.025quant`, ymax = `0.975quant`,
                            fill = V2), alpha = 0.05 ) +
  geom_path(mapping = aes(y = `0.025quant`, color = V2), size = 0.25, linetype = 2) +
  geom_path(mapping = aes(y = `0.975quant`, color = V2), size = 0.25, linetype = 2) +
  geom_line(mapping = aes(color = V2), size = 1.05) +
  labs(x = "Unusual cold days", y = "Offset (Day of year)", fill = "GDD", color = "GDD") + 
  scale_color_manual(values = c("#046865", "#FDCA9B", "#8C032C")) +
  scale_fill_manual(values = c("#046865", "#FDCA9B", "#8C032C")) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    panel.border = element_blank()
  ) +
  facet_wrap(~V3) 

cold_GDD_volt_plot

ggsave(plot = cold_GDD_volt_plot, filename = "Figures/offest_cold_GDD_volt_interaction.png",
       height = 3, width = 7)


################################################################################
################################################################################
#############Interaction Figure 2 ###############################################
################################################################################
################################################################################
# make interaction figures Prec : UnusualWarmSum

#unusualWarmSum is x axis (so pred_vals),
# Prec is legend, so start with -1
#voltinism is facet, start with U
term1.U <- data.frame(offset = rep(NA, length(pred_vals)),
                      accum_GDD_sc = rep(NA, length(pred_vals)), 
                      annualPrec = rep(-1, length(pred_vals)),
                      UnusualColdSum = rep(NA, length(pred_vals)),
                      UnusualWarmSum = pred_vals,
                      voltinism = rep("U", length(pred_vals)), 
                      dstdoy = rep(NA, length(pred_vals)),
                      overwinteringStage = rep(NA, length(pred_vals)),
                      id_cells = rep(NA, length(pred_vals)),
                      unique_name = rep(NA, length(pred_vals)))

#unusualWarmSum is x axis (so pred_vals),
# Prec is legend, so start with 0
#voltinism is facet, start with U
term1.U <- data.frame(offset = rep(NA, length(pred_vals)),
                      accum_GDD_sc = rep(NA, length(pred_vals)), 
                      annualPrec = rep(0, length(pred_vals)),
                      UnusualColdSum = rep(NA, length(pred_vals)),
                      UnusualWarmSum = pred_vals,
                      voltinism = rep("U", length(pred_vals)), 
                      dstdoy = rep(NA, length(pred_vals)),
                      overwinteringStage = rep(NA, length(pred_vals)),
                      id_cells = rep(NA, length(pred_vals)),
                      unique_name = rep(NA, length(pred_vals)))


#unusualWarmSum is x axis (so pred_vals),
# Prec is legend, so start with 1
#voltinism is facet, start with U
term1.U <- data.frame(offset = rep(NA, length(pred_vals)),
                      accum_GDD_sc = rep(NA, length(pred_vals)), 
                      annualPrec = rep(1, length(pred_vals)),
                      UnusualColdSum = rep(NA, length(pred_vals)),
                      UnusualWarmSum = pred_vals,
                      voltinism = rep("U", length(pred_vals)), 
                      dstdoy = rep(NA, length(pred_vals)),
                      overwinteringStage = rep(NA, length(pred_vals)),
                      id_cells = rep(NA, length(pred_vals)),
                      unique_name = rep(NA, length(pred_vals)))


## now change facet to "M"
term1.M <- data.frame(offset = rep(NA, length(pred_vals)),
                     accum_GDD_sc = rep(NA, length(pred_vals)), 
                     annualPrec = rep(-1, length(pred_vals)),
                     UnusualColdSum = rep(NA, length(pred_vals)),
                     UnusualWarmSum = pred_vals,
                     voltinism = rep("M", length(pred_vals)), 
                     dstdoy = rep(NA, length(pred_vals)),
                     overwinteringStage = rep(NA, length(pred_vals)),
                     id_cells = rep(NA, length(pred_vals)),
                     unique_name = rep(NA, length(pred_vals)))


#M
term2.M <- data.frame(offset = rep(NA, length(pred_vals)),
                      accum_GDD_sc = rep(NA, length(pred_vals)), 
                      annualPrec = rep(0, length(pred_vals)),
                      UnusualColdSum = rep(NA, length(pred_vals)),
                      UnusualWarmSum = pred_vals,
                      voltinism = rep("M", length(pred_vals)), 
                      dstdoy = rep(NA, length(pred_vals)),
                      overwinteringStage = rep(NA, length(pred_vals)),
                      id_cells = rep(NA, length(pred_vals)),
                      unique_name = rep(NA, length(pred_vals)))

# M
term3.M <- data.frame(offset = rep(NA, length(pred_vals)),
                      accum_GDD_sc = rep(NA, length(pred_vals)), 
                      annualPrec = rep(1, length(pred_vals)),
                      UnusualColdSum = rep(NA, length(pred_vals)),
                      UnusualWarmSum = pred_vals,
                      voltinism = rep("M", length(pred_vals)), 
                      dstdoy = rep(NA, length(pred_vals)),
                      overwinteringStage = rep(NA, length(pred_vals)),
                      id_cells = rep(NA, length(pred_vals)),
                      unique_name = rep(NA, length(pred_vals)))

pred.df.1.U <- rbind(inla_mdf, term1.U)
pred.df.2.U <- rbind(inla_mdf, term2.U)
pred.df.3.U <- rbind(inla_mdf, term3.U)

pred.df.1.M <- rbind(inla_mdf, term1.M)
pred.df.2.M <- rbind(inla_mdf, term2.M)
pred.df.3.M <- rbind(inla_mdf, term3.M)

# fit models
pglmm1.U <- pglmm(offset ~ accum_GDD_sc + annualPrec + UnusualColdSum + 
                    UnusualWarmSum + voltinism + dstdoy +
                    overwinteringStage + 
                    (1 | unique_name__) + (1 | id_cells)
                  + accum_GDD_sc:UnusualColdSum + annualPrec:UnusualWarmSum
                  + UnusualColdSum:voltinism + UnusualWarmSum:voltinism,
                  data = pred.df.1.U,
                  cov_ranef = list(unique_name = tt),
                  bayes = T)

pglmm2.U <- pglmm(offset ~ accum_GDD_sc + annualPrec + UnusualColdSum + 
                    UnusualWarmSum + voltinism + dstdoy +
                    overwinteringStage + 
                    (1 | unique_name__) + (1 | id_cells)
                  + accum_GDD_sc:UnusualColdSum + annualPrec:UnusualWarmSum
                  + UnusualColdSum:voltinism + UnusualWarmSum:voltinism,
                  data = pred.df.2.U,
                  cov_ranef = list(unique_name = tt),
                  bayes = T)


pglmm3.U <- pglmm(offset ~ accum_GDD_sc + annualPrec + UnusualColdSum + 
                    UnusualWarmSum + voltinism + dstdoy +
                    overwinteringStage + 
                    (1 | unique_name__) + (1 | id_cells)
                  + accum_GDD_sc:UnusualColdSum + annualPrec:UnusualWarmSum
                  + UnusualColdSum:voltinism + UnusualWarmSum:voltinism,
                  data = pred.df.3.U,
                  cov_ranef = list(unique_name = tt),
                  bayes = T)


pglmm1.M <- pglmm(offset ~ accum_GDD_sc + annualPrec + UnusualColdSum + 
                    UnusualWarmSum + voltinism + dstdoy +
                    overwinteringStage + 
                    (1 | unique_name__) + (1 | id_cells)
                  + accum_GDD_sc:UnusualColdSum + annualPrec:UnusualWarmSum
                  + UnusualColdSum:voltinism + UnusualWarmSum:voltinism,
                  data = pred.df.1.M,
                  cov_ranef = list(unique_name = tt),
                  bayes = T)


pglmm2.M <- pglmm(offset ~ accum_GDD_sc + annualPrec + UnusualColdSum + 
                    UnusualWarmSum + voltinism + dstdoy +
                    overwinteringStage + 
                    (1 | unique_name__) + (1 | id_cells)
                  + accum_GDD_sc:UnusualColdSum + annualPrec:UnusualWarmSum
                  + UnusualColdSum:voltinism + UnusualWarmSum:voltinism,
                  data = pred.df.2.M,
                  cov_ranef = list(unique_name = tt),
                  bayes = T)

pglmm3.M <- pglmm(offset ~ accum_GDD_sc + annualPrec + UnusualColdSum + 
                    UnusualWarmSum + voltinism + dstdoy +
                    overwinteringStage + 
                    (1 | unique_name__) + (1 | id_cells)
                  + accum_GDD_sc:UnusualColdSum + annualPrec:UnusualWarmSum
                  + UnusualColdSum:voltinism + UnusualWarmSum:voltinism,
                  data = pred.df.3.M,
                  cov_ranef = list(unique_name = tt),
                  bayes = T)


## make df of results
rdf1.1 <- pglmm1.U$inla.model$summary.linear.predictor[(nrow(mdf_phylo)+1):nrow(pred.df.1.U),] %>% 
  mutate(V1 = pred_vals,
         V2 = "-1 - Low",
         V3 = "Univoltine")

rdf2.1 <- pglmm2.U$inla.model$summary.linear.predictor[(nrow(mdf_phylo)+1):nrow(pred.df.1.U),] %>% 
  mutate(V1 = pred_vals,
         V2 = "0 - Average",
         V3 = "Univoltine")

rdf3.1 <- pglmm3.U$inla.model$summary.linear.predictor[(nrow(mdf_phylo)+1):nrow(pred.df.1.U),] %>% 
  mutate(V1 = pred_vals,
         V2 = "1 - High",
         V3 = "Univoltine")

rdf1.2 <- pglmm1.M$inla.model$summary.linear.predictor[(nrow(mdf_phylo)+1):nrow(pred.df.1.M),] %>% 
  mutate(V1 = pred_vals,
         V2 = "-1 - Low",
         V3 = "Multivoltine")

rdf2.2 <- pglmm2.M$inla.model$summary.linear.predictor[(nrow(mdf_phylo)+1):nrow(pred.df.1.M),] %>% 
  mutate(V1 = pred_vals,
         V2 = "0 - Average",
         V3 = "Multivoltine")

rdf3.2 <- pglmm3.M$inla.model$summary.linear.predictor[(nrow(mdf_phylo)+1):nrow(pred.df.1.M),] %>% 
  mutate(V1 = pred_vals,
         V2 = "1 - High",
         V3 = "Multivoltine")


rdf_total <- rbind(rdf1.1, rdf2.1, rdf3.1,
                   rdf1.2, rdf2.2, rdf3.2)

rdf_total$V2 <- factor(rdf_total$V2, levels = c("-1 - Low", "0 - Average", "1 - High"))
rdf_total$V3 <- factor(rdf_total$V3, levels = c("Univoltine", "Multivoltine"))

prec_warm_volt_plot <- ggplot(rdf_total, mapping = aes(x = V1, y = mean)) +
  geom_ribbon(mapping = aes(ymin = `0.025quant`, ymax = `0.975quant`,
                            fill = V2), alpha = 0.05 ) +
  geom_path(mapping = aes(y = `0.025quant`, color = V2), size = 0.25, linetype = 2) +
  geom_path(mapping = aes(y = `0.975quant`, color = V2), size = 0.25, linetype = 2) +
  geom_line(mapping = aes(color = V2), size = 1.05) +
  labs(x = "Unusual warm days", y = "Offset (Day of year)", 
       fill = "Annual prec", color = "Annual prec") + 
  scale_color_manual(values = c("#003249", "#8390FA", "#3BF4FB")) +
  scale_fill_manual(values = c("#003249", "#8390FA",  "#3BF4FB")) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    panel.border = element_blank()
  ) +
  facet_wrap(~V3) 

prec_warm_volt_plot

ggsave(plot = prec_warm_volt_plot, filename = "Figures/offest_GDDModel_prec_warm_volt.png",
       height = 3, width = 7)


cowplot::plot_grid(cold_GDD_volt_plot, prec_warm_volt_plot,
                   nrow = 2, ncol = 1,
                   labels = c("A", "B"))

ggsave(filename = "ManuscriptFigures/SI3_Fig2.png",
       height = 6, width = 7)
