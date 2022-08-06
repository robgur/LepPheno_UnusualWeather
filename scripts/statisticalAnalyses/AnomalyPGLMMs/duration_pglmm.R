library(phyr)
library(car)
library(sjPlot)
library(ggplot2)
library(geoR)
library(ncf)
library(fpp2)

# get data ready for pglmm models
source('scripts/pglmm_dataSetUp.R')

mdf_phylo <- mdf_phylo 

# duration model
duration_m <- pglmm(duration ~ annualTemp + annualPrec + 
                      UnusualColdSum + UnusualWarmSum + 
                      voltinism + Seas + dstdoy + ndstcol + overwinteringStage + 
                      (1 | unique_name__) + 
                      annualTemp:UnusualColdSum + annualTemp:UnusualWarmSum + 
                      annualPrec:UnusualWarmSum + 
                      UnusualColdSum:UnusualWarmSum + UnusualColdSum:voltinism,
              data = mdf_phylo, 
              cov_ranef = list(unique_name = tt),
              bayes = T)

duration_m_noAnom <- pglmm(duration ~ annualTemp + annualPrec + 
                      voltinism + Seas + dstdoy + ndstcol + overwinteringStage + 
                      (1 | unique_name__),
                    data = mdf_phylo, 
                    cov_ranef = list(unique_name = tt),
                    bayes = T)
summary(duration_m)
rr2::R2(duration_m) # 0.778

rr2::R2(duration_m_noAnom)# 0.728

phyr::plot_bayes(duration_m)

#save model results as table
im <- duration_m$inla.model
im_fix <- im$summary.fixed %>% 
  tibble::rownames_to_column("Effects") %>% 
  dplyr::select(Effects, mean, `0.025quant`, `0.975quant`)
tab_df(im_fix, title = "PGLMM duration", file = "results/duration_PGLMM.doc")

## Model Assumption Checks
# residuals
resids <- residuals(duration_m)
hist(resids)
qqnorm(resids)
qqline(resids)

# spatial autocorrelation
# spatially visualize residuals
resid_duration <- mdf_phylo %>% 
  mutate(resid = resids)

ggplot(resid_duration, mapping = aes(x = lon, y = lat, color = resid)) +
  geom_point(size = 10) +
  scale_color_viridis_c(option = "turbo")

## look into temporal autocorrelation
resid_duration_sort <- resid_duration[order(resid_duration$year),]

acf(resid_duration_sort$resid) 

dur_acf <- ggAcf(resid_duration_sort$resid) +
  theme_bw() +
  ggtitle("")

save(dur_acf, file = "Figures/duration_acf.rds")

# variogram approach
jitter_dist = 125000 # 125 km

gdf <- as.geodata(resid_duration, coords.col = c("lon", "lat"), data.col = "resid")
gdf <-  geoR::jitterDupCoords(gdf, max = jitter_dist)
v1 <- variog(gdf, trend = "1st")

plot(v1)

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


duration_correlogram_plot <- ggplot() +
  geom_line(fit_df,
            mapping = aes(x = distance, y = correlation)) +
  geom_point(fit_df,
             mapping = aes(x = distance, y = correlation, fill = p.value),
             size = 3, shape = 21, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("white", "black")) +
  labs(x = "Distance", y = "Moran's I", fill = "P-Value") +
  theme_classic()

duration_correlogram_plot # no evidence of significant spatial autocorrelation

save(duration_correlogram_plot, file = "Figures/duration_correlogram.rds")

################################################################################
################################################################################
#############Interaction Figures ###############################################
################################################################################
################################################################################
# make interaction figures rob included in slack
## plot_model(dur_m, type = "pred", terms = c("UnusualColdSum", "UnusualWarmSum","annualTemp"))

inla_mdf <- mdf_phylo %>% 
  dplyr::select(duration, 
                annualTemp, 
                annualPrec, 
                UnusualColdSum, 
                year, 
                UnusualWarmSum, 
                voltinism, 
                Seas, 
                tempSeas, 
                dstdoy,
                overwinteringStage, 
                diurnality,
                id_cells,
                unique_name)

pred_vals <- -3:3

#unusualcoldsum is x axis (so pred_vals),
#unusualwarmsum is legend, so start with -1
#annual temp is facet, start with -1
term1.1 <- data.frame(duration = rep(NA, length(pred_vals)),
                    annualTemp = rep(-1, length(pred_vals)), 
                    annualPrec = rep(NA, length(pred_vals)),
                    UnusualColdSum = pred_vals, 
                    year = rep(NA, length(pred_vals)), 
                    UnusualWarmSum = rep(-1, length(pred_vals)),
                    voltinism = rep(NA, length(pred_vals)), 
                    Seas = rep(NA, length(pred_vals)), 
                    tempSeas = rep(NA, length(pred_vals)), 
                    dstdoy = rep(NA, length(pred_vals)),
                    overwinteringStage = rep(NA, length(pred_vals)), 
                    diurnality = rep(NA, length(pred_vals)),
                    id_cells = rep(NA, length(pred_vals)),
                    unique_name = rep(NA, length(pred_vals)))
  
#unusualcoldsum is x axis (so pred_vals),
#unusualwarmsum is legend, now we on to 0
#annual temp is facet, start with -1 # keep the same
term2.1 <- data.frame(duration = rep(NA, length(pred_vals)),
                    annualTemp = rep(-1, length(pred_vals)), 
                    annualPrec = rep(NA, length(pred_vals)),
                    UnusualColdSum = pred_vals, 
                    year = rep(NA, length(pred_vals)), 
                    UnusualWarmSum = rep(0, length(pred_vals)),
                    voltinism = rep(NA, length(pred_vals)), 
                    Seas = rep(NA, length(pred_vals)), 
                    tempSeas = rep(NA, length(pred_vals)), 
                    dstdoy = rep(NA, length(pred_vals)),
                    overwinteringStage = rep(NA, length(pred_vals)), 
                    diurnality = rep(NA, length(pred_vals)),
                    id_cells = rep(NA, length(pred_vals)),
                    unique_name = rep(NA, length(pred_vals)))

#unusualcoldsum is x axis (so pred_vals),
#unusualwarmsum is legend, now we on to 1
#annual temp is facet, start with -1 # keep the same
term3.1 <- data.frame(duration = rep(NA, length(pred_vals)),
                    annualTemp = rep(-1, length(pred_vals)), 
                    annualPrec = rep(NA, length(pred_vals)),
                    UnusualColdSum = pred_vals, 
                    year = rep(NA, length(pred_vals)), 
                    UnusualWarmSum = rep(1, length(pred_vals)),
                    voltinism = rep(NA, length(pred_vals)), 
                    Seas = rep(NA, length(pred_vals)), 
                    tempSeas = rep(NA, length(pred_vals)), 
                    dstdoy = rep(NA, length(pred_vals)),
                    overwinteringStage = rep(NA, length(pred_vals)), 
                    diurnality = rep(NA, length(pred_vals)),
                    id_cells = rep(NA, length(pred_vals)),
                    unique_name = rep(NA, length(pred_vals)))

### changed all the unusalwarm with annual temp being -1, so neeed to change annual temp to 0
#unusualcoldsum is x axis (so pred_vals),
#unusualwarmsum is legend, now we back to -1
#annual temp is facet, start with 0 # keep the same
term1.2 <- data.frame(duration = rep(NA, length(pred_vals)),
                    annualTemp = rep(0, length(pred_vals)), 
                    annualPrec = rep(NA, length(pred_vals)),
                    UnusualColdSum = pred_vals, 
                    year = rep(NA, length(pred_vals)), 
                    UnusualWarmSum = rep(-1, length(pred_vals)),
                    voltinism = rep(NA, length(pred_vals)), 
                    Seas = rep(NA, length(pred_vals)), 
                    tempSeas = rep(NA, length(pred_vals)), 
                    dstdoy = rep(NA, length(pred_vals)),
                    overwinteringStage = rep(NA, length(pred_vals)), 
                    diurnality = rep(NA, length(pred_vals)),
                    id_cells = rep(NA, length(pred_vals)),
                    unique_name = rep(NA, length(pred_vals)))

#changing term 2.1 to have annual temp of 0
term2.2 <- data.frame(duration = rep(NA, length(pred_vals)),
                      annualTemp = rep(0, length(pred_vals)), 
                      annualPrec = rep(NA, length(pred_vals)),
                      UnusualColdSum = pred_vals, 
                      year = rep(NA, length(pred_vals)), 
                      UnusualWarmSum = rep(0, length(pred_vals)),
                      voltinism = rep(NA, length(pred_vals)), 
                      Seas = rep(NA, length(pred_vals)), 
                      tempSeas = rep(NA, length(pred_vals)), 
                      dstdoy = rep(NA, length(pred_vals)),
                      overwinteringStage = rep(NA, length(pred_vals)), 
                      diurnality = rep(NA, length(pred_vals)),
                      id_cells = rep(NA, length(pred_vals)),
                      unique_name = rep(NA, length(pred_vals)))

#changing term 3.1 to have annual temp of 0
term3.2 <- data.frame(duration = rep(NA, length(pred_vals)),
                      annualTemp = rep(0, length(pred_vals)), 
                      annualPrec = rep(NA, length(pred_vals)),
                      UnusualColdSum = pred_vals, 
                      year = rep(NA, length(pred_vals)), 
                      UnusualWarmSum = rep(1, length(pred_vals)),
                      voltinism = rep(NA, length(pred_vals)), 
                      Seas = rep(NA, length(pred_vals)), 
                      tempSeas = rep(NA, length(pred_vals)), 
                      dstdoy = rep(NA, length(pred_vals)),
                      overwinteringStage = rep(NA, length(pred_vals)), 
                      diurnality = rep(NA, length(pred_vals)),
                      id_cells = rep(NA, length(pred_vals)),
                      unique_name = rep(NA, length(pred_vals)))


# changing term 1.1 to have annual temp of 1
term1.3 <- data.frame(duration = rep(NA, length(pred_vals)),
                      annualTemp = rep(1, length(pred_vals)), 
                      annualPrec = rep(NA, length(pred_vals)),
                      UnusualColdSum = pred_vals, 
                      year = rep(NA, length(pred_vals)), 
                      UnusualWarmSum = rep(-1, length(pred_vals)),
                      voltinism = rep(NA, length(pred_vals)), 
                      Seas = rep(NA, length(pred_vals)), 
                      tempSeas = rep(NA, length(pred_vals)), 
                      dstdoy = rep(NA, length(pred_vals)),
                      overwinteringStage = rep(NA, length(pred_vals)), 
                      diurnality = rep(NA, length(pred_vals)),
                      id_cells = rep(NA, length(pred_vals)),
                      unique_name = rep(NA, length(pred_vals)))

#changing term 2.1 to have annual temp of 1
term2.3 <- data.frame(duration = rep(NA, length(pred_vals)),
                      annualTemp = rep(1, length(pred_vals)), 
                      annualPrec = rep(NA, length(pred_vals)),
                      UnusualColdSum = pred_vals, 
                      year = rep(NA, length(pred_vals)), 
                      UnusualWarmSum = rep(0, length(pred_vals)),
                      voltinism = rep(NA, length(pred_vals)), 
                      Seas = rep(NA, length(pred_vals)), 
                      tempSeas = rep(NA, length(pred_vals)), 
                      dstdoy = rep(NA, length(pred_vals)),
                      overwinteringStage = rep(NA, length(pred_vals)), 
                      diurnality = rep(NA, length(pred_vals)),
                      id_cells = rep(NA, length(pred_vals)),
                      unique_name = rep(NA, length(pred_vals)))

#changing term 3.1 to have annual temp of 1
term3.3 <- data.frame(duration = rep(NA, length(pred_vals)),
                      annualTemp = rep(1, length(pred_vals)), 
                      annualPrec = rep(NA, length(pred_vals)),
                      UnusualColdSum = pred_vals, 
                      year = rep(NA, length(pred_vals)), 
                      UnusualWarmSum = rep(1, length(pred_vals)),
                      voltinism = rep(NA, length(pred_vals)), 
                      Seas = rep(NA, length(pred_vals)), 
                      tempSeas = rep(NA, length(pred_vals)), 
                      dstdoy = rep(NA, length(pred_vals)),
                      overwinteringStage = rep(NA, length(pred_vals)), 
                      diurnality = rep(NA, length(pred_vals)),
                      id_cells = rep(NA, length(pred_vals)),
                      unique_name = rep(NA, length(pred_vals)))

pred.df.1.1 <- rbind(inla_mdf, term1.1)
pred.df.2.1 <- rbind(inla_mdf, term2.1)
pred.df.3.1 <- rbind(inla_mdf, term3.1)

pred.df.1.2 <- rbind(inla_mdf, term1.2)
pred.df.2.2 <- rbind(inla_mdf, term2.2)
pred.df.3.2 <- rbind(inla_mdf, term3.2)

pred.df.1.3 <- rbind(inla_mdf, term1.3)
pred.df.2.3 <- rbind(inla_mdf, term2.3)
pred.df.3.3 <- rbind(inla_mdf, term3.3)

pglmm1.1 <- pglmm(duration ~ annualTemp + annualPrec + UnusualColdSum + year + 
                  UnusualWarmSum + voltinism + Seas + tempSeas + dstdoy + 
                  overwinteringStage + diurnality + 
                  (1 | unique_name__) + (1 | id_cells) + 
                  annualTemp:UnusualColdSum + annualTemp:year + 
                  annualPrec:UnusualColdSum + annualPrec:UnusualWarmSum +
                  UnusualColdSum:UnusualWarmSum + year:voltinism,
                data = pred.df.1.1, 
                cov_ranef = list(unique_name = tt),
                bayes = T)

pglmm2.1 <- pglmm(duration ~ annualTemp + annualPrec + UnusualColdSum + year + 
                  UnusualWarmSum + voltinism + Seas + tempSeas + dstdoy + 
                  overwinteringStage + diurnality + 
                  (1 | unique_name__) + (1 | id_cells) + 
                  annualTemp:UnusualColdSum + annualTemp:year + 
                  annualPrec:UnusualColdSum + annualPrec:UnusualWarmSum +
                  UnusualColdSum:UnusualWarmSum + year:voltinism,
                data = pred.df.2.1, 
                cov_ranef = list(unique_name = tt),
                bayes = T)

pglmm3.1 <- pglmm(duration ~ annualTemp + annualPrec + UnusualColdSum + year + 
                  UnusualWarmSum + voltinism + Seas + tempSeas + dstdoy + 
                  overwinteringStage + diurnality + 
                  (1 | unique_name__) + (1 | id_cells) + 
                  annualTemp:UnusualColdSum + annualTemp:year + 
                  annualPrec:UnusualColdSum + annualPrec:UnusualWarmSum +
                  UnusualColdSum:UnusualWarmSum + year:voltinism,
                data = pred.df.3.1, 
                cov_ranef = list(unique_name = tt),
                bayes = T)

pglmm1.2 <- pglmm(duration ~ annualTemp + annualPrec + UnusualColdSum + year + 
                    UnusualWarmSum + voltinism + Seas + tempSeas + dstdoy + 
                    overwinteringStage + diurnality + 
                    (1 | unique_name__) + (1 | id_cells) + 
                    annualTemp:UnusualColdSum + annualTemp:year + 
                    annualPrec:UnusualColdSum + annualPrec:UnusualWarmSum +
                    UnusualColdSum:UnusualWarmSum + year:voltinism,
                  data = pred.df.1.2, 
                  cov_ranef = list(unique_name = tt),
                  bayes = T)

pglmm2.2 <- pglmm(duration ~ annualTemp + annualPrec + UnusualColdSum + year + 
                    UnusualWarmSum + voltinism + Seas + tempSeas + dstdoy + 
                    overwinteringStage + diurnality + 
                    (1 | unique_name__) + (1 | id_cells) + 
                    annualTemp:UnusualColdSum + annualTemp:year + 
                    annualPrec:UnusualColdSum + annualPrec:UnusualWarmSum +
                    UnusualColdSum:UnusualWarmSum + year:voltinism,
                  data = pred.df.2.2, 
                  cov_ranef = list(unique_name = tt),
                  bayes = T)

pglmm3.2 <- pglmm(duration ~ annualTemp + annualPrec + UnusualColdSum + year + 
                    UnusualWarmSum + voltinism + Seas + tempSeas + dstdoy + 
                    overwinteringStage + diurnality + 
                    (1 | unique_name__) + (1 | id_cells) + 
                    annualTemp:UnusualColdSum + annualTemp:year + 
                    annualPrec:UnusualColdSum + annualPrec:UnusualWarmSum +
                    UnusualColdSum:UnusualWarmSum + year:voltinism,
                  data = pred.df.3.2, 
                  cov_ranef = list(unique_name = tt),
                  bayes = T)

pglmm1.3 <- pglmm(duration ~ annualTemp + annualPrec + UnusualColdSum + year + 
                    UnusualWarmSum + voltinism + Seas + tempSeas + dstdoy + 
                    overwinteringStage + diurnality + 
                    (1 | unique_name__) + (1 | id_cells) + 
                    annualTemp:UnusualColdSum + annualTemp:year + 
                    annualPrec:UnusualColdSum + annualPrec:UnusualWarmSum +
                    UnusualColdSum:UnusualWarmSum + year:voltinism,
                  data = pred.df.1.3, 
                  cov_ranef = list(unique_name = tt),
                  bayes = T)

pglmm2.3 <- pglmm(duration ~ annualTemp + annualPrec + UnusualColdSum + year + 
                    UnusualWarmSum + voltinism + Seas + tempSeas + dstdoy + 
                    overwinteringStage + diurnality + 
                    (1 | unique_name__) + (1 | id_cells) + 
                    annualTemp:UnusualColdSum + annualTemp:year + 
                    annualPrec:UnusualColdSum + annualPrec:UnusualWarmSum +
                    UnusualColdSum:UnusualWarmSum + year:voltinism,
                  data = pred.df.2.3, 
                  cov_ranef = list(unique_name = tt),
                  bayes = T)

pglmm3.3 <- pglmm(duration ~ annualTemp + annualPrec + UnusualColdSum + year + 
                    UnusualWarmSum + voltinism + Seas + tempSeas + dstdoy + 
                    overwinteringStage + diurnality + 
                    (1 | unique_name__) + (1 | id_cells) + 
                    annualTemp:UnusualColdSum + annualTemp:year + 
                    annualPrec:UnusualColdSum + annualPrec:UnusualWarmSum +
                    UnusualColdSum:UnusualWarmSum + year:voltinism,
                  data = pred.df.3.3, 
                  cov_ranef = list(unique_name = tt),
                  bayes = T)

## make df of results
rdf1.1 <- pglmm1.1$inla.model$summary.linear.predictor[(nrow(mdf_phylo)+1):nrow(pred.df.1.1),] %>% 
  mutate(V1 = pred_vals,
         V2 = "Low",
         V3 = "Cold")

rdf2.1 <- pglmm2.1$inla.model$summary.linear.predictor[(nrow(mdf_phylo)+1):nrow(pred.df.1.1),] %>% 
  mutate(V1 = pred_vals,
         V2 = "Mid",
         V3 = "Cold")

rdf3.1 <- pglmm3.1$inla.model$summary.linear.predictor[(nrow(mdf_phylo)+1):nrow(pred.df.1.1),] %>% 
  mutate(V1 = pred_vals,
         V2 = "High",
         V3 = "Cold")

rdf1.2 <- pglmm1.2$inla.model$summary.linear.predictor[(nrow(mdf_phylo)+1):nrow(pred.df.1.1),] %>% 
  mutate(V1 = pred_vals,
         V2 = "Low",
         V3 = "Average")

rdf2.2 <- pglmm2.2$inla.model$summary.linear.predictor[(nrow(mdf_phylo)+1):nrow(pred.df.1.1),] %>% 
  mutate(V1 = pred_vals,
         V2 = "Mid",
         V3 = "Average")

rdf3.2 <- pglmm3.2$inla.model$summary.linear.predictor[(nrow(mdf_phylo)+1):nrow(pred.df.1.1),] %>% 
  mutate(V1 = pred_vals,
         V2 = "High",
         V3 = "Average")

rdf1.3 <- pglmm1.3$inla.model$summary.linear.predictor[(nrow(mdf_phylo)+1):nrow(pred.df.1.1),] %>% 
  mutate(V1 = pred_vals,
         V2 = "Low",
         V3 = "Warm")

rdf2.3 <- pglmm2.3$inla.model$summary.linear.predictor[(nrow(mdf_phylo)+1):nrow(pred.df.1.1),] %>% 
  mutate(V1 = pred_vals,
         V2 = "Mid",
         V3 = "Warm")

rdf3.3 <- pglmm3.3$inla.model$summary.linear.predictor[(nrow(mdf_phylo)+1):nrow(pred.df.1.1),] %>% 
  mutate(V1 = pred_vals,
         V2 = "High",
         V3 = "Warm")


rdf_total <- rbind(rdf1.1, rdf2.1, rdf3.1,
                   rdf1.2, rdf2.2, rdf3.2,
                   rdf1.3, rdf2.3, rdf3.3)

rdf_total$V2 <- factor(rdf_total$V2, levels = c("Low", "Mid", "High"))
rdf_total$V3 <- factor(rdf_total$V3, levels = c("Cold", "Average", "Warm"))

cold_warm_temp_plot <- ggplot(rdf_total, mapping = aes(x = V1, y = mean)) +
  geom_ribbon(mapping = aes(ymin = `0.025quant`, ymax = `0.975quant`,
                            fill = V2), alpha = 0.05 ) +
  geom_path(mapping = aes(y = `0.025quant`, color = V2), size = 0.25, linetype = 2) +
  geom_path(mapping = aes(y = `0.975quant`, color = V2), size = 0.25, linetype = 2) +
  geom_line(mapping = aes(color = V2), size = 1.05) +
  labs(x = "Cold day anomaly", y = "Duration", fill = "Warm day anomaly", color = "Warm day anomaly") + 
  scale_color_manual(values = c("#0DFFCA", "#47B4C7", "#2A4B61")) +
  scale_fill_manual(values = c("#0DFFCA", "#47B4C7", "#2A4B61")) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    panel.border = element_blank()
  ) +
  facet_wrap(~V3) 

cold_warm_temp_plot

ggsave(plot = cold_warm_temp_plot, filename = "Figures/dur_cold_warm_temp_interaction.png",
       height = 3, width = 7)

#save(dpt, file = "Fig_RData/dur_prec_temp.Rdata") # have not saved RData but this is what I do to save dataframe for plotting multipanel plots