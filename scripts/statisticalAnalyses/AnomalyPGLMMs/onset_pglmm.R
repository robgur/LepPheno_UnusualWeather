library(phyr)
library(car)
library(sjPlot)
library(ggplot2)
library(geoR)
library(ncf)
library(fpp2)

# get data ready for pglmm models
source('scripts/pglmm_dataSetUp.R')

# onset model
onset_m <- pglmm(onset ~ annualTemp + annualPrec + UnusualColdDaysOnset + 
                   UnusualWarmDaysOnset + voltinism + Seas + tempSeas +
                   dstdoy + overwinteringStage + (1 | unique_name__) + (1 | id_cells) +
                   annualTemp:annualPrec + annualTemp:UnusualColdDaysOnset + 
                   annualTemp:UnusualWarmDaysOnset +
                   annualTemp:voltinism + annualPrec:UnusualColdDaysOnset,
                data = mdf_phylo,
                cov_ranef = list(unique_name = tt),
               bayes = T)

summary(onset_m)
rr2::R2(onset_m) #0.780

phyr::plot_bayes(onset_m)

#save model results as table
im <- onset_m$inla.model
im_fix <- im$summary.fixed %>% 
  tibble::rownames_to_column("Effects") %>% 
  dplyr::select(Effects, mean, `0.025quant`, `0.975quant`)
tab_df(im_fix, title = "PGLMM Onset", file = "results/onset_PGLMM.doc")

## Model Assumption Checks
# residuals
resids <- residuals(onset_m)
hist(resids)
qqnorm(resids)
qqline(resids)

# spatial autocorrelation
# spatially visualize residuals
resid_onset <- mdf_phylo %>% 
  mutate(resid = resids)

ggplot(resid_onset, mapping = aes(x = lon, y = lat, color = resid)) +
  geom_point(size = 10) +
  scale_color_viridis_c(option = "turbo")

## look into temporal autocorrelation
resid_onset_sort <- resid_onset[order(resid_onset$year),]

acfPlot <- acf(resid_onset_sort$resid) 

on_acf <- ggAcf(resid_onset_sort$resid) +
  theme_bw() +
  ggtitle("")

save(on_acf, file = "Figures/onset_acf.rds")

ggplot(resid_onset, mapping = aes(x = year, y = resid)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(x = "Year", y = "Residual")+
  theme_bw()



# variogram approach
jitter_dist = 125000 # 125 km

gdf <- as.geodata(resid_onset, coords.col = c("lon", "lat"), data.col = "resid")
gdf <-  geoR::jitterDupCoords(gdf, max = jitter_dist)
v1 <- variog(gdf, trend = "1st")

plot(v1) # spatial dependence is greater at greater spatial lags

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

save(onset_correlogram_plot, file = "Figures/onset_correlogram.rds")


################################################################################
################################################################################
#############Interaction Figures ###############################################
################################################################################
################################################################################
# make interaction figures rob included in slack
inla_mdf <- mdf_phylo %>% 
  dplyr::select(onset, 
                annualTemp, 
                annualPrec, 
                UnusualColdDaysOnset, 
                year, 
                UnusualWarmDaysOnset, 
                voltinism, 
                Seas, 
                precSeas, 
                ndstcol,
                overwinteringStage, 
                diurnality,
                id_cells,
                unique_name)

pred_vals <- -3:3

#unusualcoldsum is x axis (so pred_vals),
#unusualwarmsum is legend, so start with -1
#voltinism is facet, start with U
term1.U <- data.frame(onset = rep(NA, length(pred_vals)),
                      annualTemp = rep(-1, length(pred_vals)), 
                      annualPrec = rep(NA, length(pred_vals)),
                      UnusualColdDaysOnset = rep(NA, length(pred_vals)), 
                      year = rep(NA, length(pred_vals)), 
                      UnusualWarmDaysOnset = pred_vals,
                      voltinism = rep("U", length(pred_vals)), 
                      Seas = rep(NA, length(pred_vals)), 
                      precSeas = rep(NA, length(pred_vals)), 
                      ndstcol = rep(NA, length(pred_vals)),
                      overwinteringStage = rep(NA, length(pred_vals)), 
                      diurnality = rep(NA, length(pred_vals)),
                      id_cells = rep(NA, length(pred_vals)),
                      unique_name = rep(NA, length(pred_vals)))

#unusualcoldsum is x axis (so pred_vals),
#unusualwarmsum is legend, now we on to 0
#voltinism is facet so keep at "U
term2.U <- data.frame(onset = rep(NA, length(pred_vals)),
                      annualTemp = rep(0, length(pred_vals)), 
                      annualPrec = rep(NA, length(pred_vals)),
                      UnusualColdDaysOnset = rep(NA, length(pred_vals)), 
                      year = rep(NA, length(pred_vals)), 
                      UnusualWarmDaysOnset = pred_vals,
                      voltinism = rep("U", length(pred_vals)), 
                      Seas = rep(NA, length(pred_vals)), 
                      precSeas = rep(NA, length(pred_vals)), 
                      ndstcol = rep(NA, length(pred_vals)),
                      overwinteringStage = rep(NA, length(pred_vals)), 
                      diurnality = rep(NA, length(pred_vals)),
                      id_cells = rep(NA, length(pred_vals)),
                      unique_name = rep(NA, length(pred_vals)))

#unusualcoldsum is x axis (so pred_vals),
#unusualwarmsum is legend, now we on to 1
#voltinism is facet so keep at "U
term3.U <- data.frame(onset = rep(NA, length(pred_vals)),
                      annualTemp = rep(1, length(pred_vals)), 
                      annualPrec = rep(NA, length(pred_vals)),
                      UnusualColdDaysOnset = rep(NA, length(pred_vals)), 
                      year = rep(NA, length(pred_vals)), 
                      UnusualWarmDaysOnset =pred_vals,
                      voltinism = rep("U", length(pred_vals)), 
                      Seas = rep(NA, length(pred_vals)), 
                      precSeas = rep(NA, length(pred_vals)), 
                      ndstcol = rep(NA, length(pred_vals)),
                      overwinteringStage = rep(NA, length(pred_vals)), 
                      diurnality = rep(NA, length(pred_vals)),
                      id_cells = rep(NA, length(pred_vals)),
                      unique_name = rep(NA, length(pred_vals)))

## now change facet to "M"
term1.M <- data.frame(onset = rep(NA, length(pred_vals)),
                      annualTemp = rep(-1, length(pred_vals)), 
                      annualPrec = rep(NA, length(pred_vals)),
                      UnusualColdDaysOnset = rep(NA, length(pred_vals)), 
                      year = rep(NA, length(pred_vals)), 
                      UnusualWarmDaysOnset =pred_vals,
                      voltinism = rep("M", length(pred_vals)), 
                      Seas = rep(NA, length(pred_vals)), 
                      precSeas = rep(NA, length(pred_vals)), 
                      ndstcol = rep(NA, length(pred_vals)),
                      overwinteringStage = rep(NA, length(pred_vals)), 
                      diurnality = rep(NA, length(pred_vals)),
                      id_cells = rep(NA, length(pred_vals)),
                      unique_name = rep(NA, length(pred_vals)))

#M
term2.M <- data.frame(onset = rep(NA, length(pred_vals)),
                      annualTemp = rep(0, length(pred_vals)), 
                      annualPrec = rep(NA, length(pred_vals)),
                      UnusualColdDaysOnset = rep(NA, length(pred_vals)), 
                      year = rep(NA, length(pred_vals)), 
                      UnusualWarmDaysOnset =pred_vals,
                      voltinism = rep("M", length(pred_vals)), 
                      Seas = rep(NA, length(pred_vals)), 
                      precSeas = rep(NA, length(pred_vals)), 
                      ndstcol = rep(NA, length(pred_vals)),
                      overwinteringStage = rep(NA, length(pred_vals)), 
                      diurnality = rep(NA, length(pred_vals)),
                      id_cells = rep(NA, length(pred_vals)),
                      unique_name = rep(NA, length(pred_vals)))

# M
term3.M <- data.frame(onset = rep(NA, length(pred_vals)),
                      annualTemp = rep(1, length(pred_vals)), 
                      annualPrec = rep(NA, length(pred_vals)),
                      UnusualColdDaysOnset = rep(NA, length(pred_vals)), 
                      year = rep(NA, length(pred_vals)), 
                      UnusualWarmDaysOnset =pred_vals,
                      voltinism = rep("M", length(pred_vals)), 
                      Seas = rep(NA, length(pred_vals)), 
                      precSeas = rep(NA, length(pred_vals)), 
                      ndstcol = rep(NA, length(pred_vals)),
                      overwinteringStage = rep(NA, length(pred_vals)), 
                      diurnality = rep(NA, length(pred_vals)),
                      id_cells = rep(NA, length(pred_vals)),
                      unique_name = rep(NA, length(pred_vals)))

pred.df.1.U <- rbind(inla_mdf, term1.U)
pred.df.2.U <- rbind(inla_mdf, term2.U)
pred.df.3.U <- rbind(inla_mdf, term3.U)

pred.df.1.M <- rbind(inla_mdf, term1.M)
pred.df.2.M <- rbind(inla_mdf, term2.M)
pred.df.3.M <- rbind(inla_mdf, term3.M)


# fit models
pglmm1.U <- pglmm(onset ~ annualTemp + annualPrec + UnusualColdDaysOnset +
                    year + UnusualWarmDaysOnset + voltinism + Seas + precSeas + 
                    ndstcol + overwinteringStage + 
                    (1 | unique_name__) + (1 | id_cells) + 
                    annualTemp:annualPrec + annualTemp:UnusualColdDaysOnset + 
                    annualTemp:year + annualTemp:UnusualWarmDaysOnset + 
                    annualTemp:voltinism + annualPrec:UnusualColdDaysOnset, 
                  data = pred.df.1.U,
                  cov_ranef = list(unique_name = tt),
                  bayes = T)

pglmm2.U <- pglmm(onset ~ annualTemp + annualPrec + UnusualColdDaysOnset +
                    year + UnusualWarmDaysOnset + voltinism + Seas + precSeas + 
                    ndstcol + overwinteringStage + 
                    (1 | unique_name__) + (1 | id_cells) + 
                    annualTemp:annualPrec + annualTemp:UnusualColdDaysOnset + 
                    annualTemp:year + annualTemp:UnusualWarmDaysOnset + 
                    annualTemp:voltinism + annualPrec:UnusualColdDaysOnset, 
                  data = pred.df.2.U,
                  cov_ranef = list(unique_name = tt),
                  bayes = T)

pglmm3.U <-pglmm(onset ~ annualTemp + annualPrec + UnusualColdDaysOnset +
                   year + UnusualWarmDaysOnset + voltinism + Seas + precSeas + 
                   ndstcol + overwinteringStage + 
                   (1 | unique_name__) + (1 | id_cells) + 
                   annualTemp:annualPrec + annualTemp:UnusualColdDaysOnset + 
                   annualTemp:year + annualTemp:UnusualWarmDaysOnset + 
                   annualTemp:voltinism + annualPrec:UnusualColdDaysOnset, 
                 data = pred.df.3.U,
                 cov_ranef = list(unique_name = tt),
                 bayes = T)

pglmm1.M <- pglmm(onset ~ annualTemp + annualPrec + UnusualColdDaysOnset +
                    year + UnusualWarmDaysOnset + voltinism + Seas + precSeas + 
                    ndstcol + overwinteringStage + 
                    (1 | unique_name__) + (1 | id_cells) + 
                    annualTemp:annualPrec + annualTemp:UnusualColdDaysOnset + 
                    annualTemp:year + annualTemp:UnusualWarmDaysOnset + 
                    annualTemp:voltinism + annualPrec:UnusualColdDaysOnset, 
                  data = pred.df.1.M,
                  cov_ranef = list(unique_name = tt),
                  bayes = T)

pglmm2.M <- pglmm(onset ~ annualTemp + annualPrec + UnusualColdDaysOnset +
                    year + UnusualWarmDaysOnset + voltinism + Seas + precSeas + 
                    ndstcol + overwinteringStage + 
                    (1 | unique_name__) + (1 | id_cells) + 
                    annualTemp:annualPrec + annualTemp:UnusualColdDaysOnset + 
                    annualTemp:year + annualTemp:UnusualWarmDaysOnset + 
                    annualTemp:voltinism + annualPrec:UnusualColdDaysOnset, 
                  data = pred.df.2.M,
                  cov_ranef = list(unique_name = tt),
                  bayes = T)

pglmm3.M <- pglmm(onset ~ annualTemp + annualPrec + UnusualColdDaysOnset +
                    year + UnusualWarmDaysOnset + voltinism + Seas + precSeas + 
                    ndstcol + overwinteringStage + 
                    (1 | unique_name__) + (1 | id_cells) + 
                    annualTemp:annualPrec + annualTemp:UnusualColdDaysOnset + 
                    annualTemp:year + annualTemp:UnusualWarmDaysOnset + 
                    annualTemp:voltinism + annualPrec:UnusualColdDaysOnset, 
                  data = pred.df.3.M,
                  cov_ranef = list(unique_name = tt),
                  bayes = T)

## make df of results
rdf1.1 <- pglmm1.U$inla.model$summary.linear.predictor[(nrow(mdf_phylo)+1):nrow(pred.df.1.U),] %>% 
  mutate(V1 = pred_vals,
         V2 = "Cold",
         V3 = "Univoltine")

rdf2.1 <- pglmm2.U$inla.model$summary.linear.predictor[(nrow(mdf_phylo)+1):nrow(pred.df.1.U),] %>% 
  mutate(V1 = pred_vals,
         V2 = "Average",
         V3 = "Univoltine")

rdf3.1 <- pglmm3.U$inla.model$summary.linear.predictor[(nrow(mdf_phylo)+1):nrow(pred.df.1.U),] %>% 
  mutate(V1 = pred_vals,
         V2 = "Warm",
         V3 = "Univoltine")

rdf1.2 <- pglmm1.M$inla.model$summary.linear.predictor[(nrow(mdf_phylo)+1):nrow(pred.df.1.M),] %>% 
  mutate(V1 = pred_vals,
         V2 = "Cold",
         V3 = "Multivoltine")

rdf2.2 <- pglmm2.M$inla.model$summary.linear.predictor[(nrow(mdf_phylo)+1):nrow(pred.df.1.M),] %>% 
  mutate(V1 = pred_vals,
         V2 = "Average",
         V3 = "Multivoltine")

rdf3.2 <- pglmm3.M$inla.model$summary.linear.predictor[(nrow(mdf_phylo)+1):nrow(pred.df.1.M),] %>% 
  mutate(V1 = pred_vals,
         V2 = "Warm",
         V3 = "Multivoltine")


rdf_total <- rbind(rdf1.1, rdf2.1, rdf3.1,
                   rdf1.2, rdf2.2, rdf3.2)

rdf_total$V2 <- factor(rdf_total$V2, levels = c("Cold", "Average", "Warm"))
rdf_total$V3 <- factor(rdf_total$V3, levels = c("Univoltine", "Multivoltine"))

warm_temp_volt_plot <- ggplot(rdf_total, mapping = aes(x = V1, y = mean)) +
  geom_ribbon(mapping = aes(ymin = `0.025quant`, ymax = `0.975quant`,
                            fill = V2), alpha = 0.05 ) +
  geom_path(mapping = aes(y = `0.025quant`, color = V2), size = 0.25, linetype = 2) +
  geom_path(mapping = aes(y = `0.975quant`, color = V2), size = 0.25, linetype = 2) +
  geom_line(mapping = aes(color = V2), size = 1.05) +
  labs(x = "Warm day anomaly", y = "Onset", fill = "Annual temp", color = "Annual temp") + 
  scale_color_manual(values = c("#046865", "#FDCA9B", "#8C032C")) +
  scale_fill_manual(values = c("#046865", "#FDCA9B", "#8C032C")) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    panel.border = element_blank()
  ) +
  facet_wrap(~V3) 

warm_temp_volt_plot

ggsave(plot = warm_temp_volt_plot, filename = "Figures/onset_warm_temp_volt_interaction.png",
       height = 3, width = 7)
