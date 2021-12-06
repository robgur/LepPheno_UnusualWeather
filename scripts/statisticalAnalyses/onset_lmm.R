library(tidyverse)
library(lme4)
library(lmerTest)
library(car)
library(sjPlot)
library(MuMIn)

## read in data
mdf <- read.csv("data/LMM_Data/mdf_pruned_hopkins.csv") %>% 
  dplyr::filter(overwinteringStage != "None")

mdf <- mdf %>% 
  mutate(annualTemp = scale(annualTemp),
         tempSeas = scale(tempSeas),
         annualPrec = scale(annualPrec),
         precSeas = scale(precSeas),
         year = scale(year),
         dstdoy = scale(dstdoy),
         ndstcol = scale(ndstcol)) %>% 
  mutate(dur = q95 - q5,
         dur_CI = mean(c(q5_CI, q95_CI))) %>% 
  mutate(on_w = 1/(q5_CI + 1),
         off_w = 1/(q95_CI + 1),
         peak_w = 1/(q95_CI +1),
         dur_w = 1/(dur_CI + 1)) %>% 
  mutate(voltinism = if_else(condition = voltinism == "F",
                             true = "M",
                             false = voltinism))


### Onset modeling
onset_m <- lmer(q5 ~ annualTemp + tempSeas + annualPrec + precSeas + year +
                  dstdoy + ndstcol + 
                  voltinism + overwinteringStage + diurnality + 
                  annualTemp:voltinism +
                  annualTemp:overwinteringStage +
                  annualTemp:diurnality+
                  year:voltinism +
                  year:overwinteringStage +
                  year:diurnality +
                  annualTemp:annualPrec + 
                  annualTemp:year + 
                  tempSeas:year +
                  annualPrec:year +
                  precSeas:year +
                  (1|validName) + (1|id_cells),
                data = mdf,
                REML = FALSE, 
                lmerControl(optimizer = "bobyqa"),
                weights = on_w)

onset_s <- step(onset_m)
onset_s

onset_top <- lmer(q5 ~ annualTemp + tempSeas + annualPrec + precSeas +
                    year + dstdoy +
                    voltinism + overwinteringStage + diurnality + 
                    (1 | validName) + (1 | id_cells) + 
                    annualTemp:voltinism +
                    annualTemp:diurnality + 
                    year:overwinteringStage +
                    year:diurnality +
                    annualTemp:annualPrec +
                    annualTemp:year +
                    tempSeas:year +
                    annualPrec:year +
                    precSeas:year,
                  data = mdf,
                  REML = FALSE, 
                  lmerControl(optimizer = "bobyqa"),
                  weights = on_w)

vif(onset_top) # vif too high on year effect, removing,
# diurnality:year

onset_top <- lmer(q5 ~ annualTemp + tempSeas + annualPrec + precSeas +
                     dstdoy +
                    voltinism + overwinteringStage + diurnality + 
                    (1 | validName) + (1 | id_cells) + 
                    annualTemp:voltinism +
                    annualTemp:diurnality + 
                    year:overwinteringStage +
                    annualTemp:annualPrec +
                    annualTemp:year +
                    tempSeas:year +
                    annualPrec:year +
                    precSeas:year,
                  data = mdf,
                  REML = FALSE, 
                  lmerControl(optimizer = "bobyqa"),
                  weights = on_w)
car::vif(onset_top) ## still too high for annualTemp:diurnality
onset_top <- lmer(q5 ~ annualTemp + tempSeas + annualPrec + precSeas +
                    dstdoy +
                    voltinism + overwinteringStage + diurnality + 
                    (1 | validName) + (1 | id_cells) + 
                    annualTemp:voltinism +
                    year:overwinteringStage +
                    annualTemp:annualPrec +
                    annualTemp:year +
                    tempSeas:year +
                    annualPrec:year +
                    precSeas:year,
                  data = mdf,
                  REML = FALSE, 
                  lmerControl(optimizer = "bobyqa"),
                  weights = on_w)

vif(onset_top)

summary(onset_top)

## now do this without the weights
onset_m <- lmer(q5 ~ annualTemp + tempSeas + annualPrec + precSeas + year +
                  dstdoy + ndstcol + 
                  voltinism + overwinteringStage + diurnality + 
                  annualTemp:voltinism +
                  annualTemp:overwinteringStage +
                  annualTemp:diurnality+
                  year:voltinism +
                  year:overwinteringStage +
                  year:diurnality +
                  annualTemp:annualPrec + 
                  annualTemp:year + 
                  tempSeas:year +
                  annualPrec:year +
                  precSeas:year +
                  (1|validName) + (1|id_cells),
                data = mdf,
                REML = FALSE, 
                lmerControl(optimizer = "bobyqa"))

onset_s <- step(onset_m)
onset_s

onset_top_noW <- lmer(q5 ~ annualTemp + annualPrec + precSeas + year +
                        dstdoy +
                        overwinteringStage + 
                        annualTemp:overwinteringStage +
                        year:overwinteringStage +
                        annualTemp:annualPrec + 
                        annualPrec:year +
                        (1|validName) + (1|id_cells),
                      data = mdf,
                      REML = FALSE, 
                      lmerControl(optimizer = "bobyqa"))

car::vif(onset_top_noW)


## remove year
onset_top_noW <- lmer(q5 ~ annualTemp + annualPrec + precSeas +
                        dstdoy +
                        overwinteringStage + 
                        annualTemp:overwinteringStage +
                        year:overwinteringStage +
                        annualTemp:annualPrec + 
                        annualPrec:year +
                        (1|validName) + (1|id_cells),
                      data = mdf,
                      REML = FALSE, 
                      lmerControl(optimizer = "bobyqa"))

car::vif(onset_top_noW)
r.squaredGLMM(onset_top_noW)

## onset doing climate models first
onset_clim <- lmer(q5 ~ annualTemp + tempSeas + annualPrec + precSeas + 
                     dstdoy + ndstcol + 
                     annualTemp:annualPrec + 
                     (1|validName) + (1|id_cells),
                   data = mdf,
                   REML = FALSE, 
                   lmerControl(optimizer = "bobyqa"))

onset_s <- step(onset_clim)
onset_s

# now add in more traits
model_climFirst <- lmer(q5 ~ annualTemp + annualPrec + 
                          dstdoy + 
                          annualTemp:annualPrec + 
                          year +
                          voltinism + overwinteringStage + diurnality + 
                          annualTemp:voltinism +
                          annualTemp:overwinteringStage +
                          annualTemp:diurnality+
                          year:voltinism +
                          year:overwinteringStage +
                          year:diurnality +
                          annualTemp:year + 
                          annualPrec:year +
                          (1|validName) + (1|id_cells),
                        data = mdf,
                        REML = FALSE, 
                        lmerControl(optimizer = "bobyqa"))


model_climFirst_s <- step(model_climFirst)
model_climFirst_s

top_model_climFirst <- lmer(q5 ~ annualTemp + annualPrec + 
                              dstdoy + 
                              year +
                              overwinteringStage + 
                              annualTemp:annualPrec + 
                              annualTemp:overwinteringStage +
                              year:overwinteringStage +
                              annualPrec:year +
                              (1|validName) + (1|id_cells),
                            data = mdf,
                            REML = FALSE, 
                            lmerControl(optimizer = "bobyqa"))

car::vif(top_model_climFirst) # remove year

top_model_climFirst <- lmer(q5 ~ annualTemp + annualPrec + 
                              dstdoy + 
                              overwinteringStage + 
                              annualTemp:annualPrec + 
                              annualTemp:overwinteringStage +
                              year:overwinteringStage +
                              annualPrec:year +
                              (1|validName) + (1|id_cells),
                            data = mdf,
                            REML = FALSE, 
                            lmerControl(optimizer = "bobyqa"))

car::vif(top_model_climFirst)

AICc(top_model_climFirst, onset_top_noW, onset_top)
Weights(AICc(top_model_climFirst, onset_top_noW, onset_top))

## top onset model 
summary(onset_top_noW)

plot_model(onset_top_noW, type = "pred", terms = c("annualTemp", "overwinteringStage"))
plot_model(onset_top_noW, type = "pred", terms = c("year", "overwinteringStage"), ci.lvl = NA)
plot_model(onset_top_noW, type = "pred", terms = c("annualTemp", "annualPrec"))
plot_model(onset_top_noW, type = "pred", terms = c("year", "annualPrec"),ci.lvl = NA)

tab_model(onset_top_noW, file = "tables/onset_LMM.doc")

## phylogenetic model time
library(phyr)
library(ape)
# function to capitalize spp names
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

tt <- read.tree("data/phylogeny/insect_tree_wBranches.tre")
tt$tip.label <- stringr::str_replace(tt$tip.label, pattern = "_", " ")
tt$tip.label <- firstup(tt$tip.label)
phylo_list <- read.csv('data/phylogeny/phy_spp_names.csv')

phylo_list <- phylo_list %>% 
  mutate(search_string2 = firstup(phylo_list$search_string)) %>% 
  dplyr::select(search_string, unique_name, ott_id, search_string2)

mdf_phylo <- left_join(mdf, phylo_list, by = c("validName" = "search_string2"))

mdf_phylo <- mdf_phylo %>% 
  filter(validName != "Polites sonora") %>% 
  filter(validName != "Nadata gibossa")

pglmm_onset <- pglmm(formula = q5 ~ annualTemp + annualPrec + precSeas +
                       dstdoy +
                       overwinteringStage + 
                       annualTemp:overwinteringStage +
                       year:overwinteringStage +
                       annualTemp:annualPrec + 
                       annualPrec:year +
                       (1|unique_name__) + (1|id_cells),
              data = mdf_phylo, 
              cov_ranef = list(unique_name = tt), 
              bayes = TRUE)

summary(onset_top_noW)
summary(pglmm_onset)

im <- pglmm_onset$inla.model
im_fix <- im$summary.fixed %>% 
  tibble::rownames_to_column("Effects") %>% 
  dplyr::select(Effects, mean, `0.025quant`, `0.975quant`)
tab_df(im_fix, title = "PGLMM Onset", file = "tables/onset_PGLMM.doc")

## Model Assumption Checks
# residuals
resids <- residuals(pglmm_onset)
hist(resids)
qqnorm(resids)
qqline(resids)
shapiro.test(resids)

# spatial autocorrelation
library(geoR)
jitter_dist = 125000 # 125 km
resid_onset <- mdf_phylo %>% 
  mutate(resid = resids)

gdf <- as.geodata(resid_onset, coords.col = c("lon", "lat"), data.col = "resid")
gdf <-  geoR::jitterDupCoords(gdf, max = jitter_dist)
v1 <- variog(gdf, trend = "1st")

plot(v1)

## Goodness of fit
rr2::R2(pglmm_onset) #0.746

## Examine spatial autocorrelation a second way
library(ncf)
rdf <- mutate(mdf_phylo, residuals = resids)

fit3 = correlog(x = rdf$lon, y = rdf$lat, z = rdf$residuals,
                increment = 250000, latlon = F, resamp = 100)
plot(fit3)

fit3_df <- data.frame(distance = fit3$mean.of.class, 
                      correlation = fit3$correlation,
                      p.value = if_else(fit3$p < 0.025, 
                                        true = "Sig",
                                        false = "Not Sig")) %>% 
  filter(distance < 5000000 & distance > 1)


onset_correlogram_plot <- ggplot() +
  geom_line(fit3_df,
            mapping = aes(x = distance, y = correlation)) +
  geom_point(fit3_df,
             mapping = aes(x = distance, y = correlation, fill = p.value),
             size = 3, shape = 21, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("white", "black")) +
  labs(x = "Distance", y = "Moran's I", fill = "P-Value") +
  theme_classic()

onset_correlogram_plot # no autocorrelation found at closest spatial lag, so spatial model not made

# save plot as Rdata to make multipanel Moran's I plot for SI
save(onset_correlogram_plot, file = "figures/onset_correlogram_plot.Rdata")
