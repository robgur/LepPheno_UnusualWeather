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


### Peak modeling
peak_m <- lmer(q50 ~ annualTemp + tempSeas + annualPrec + precSeas + year +
                  dstdoy + ndstcol + 
                  voltinism + overwinteringStage + diurnality + Seas +
                  annualTemp:voltinism +
                  annualTemp:overwinteringStage +
                  annualTemp:diurnality+
                  annualTemp:Seas +
                  year:voltinism +
                  year:overwinteringStage +
                  year:diurnality +
                  year:Seas +
                  annualTemp:annualPrec + 
                  annualTemp:year + 
                  tempSeas:year +
                  annualPrec:year +
                  precSeas:year +
                  (1|validName) + (1|id_cells),
                data = mdf,
                REML = FALSE, 
                lmerControl(optimizer = "bobyqa"),
                weights = peak_w)

peak_s <- step(peak_m)
peak_s

top_peak_m <- lmer(q50 ~ annualTemp + annualPrec +  year +
                     dstdoy + 
                     voltinism + overwinteringStage + Seas + 
                     annualTemp:overwinteringStage +
                     annualTemp:Seas +
                     year:overwinteringStage +
                     annualPrec:year +
                     precSeas:year +
                     (1|validName) + (1|id_cells),
                   data = mdf,
                   REML = FALSE, 
                   lmerControl(optimizer = "bobyqa"),
                   weights = peak_w)

car::vif(top_peak_m) # remove year

top_peak_m <- lmer(q50 ~ annualTemp + annualPrec +
                     dstdoy + 
                     voltinism + overwinteringStage + Seas + 
                     annualTemp:overwinteringStage +
                     annualTemp:Seas +
                     year:overwinteringStage +
                     annualPrec:year +
                     precSeas:year +
                     (1|validName) + (1|id_cells),
                   data = mdf,
                   REML = FALSE, 
                   lmerControl(optimizer = "bobyqa"),
                   weights = peak_w)

car::vif(top_peak_m) 

# check if model is stable
top_peak_m_s <- step(top_peak_m)
top_peak_m_s

### next we will make the no weight models
peak_m_noW <- lmer(q50 ~ annualTemp + tempSeas + annualPrec + precSeas + year +
                 dstdoy + ndstcol + 
                 voltinism + overwinteringStage + diurnality + Seas +
                 annualTemp:voltinism +
                 annualTemp:overwinteringStage +
                 annualTemp:diurnality+
                 annualTemp:Seas +
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

peak_s <- step(peak_m_noW)
peak_s

top_peak_noW <- lmer(q50 ~ annualTemp + annualPrec + year +
                       voltinism + overwinteringStage + Seas +
                       annualTemp:overwinteringStage +
                       annualTemp:Seas +
                       year:overwinteringStage +
                       annualTemp:annualPrec +
                       annualPrec:year +
                       precSeas:year +
                       (1|validName) + (1|id_cells),
                     data = mdf,
                     REML = FALSE, 
                     lmerControl(optimizer = "bobyqa"))

car::vif(top_peak_noW) # no year

top_peak_noW <- lmer(q50 ~ annualTemp + annualPrec + 
                       voltinism + overwinteringStage + Seas +
                       annualTemp:overwinteringStage +
                       annualTemp:Seas +
                       year:overwinteringStage +
                       annualTemp:annualPrec +
                       annualPrec:year +
                       precSeas:year +
                       (1|validName) + (1|id_cells),
                     data = mdf,
                     REML = FALSE, 
                     lmerControl(optimizer = "bobyqa"))

car::vif(top_peak_noW)

## check if top model is stable
top_peak_noW_s <- step(top_peak_noW)
top_peak_noW_s #yep all good

### finally the climate only modeling approach
peak_clim <- lmer(q50 ~ annualTemp + tempSeas + annualPrec + precSeas + 
                     dstdoy + ndstcol + 
                     annualTemp:annualPrec + 
                     (1|validName) + (1|id_cells),
                   data = mdf,
                   REML = FALSE, 
                   lmerControl(optimizer = "bobyqa"))

peak_clim_s <- step(peak_clim)
peak_clim_s

## now add the traits and year variables
peak_clim_traits <- lmer(q50 ~ annualTemp + annualPrec + 
                           annualTemp:annualPrec + 
                           voltinism + overwinteringStage + diurnality + Seas +
                           annualTemp:voltinism +
                           annualTemp:overwinteringStage +
                           annualTemp:diurnality+
                           annualTemp:Seas +
                           year:voltinism +
                           year:overwinteringStage +
                           year:diurnality +
                           annualTemp:year + 
                           annualPrec:year +
                           (1|validName) + (1|id_cells),
                         data = mdf,
                         REML = FALSE, 
                         lmerControl(optimizer = "bobyqa"))

peak_clim_traits_s <- step(peak_clim_traits)
peak_clim_traits_s

top_clim_traits_m <- lmer(q50 ~ annualTemp + annualPrec + 
                            voltinism + overwinteringStage + Seas +
                            annualTemp:annualPrec + 
                            annualTemp:overwinteringStage +
                            annualTemp:Seas +
                            overwinteringStage:year +
                            (1|validName) + (1|id_cells),
                          data = mdf,
                          REML = FALSE, 
                          lmerControl(optimizer = "bobyqa"))

car::vif(top_clim_traits_m)


AICc(top_clim_traits_m, top_peak_noW, top_peak_m)
Weights(AICc(top_clim_traits_m, top_peak_noW, top_peak_m))

## top peak model 
r.squaredGLMM(top_peak_noW)
summary(top_peak_noW)

plot_model(top_peak_noW, type = "pred", terms = c("annualTemp", "overwinteringStage"))
plot_model(top_peak_noW, type = "pred", terms = c("annualTemp", "Seas"))
plot_model(top_peak_noW, type = "pred", terms = c("year", "overwinteringStage"), ci.lvl = NA)
plot_model(top_peak_noW, type = "pred", terms = c("annualTemp", "annualPrec"))
plot_model(top_peak_noW, type = "pred", terms = c("year", "precSeas"), ci.lvl = NA)
plot_model(top_peak_noW, type = "pred", terms = c("year", "annualPrec"), ci.lvl = NA)


tab_model(top_peak_noW, file = "tables/middle_LMM.doc")

# phylogenetic model time
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

pglmm_peak <- pglmm(formula = q50 ~ annualTemp + annualPrec + 
                       voltinism + overwinteringStage + Seas +
                       annualTemp:overwinteringStage +
                       annualTemp:Seas +
                       year:overwinteringStage +
                       annualTemp:annualPrec +
                       annualPrec:year +
                       precSeas:year +
                       (1|unique_name__) + (1|id_cells),
                     data = mdf_phylo, 
                     cov_ranef = list(unique_name = tt), 
                     bayes = TRUE)

summary(pglmm_peak)
rr2::R2(pglmm_peak) #0.753
im <- pglmm_peak$inla.model
im_fix <- im$summary.fixed %>% 
  tibble::rownames_to_column("Effects") %>% 
  dplyr::select(Effects, mean, `0.025quant`, `0.975quant`)
tab_df(im_fix, title = "PGLMM Middle", file = "tables/middle_PGLMM.doc")


## Model Assumption Checks
# residuals
resids <- residuals(pglmm_peak)
hist(resids)
qqnorm(resids)
qqline(resids)
shapiro.test(resids)

# spatial autocorrelation
library(geoR)
jitter_dist = 125000 # 125 km
resid_peak <- mdf_phylo %>% 
  mutate(resid = resids)

gdf <- as.geodata(resid_peak, coords.col = c("lon", "lat"), data.col = "resid")
gdf <-  geoR::jitterDupCoords(gdf, max = jitter_dist)
v1 <- variog(gdf, trend = "1st")

plot(v1)

## Goodness of fit
rr2::R2(pglmm_peak) #0.746

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


peak_correlogram_plot <- ggplot() +
  geom_line(fit3_df,
            mapping = aes(x = distance, y = correlation)) +
  geom_point(fit3_df,
             mapping = aes(x = distance, y = correlation, fill = p.value),
             size = 3, shape = 21, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("white", "black")) +
  labs(x = "Distance", y = "Moran's I", fill = "P-Value") +
  theme_classic()

peak_correlogram_plot # no autocorrelation found at closest spatial lag, so spatial model not made

# save plot as Rdata to make multipanel Moran's I plot for SI
save(peak_correlogram_plot, file = "figures/peak_correlogram_plot.Rdata")
