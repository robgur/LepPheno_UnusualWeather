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


### Duration modeling
dur_m <- lmer(dur ~ annualTemp + tempSeas + annualPrec + precSeas + year +
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
              weights = dur_w)

dur_m_s <- step(dur_m)
dur_m_s

top_dur <- lmer(dur ~ annualTemp + annualPrec + year +
                  dstdoy +
                  voltinism + overwinteringStage + diurnality + Seas +
                  annualTemp:voltinism +
                  annualTemp:overwinteringStage +
                  annualTemp:Seas +
                  year:diurnality + 
                  year:Seas +
                  annualTemp:annualPrec +
                  (1|validName) + (1|id_cells),
                data = mdf,
                REML = FALSE, 
                lmerControl(optimizer = "bobyqa"),
                weights = dur_w)

car::vif(top_dur) # remove year

top_dur <- lmer(dur ~ annualTemp + annualPrec + 
                  dstdoy +
                  voltinism + overwinteringStage + diurnality + Seas +
                  annualTemp:voltinism +
                  annualTemp:overwinteringStage +
                  annualTemp:Seas +
                  year:diurnality + 
                  year:Seas +
                  annualTemp:annualPrec +
                  (1|validName) + (1|id_cells),
                data = mdf,
                REML = FALSE, 
                lmerControl(optimizer = "bobyqa"),
                weights = dur_w)

car::vif(top_dur) ## still too high remove annualTemp:Seas

top_dur <- lmer(dur ~ annualTemp + annualPrec + 
                  dstdoy +
                  voltinism + overwinteringStage + diurnality + Seas +
                  annualTemp:voltinism +
                  annualTemp:overwinteringStage +
                  year:diurnality + 
                  year:Seas +
                  annualTemp:annualPrec +
                  (1|validName) + (1|id_cells),
                data = mdf,
                REML = FALSE, 
                lmerControl(optimizer = "bobyqa"),
                weights = dur_w)

car::vif(top_dur) ## still too high remove annualTemp:Seas
#stable top model?
top_dur_s <- step(top_dur)
top_dur_s ## nope new top model

top_dur <- lmer(dur ~ annualTemp + annualPrec + 
                  dstdoy +
                  voltinism + overwinteringStage + diurnality + Seas +
                  annualTemp:voltinism +
                  annualTemp:overwinteringStage +
                  year:Seas +
                  annualTemp:annualPrec +
                  (1|validName) + (1|id_cells),
                data = mdf,
                REML = FALSE, 
                lmerControl(optimizer = "bobyqa"),
                weights = dur_w)

vif(top_dur)
# now stable?
top_dur_s2 <- step(top_dur)
top_dur_s2 # yep good

## no weights models
dur_m_noW <- lmer(dur ~ annualTemp + tempSeas + annualPrec + precSeas + year +
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
              lmerControl(optimizer = "bobyqa"))

dur_s_noW <- step(dur_m_noW)
dur_s_noW

dur_top_noW <- lmer(dur ~ annualTemp + annualPrec + year +
                      dstdoy +
                      voltinism + overwinteringStage + diurnality + Seas +
                      annualTemp:voltinism +
                      annualTemp:overwinteringStage +
                      annualTemp:Seas + 
                      year:diurnality +
                      year:Seas +
                      annualTemp:annualPrec +
                      (1|validName) + (1|id_cells),
                    data = mdf,
                    REML = FALSE, 
                    lmerControl(optimizer = "bobyqa"))

car::vif(dur_top_noW) # remove year

dur_top_noW <- lmer(dur ~ annualTemp + annualPrec +
                      dstdoy +
                      voltinism + overwinteringStage + diurnality + Seas +
                      annualTemp:voltinism +
                      annualTemp:overwinteringStage +
                      annualTemp:Seas + 
                      year:diurnality +
                      year:Seas +
                      annualTemp:annualPrec +
                      (1|validName) + (1|id_cells),
                    data = mdf,
                    REML = FALSE, 
                    lmerControl(optimizer = "bobyqa"))

car::vif(dur_top_noW) # remove annualTemp:Seas

dur_top_noW <- lmer(dur ~ annualTemp + annualPrec +
                      dstdoy +
                      voltinism + overwinteringStage + diurnality + Seas +
                      annualTemp:voltinism +
                      annualTemp:overwinteringStage +
                      year:diurnality +
                      year:Seas +
                      annualTemp:annualPrec +
                      (1|validName) + (1|id_cells),
                    data = mdf,
                    REML = FALSE, 
                    lmerControl(optimizer = "bobyqa"))

car::vif(dur_top_noW)
# stable top model?
dur_top_noW_s <- step(dur_top_noW)
dur_top_noW_s # nope, new top model

dur_top_noW <- lmer(dur ~ annualTemp + annualPrec +
                      dstdoy +
                      voltinism + overwinteringStage + diurnality + Seas +
                      annualTemp:voltinism +
                      annualTemp:overwinteringStage +
                      year:Seas +
                      annualTemp:annualPrec +
                      (1|validName) + (1|id_cells),
                    data = mdf,
                    REML = FALSE, 
                    lmerControl(optimizer = "bobyqa"))

vif(dur_top_noW)
# now stable?
dur_top_noW_s2 <- step(dur_top_noW)
dur_top_noW_s2

# Clim first models now
## climate first model
dur_clim <- lmer(dur ~ annualTemp + tempSeas + annualPrec + precSeas + 
                   dstdoy + ndstcol + 
                   annualTemp:annualPrec + 
                   (1|validName) + (1|id_cells),
                 data = mdf,
                 REML = FALSE, 
                 lmerControl(optimizer = "bobyqa"))

dur_clim_s <- step(dur_clim)
dur_clim_s

## add in traits now
dur_clim_traits <- lmer(dur ~ annualTemp + annualPrec + 
                   dstdoy +
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
                     annualPrec:year +
                   (1|validName) + (1|id_cells),
                 data = mdf,
                 REML = FALSE, 
                 lmerControl(optimizer = "bobyqa"))

dur_clim_traits_s <- step(dur_clim_traits)
dur_clim_traits_s

dur_clim_top <- lmer(dur ~ annualTemp + annualPrec +
                          dstdoy +
                          voltinism + overwinteringStage + diurnality + Seas +
                          annualTemp:voltinism +
                          annualTemp:overwinteringStage +
                          annualTemp:Seas +
                          Seas:year +
                          annualTemp:annualPrec +
                          (1|validName) + (1|id_cells),
                        data = mdf,
                        REML = FALSE, 
                        lmerControl(optimizer = "bobyqa"))

vif(dur_clim_top)
#stable top model?
dur_clim_top_s <- step(dur_clim_top)
dur_clim_top_s # yep, good

AICc(dur_clim_top, dur_top_noW, top_dur) 
Weights(AICc(dur_clim_top, dur_top_noW, top_dur))

summary(dur_clim_top)

plot_model(dur_clim_top, type = "pred", terms = c("annualTemp", "voltinism"))
plot_model(dur_clim_top, type = "pred", terms = c("annualTemp", "overwinteringStage"))
plot_model(dur_clim_top, type = "pred", terms = c("annualTemp", "Seas"))
plot_model(dur_clim_top, type = "pred", terms = c("year", "Seas"), ci.lvl = NA)
plot_model(dur_clim_top, type = "pred", terms = c("annualTemp", "annualPrec"))

tab_model(dur_clim_top, file = "tables/duration_LMM.doc")

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

pglmm_dur <- pglmm(formula = dur ~ annualTemp + annualPrec +
                                   dstdoy +
                                   voltinism + overwinteringStage + diurnality + Seas +
                                   annualTemp:voltinism +
                                   annualTemp:overwinteringStage +
                                   annualTemp:Seas +
                                   Seas:year +
                                   annualTemp:annualPrec +
                                   (1|unique_name__) + (1|id_cells),
                      data = mdf_phylo, 
                      cov_ranef = list(unique_name = tt), 
                      bayes = TRUE)

summary(pglmm_dur)

rr2::R2(pglmm_dur) #0.72
im <- pglmm_dur$inla.model
im_fix <- im$summary.fixed %>% 
  tibble::rownames_to_column("Effects") %>% 
  dplyr::select(Effects, mean, `0.025quant`, `0.975quant`)
tab_df(im_fix, title = "PGLMM Duration", file = "tables/duration_PGLMM.doc")

## Model Assumption Checks
# residuals
resids <- residuals(pglmm_dur)
hist(resids)
qqnorm(resids)
qqline(resids)
shapiro.test(resids)

# spatial autocorrelation
library(geoR)
jitter_dist = 125000 # 125 km
resid_dur <- mdf_phylo %>% 
  mutate(resid = resids)

gdf <- as.geodata(resid_dur, coords.col = c("lon", "lat"), data.col = "resid")
gdf <-  geoR::jitterDupCoords(gdf, max = jitter_dist)
v1 <- variog(gdf, trend = "1st")

plot(v1)

## Goodness of fit
rr2::R2(pglmm_dur) #0.717

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


dur_correlogram_plot <- ggplot() +
  geom_line(fit3_df,
            mapping = aes(x = distance, y = correlation)) +
  geom_point(fit3_df,
             mapping = aes(x = distance, y = correlation, fill = p.value),
             size = 3, shape = 21, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("white", "black")) +
  labs(x = "Distance", y = "Moran's I", fill = "P-Value") +
  theme_classic()

dur_correlogram_plot # no autocorrelation found at closest spatial lag, so spatial model not made

# save plot as Rdata to make multipanel Moran's I plot for SI
save(dur_correlogram_plot, file = "figures/dur_correlogram_plot.Rdata")
