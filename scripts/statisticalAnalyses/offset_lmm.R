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


### Offset modeling
off_m <- lmer(q95 ~ annualTemp + tempSeas + annualPrec + precSeas + year +
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
               lmerControl(optimizer = "bobyqa"),
               weights = off_w)

off_s <- step(off_m)
off_s

off_top_w <- lmer(q95 ~ annualTemp + tempSeas + precSeas + year +
                  dstdoy + 
                  voltinism + overwinteringStage + diurnality + Seas +
                  annualTemp:voltinism +
                  annualTemp:overwinteringStage +
                  annualTemp:Seas +
                  year:overwinteringStage +
                  annualTemp:year + 
                  tempSeas:year +
                  (1|validName) + (1|id_cells),
                data = mdf,
                REML = FALSE, 
                lmerControl(optimizer = "bobyqa"),
                weights = off_w)

car::vif(off_top_w) # remove year

off_top_w <- lmer(q95 ~ annualTemp + tempSeas + precSeas + 
                    dstdoy + 
                    voltinism + overwinteringStage + diurnality + Seas +
                    annualTemp:voltinism +
                    annualTemp:overwinteringStage +
                    annualTemp:Seas +
                    year:overwinteringStage +
                    annualTemp:year + 
                    tempSeas:year +
                    (1|validName) + (1|id_cells),
                  data = mdf,
                  REML = FALSE, 
                  lmerControl(optimizer = "bobyqa"),
                  weights = off_w)

car::vif(off_top_w) 
# is top model stable?
off_top_w_s <- step(off_top_w)
off_top_w_s # yes

### now no weights
off_m_noW <- lmer(q95 ~ annualTemp + tempSeas + annualPrec + precSeas + year +
                dstdoy + ndstcol + 
                voltinism + overwinteringStage + diurnality + Seas +
                annualTemp:voltinism +
                annualTemp:overwinteringStage +
                annualTemp:diurnality +
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

off_s_noW <- step(off_m_noW)
off_s_noW

off_top_noW <- lmer(q95 ~ annualTemp + annualPrec + 
                          dstdoy + 
                          voltinism + overwinteringStage + diurnality + Seas +
                          annualTemp:voltinism +
                          annualTemp:overwinteringStage +
                          annualTemp:Seas + 
                          (1|validName) + (1|id_cells),
                          data = mdf,
                          REML = FALSE, 
                          lmerControl(optimizer = "bobyqa"))

car::vif(off_top_noW)
# stable?
off_top_noW_s <- step(off_top_noW)
off_top_noW_s

## climate first model
off_clim <- lmer(q95 ~ annualTemp + tempSeas + annualPrec + precSeas + 
                   dstdoy + ndstcol + 
                   annualTemp:annualPrec + 
                   (1|validName) + (1|id_cells),
                 data = mdf,
                 REML = FALSE, 
                 lmerControl(optimizer = "bobyqa"))

off_clim_s <- step(off_clim)
off_clim_s

# now add traits and stuff
off_clim_traits <- lmer(q95 ~ annualTemp + annualPrec +
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
                          annualTemp:year + 
                          annualPrec:year +
                          (1|validName) + (1|id_cells),
                        data = mdf,
                        REML = FALSE, 
                        lmerControl(optimizer = "bobyqa"))

off_clim_traits_s <- step(off_clim_traits)
off_clim_traits_s

off_clim_top <- lmer(q95 ~ annualTemp + annualPrec +
                       dstdoy + 
                       voltinism + overwinteringStage + diurnality +  Seas +
                       annualTemp:voltinism +
                       annualTemp:overwinteringStage +
                       annualTemp:Seas +
                       (1|validName) + (1|id_cells),
                     data = mdf,
                     REML = FALSE, 
                     lmerControl(optimizer = "bobyqa"))

car::vif(off_clim_top)
#stable?
off_clim_top_s <- step(off_clim_top)
off_clim_top_s

AICc(off_clim_top, off_top_noW, off_top_w) ## off top_noW and off_top_w had same model
Weights(AICc(off_clim_top, off_top_noW, off_top_w)) 

summary(off_top_noW)
plot_model(off_top_noW, type = "pred", terms = c("annualTemp", "voltinism"))
plot_model(off_top_noW, type = "pred", terms = c("annualTemp", "overwinteringStage"))
plot_model(off_top_noW, type = "pred", terms = c("annualTemp", "Seas"))

tab_model(off_top_noW, file = "tables/offset_LMM.doc")

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

pglmm_offset <- pglmm(formula = q95 ~ annualTemp + annualPrec + 
                        dstdoy + 
                        voltinism + overwinteringStage + diurnality + Seas +
                        annualTemp:voltinism +
                        annualTemp:overwinteringStage +
                        annualTemp:Seas + 
                        (1|unique_name__) + (1|id_cells),
                    data = mdf_phylo, 
                    cov_ranef = list(unique_name = tt), 
                    bayes = TRUE)

summary(pglmm_offset)
rr2::R2(pglmm_offset) #0.753
im <- pglmm_offset$inla.model
im_fix <- im$summary.fixed %>% 
  tibble::rownames_to_column("Effects") %>% 
  dplyr::select(Effects, mean, `0.025quant`, `0.975quant`)
tab_df(im_fix, title = "PGLMM Offset", file = "tables/offset_PGLMM.doc")

## Model Assumption Checks
# residuals
resids <- residuals(pglmm_offset)
hist(resids)
qqnorm(resids)
qqline(resids)
shapiro.test(resids)

# spatial autocorrelation
library(geoR)
jitter_dist = 125000 # 125 km
resid_offset <- mdf_phylo %>% 
  mutate(resid = resids)

gdf <- as.geodata(resid_offset, coords.col = c("lon", "lat"), data.col = "resid")
gdf <-  geoR::jitterDupCoords(gdf, max = jitter_dist)
v1 <- variog(gdf, trend = "1st")

plot(v1)

## Goodness of fit
rr2::R2(pglmm_offset) #0.746

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


offset_correlogram_plot <- ggplot() +
  geom_line(fit3_df,
            mapping = aes(x = distance, y = correlation)) +
  geom_point(fit3_df,
             mapping = aes(x = distance, y = correlation, fill = p.value),
             size = 3, shape = 21, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("white", "black")) +
  labs(x = "Distance", y = "Moran's I", fill = "P-Value") +
  theme_classic()

offset_correlogram_plot # no autocorrelation found at closest spatial lag, so spatial model not made

# save plot as Rdata to make multipanel Moran's I plot for SI
save(offset_correlogram_plot, file = "figures/offset_correlogram_plot.Rdata")
