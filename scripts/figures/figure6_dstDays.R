library(phyr)
library(gridExtra)
library(INLA)
library(tidyverse)
library(lme4)
library(lmerTest)
library(car)
library(sjPlot)
library(MuMIn)

## read in data
mdf <- read.csv("data/LMM_Data/mdf_removeOutliersResiduals_wSeasonalityTrait.csv") %>% 
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

# drop species from phylogeny that aren't in analysis
tree_sp <- tt$tip.label
sppNotInAnalysis <- data.frame(Species = tree_sp) %>% 
  filter(!Species %in% mdf_phylo$unique_name)

tt <- ape::drop.tip(tt, tip = sppNotInAnalysis$Species)

## interaction between voltinism and annual temp
inla_mdf <- mdf_phylo %>% 
  dplyr::select(dur, q95, q5, annualTemp, tempSeas, annualPrec, precSeas,
                dstdoy, voltinism, overwinteringStage,
                diurnality, Seas, unique_name, id_cells, year)

pred_vals <- -3:3

term1 <- data.frame(
  dur = rep(NA, length(pred_vals)),
  q95 = rep(NA, length(pred_vals)),
  q5 = rep(NA, length(pred_vals)),
  annualTemp = rep(NA, length(pred_vals)),
  tempSeas = rep(NA, length(pred_vals)),
  annualPrec = rep(NA, length(pred_vals)),
  precSeas = rep(NA, length(pred_vals)),
  dstdoy = pred_vals,
  voltinism = rep(NA, length(pred_vals)),
  overwinteringStage = rep(NA, length(pred_vals)),
  diurnality = rep(NA, length(pred_vals)),
  Seas = rep(NA, length(pred_vals)),
  unique_name = rep(NA, length(pred_vals)),
  id_cells = rep(NA, length(pred_vals)),
  year =rep(NA, length(pred_vals)))

pred.df.1 <- rbind(inla_mdf, term1)

pglmm1 <-  pglmm(formula = dur ~ annualTemp + annualPrec + precSeas +
                   dstdoy +
                   voltinism + overwinteringStage + diurnality + Seas +
                   annualTemp:voltinism +
                   annualTemp:Seas + 
                   annualTemp:annualPrec +
                   (1|unique_name__) + (1|id_cells),
                 data = pred.df.1, 
                 cov_ranef = list(unique_name = tt), 
                 bayes = TRUE)


## now termination
pglmm2 <-pglmm(formula = q95 ~ annualTemp + tempSeas + precSeas +
                 dstdoy + 
                 voltinism + overwinteringStage + diurnality + Seas +
                 annualTemp:voltinism +
                 annualTemp:overwinteringStage +
                 annualTemp:Seas +
                 year:overwinteringStage +
                 annualTemp:year + 
                 tempSeas:year +
                 precSeas:year +
                 (1|unique_name__) + (1|id_cells),
               data = pred.df.1, 
               cov_ranef = list(unique_name = tt), 
               bayes = TRUE)


## now onset
pglmm3 <-   pglmm(formula = q5 ~ annualTemp + annualPrec + precSeas +
                    dstdoy +
                    overwinteringStage + 
                    annualTemp:overwinteringStage +
                    year:overwinteringStage +
                    annualTemp:annualPrec + 
                    annualPrec:year +
                    (1|unique_name__) + (1|id_cells),
                  data = pred.df.1, 
                  cov_ranef = list(unique_name = tt), 
                  bayes = TRUE)

## make df of results
rdf1 <- pglmm1$inla.model$summary.linear.predictor[(nrow(mdf_phylo)+1):nrow(pred.df.1),] %>% 
  mutate(V1 = pred_vals)

rdf2 <- pglmm2$inla.model$summary.linear.predictor[(nrow(mdf_phylo)+1):nrow(pred.df.1),] %>% 
  mutate(V1 = pred_vals)

rdf3 <- pglmm3$inla.model$summary.linear.predictor[(nrow(mdf_phylo)+1):nrow(pred.df.1),] %>% 
  mutate(V1 = pred_vals)

duration <- ggplot(rdf1, mapping = aes(x = V1, y = mean)) +
  geom_ribbon(mapping = aes(ymin = `0.025quant`, ymax = `0.975quant`), 
              alpha = 0.15) +
  geom_line(mapping = aes(), size = 1.05) +
  labs(x = "Distinct observation days", y = "Duration") + 
  theme_classic()

duration

termination <- ggplot(rdf2, mapping = aes(x = V1, y = mean)) +
  geom_ribbon(mapping = aes(ymin = `0.025quant`, ymax = `0.975quant`), 
              alpha = 0.15) +
  geom_line(mapping = aes(), size = 1.05) +
  labs(x = "Distinct observation days", y = "Termination") + 
  theme_classic()

termination

onset <- ggplot(rdf3, mapping = aes(x = V1, y = mean)) +
  geom_ribbon(mapping = aes(ymin = `0.025quant`, ymax = `0.975quant`), 
              alpha = 0.15) +
  geom_line(mapping = aes(), size = 1.05) +
  labs(x = "Distinct observation days", y = "Onset") + 
  theme_classic()

onset

cp <- egg::ggarrange(onset, termination, duration, labels = c("A","B", "C"),
                     ncol = 1)

cp

ggsave(cp, file = "figures/dstdaysCoef_revision.png", dpi = 600,
       width = 3, height = 6)

