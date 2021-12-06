library(tidyverse)
library(phyr)
library(gridExtra)
library(INLA)

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
  mutate(q95 = q95 - q5,
         q95_CI = mean(c(q5_CI, q95_CI))) %>% 
  mutate(on_w = 1/(q5_CI + 1),
         off_w = 1/(q95_CI + 1),
         peak_w = 1/(q95_CI +1),
         q95_w = 1/(q95_CI + 1)) %>% 
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

## interaction between voltinism and annual temp
inla_mdf <- mdf_phylo %>% 
  dplyr::select(q50, annualTemp, precSeas, annualPrec, dstdoy, voltinism, overwinteringStage,
                diurnality, Seas, unique_name, id_cells, year)

pred_vals <- -3:3

term1 <- data.frame(
  q50 = rep(NA, length(pred_vals)),
  annualTemp = rep(NA, length(pred_vals)),
  annualPrec = rep(NA, length(pred_vals)),
  precSeas = rep(NA, length(pred_vals)),
  dstdoy = rep(NA, length(pred_vals)),
  voltinism = rep(NA, length(pred_vals)),
  overwinteringStage = rep("E", length(pred_vals)),
  diurnality = rep(NA, length(pred_vals)),
  Seas = rep(NA, length(pred_vals)),
  unique_name = rep(NA, length(pred_vals)),
  id_cells = rep(NA, length(pred_vals)),
  year = pred_vals)

term2 <- data.frame(
  q50 = rep(NA, length(pred_vals)),
  annualTemp = rep(NA, length(pred_vals)),
  annualPrec = rep(NA, length(pred_vals)),
  precSeas = rep(NA, length(pred_vals)),
  dstdoy = rep(NA, length(pred_vals)),
  voltinism = rep(NA, length(pred_vals)),
  overwinteringStage = rep("L", length(pred_vals)),
  diurnality = rep(NA, length(pred_vals)),
  Seas = rep(NA, length(pred_vals)),
  unique_name = rep(NA, length(pred_vals)),
  id_cells = rep(NA, length(pred_vals)),
  year = pred_vals)

term3 <- data.frame(
  q50 = rep(NA, length(pred_vals)),
  annualTemp = rep(NA, length(pred_vals)),
  annualPrec = rep(NA, length(pred_vals)),
  precSeas = rep(NA, length(pred_vals)),
  dstdoy = rep(NA, length(pred_vals)),
  voltinism = rep(NA, length(pred_vals)),
  overwinteringStage = rep("M", length(pred_vals)),
  diurnality = rep(NA, length(pred_vals)),
  Seas = rep(NA, length(pred_vals)),
  unique_name = rep(NA, length(pred_vals)),
  id_cells = rep(NA, length(pred_vals)),
  year = pred_vals)

term4 <- data.frame(
  q50 = rep(NA, length(pred_vals)),
  annualTemp = rep(NA, length(pred_vals)),
  annualPrec = rep(NA, length(pred_vals)),
  precSeas = rep(NA, length(pred_vals)),
  dstdoy = rep(NA, length(pred_vals)),
  voltinism = rep(NA, length(pred_vals)),
  overwinteringStage = rep("P", length(pred_vals)),
  diurnality = rep(NA, length(pred_vals)),
  Seas = rep(NA, length(pred_vals)),
  unique_name = rep(NA, length(pred_vals)),
  id_cells = rep(NA, length(pred_vals)),
  year = pred_vals)

term5 <- data.frame(
  q50 = rep(NA, length(pred_vals)),
  annualTemp = rep(NA, length(pred_vals)),
  annualPrec = rep(NA, length(pred_vals)),
  precSeas = rep(NA, length(pred_vals)),
  dstdoy = rep(NA, length(pred_vals)),
  voltinism = rep(NA, length(pred_vals)),
  overwinteringStage = rep("A", length(pred_vals)),
  diurnality = rep(NA, length(pred_vals)),
  Seas = rep(NA, length(pred_vals)),
  unique_name = rep(NA, length(pred_vals)),
  id_cells = rep(NA, length(pred_vals)),
  year = pred_vals)


pred.df.1 <- rbind(inla_mdf, term1)
pred.df.2 <- rbind(inla_mdf, term2)
pred.df.3 <- rbind(inla_mdf, term3)
pred.df.4 <- rbind(inla_mdf, term4)
pred.df.5 <- rbind(inla_mdf, term5)

pglmm1 <- pglmm(formula = q50 ~ annualTemp + annualPrec + 
                  voltinism + overwinteringStage + Seas +
                  annualTemp:overwinteringStage +
                  annualTemp:Seas +
                  year:overwinteringStage +
                  annualTemp:annualPrec +
                  annualPrec:year +
                  precSeas:year +
                  (1|unique_name__) + (1|id_cells),
                data = pred.df.1, 
                cov_ranef = list(unique_name = tt), 
                bayes = TRUE)

pglmm2 <- pglmm(formula = q50 ~ annualTemp + annualPrec + 
                  voltinism + overwinteringStage + Seas +
                  annualTemp:overwinteringStage +
                  annualTemp:Seas +
                  year:overwinteringStage +
                  annualTemp:annualPrec +
                  annualPrec:year +
                  precSeas:year +
                  (1|unique_name__) + (1|id_cells),
                data = pred.df.2, 
                cov_ranef = list(unique_name = tt), 
                bayes = TRUE)

pglmm3 <- pglmm(formula = q50 ~ annualTemp + annualPrec + 
                  voltinism + overwinteringStage + Seas +
                  annualTemp:overwinteringStage +
                  annualTemp:Seas +
                  year:overwinteringStage +
                  annualTemp:annualPrec +
                  annualPrec:year +
                  precSeas:year +
                  (1|unique_name__) + (1|id_cells),
                data = pred.df.3, 
                cov_ranef = list(unique_name = tt), 
                bayes = TRUE)

pglmm4 <- pglmm(formula = q50 ~ annualTemp + annualPrec + 
                  voltinism + overwinteringStage + Seas +
                  annualTemp:overwinteringStage +
                  annualTemp:Seas +
                  year:overwinteringStage +
                  annualTemp:annualPrec +
                  annualPrec:year +
                  precSeas:year +
                  (1|unique_name__) + (1|id_cells),
                data = pred.df.4, 
                cov_ranef = list(unique_name = tt), 
                bayes = TRUE)

pglmm5 <- pglmm(formula = q50 ~ annualTemp + annualPrec + 
                  voltinism + overwinteringStage + Seas +
                  annualTemp:overwinteringStage +
                  annualTemp:Seas +
                  year:overwinteringStage +
                  annualTemp:annualPrec +
                  annualPrec:year +
                  precSeas:year +
                  (1|unique_name__) + (1|id_cells),
                data = pred.df.5, 
                cov_ranef = list(unique_name = tt), 
                bayes = TRUE)

## make df of results
rdf1 <- pglmm1$inla.model$summary.linear.predictor[(nrow(mdf_phylo)+1):nrow(pred.df.1),] %>% 
  mutate(V1 = pred_vals,
         V2 = "Egg")

rdf2 <- pglmm2$inla.model$summary.linear.predictor[(nrow(mdf_phylo)+1):nrow(pred.df.1),] %>% 
  mutate(V1 = pred_vals,
         V2 = "Larvae")

rdf3 <- pglmm3$inla.model$summary.linear.predictor[(nrow(mdf_phylo)+1):nrow(pred.df.1),] %>% 
  mutate(V1 = pred_vals,
         V2 = "Migratory")

rdf4 <- pglmm4$inla.model$summary.linear.predictor[(nrow(mdf_phylo)+1):nrow(pred.df.1),] %>% 
  mutate(V1 = pred_vals,
         V2 = "Pupae")

rdf5 <- pglmm5$inla.model$summary.linear.predictor[(nrow(mdf_phylo)+1):nrow(pred.df.1),] %>% 
  mutate(V1 = pred_vals,
         V2 = "Adult")

rdf_total <- rbind(rdf1, rdf2, rdf3, rdf4, rdf5)

moy <- ggplot(rdf_total, mapping = aes(x = V1, y = mean)) +
  geom_ribbon(mapping = aes(ymin = `0.025quant`, ymax = `0.975quant`,
                            fill = V2), 
              alpha = 0.15) +
  geom_line(mapping = aes(color = V2), size = 1.05) +
  labs(x = "Year", y = "Midpoint", fill = "Diapause stage", color = "Diapause stage") + 
  scale_fill_viridis_d() +
  scale_color_viridis_d() +
  theme_bw()

moy

save(moy, file = "singleInteractionFigures/midpoint_overwintering_year.Rdata")
