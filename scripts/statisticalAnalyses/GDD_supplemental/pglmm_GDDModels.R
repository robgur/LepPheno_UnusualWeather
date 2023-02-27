library(phyr)
library(car)
library(ggplot2)
library(ncf)
library(fpp2)

# get data ready for pglmm models
source('scripts/pglmm_dataSetUp_GDDModels.R')

#___________________________duraton
#model_no_anomaly_reduced_vif_ok
dur_m_noAnom_reduced <- pglmm(duration ~ accum_GDD_sc + annualPrec + 
                               voltinism + dstdoy + overwinteringStage
                             + (1 | unique_name__) + accum_GDD_sc:voltinism,
                             data = mdf_phylo,
                             cov_ranef = list(unique_name = tt),
                             bayes = T)


#model_with_anomaly_reduced_vif_ok
dur_m_reduce <- pglmm(duration ~ accum_GDD_sc + annualPrec + 
                       UnusualColdSum + UnusualWarmSum + voltinism +
                       dstdoy + overwinteringStage + 
                       (1 | unique_name__) + (1 | id_cells) + 
                       accum_GDD_sc:voltinism 
                     + annualPrec:UnusualWarmSum + UnusualColdSum:voltinism, 
                     data = mdf_phylo,
                     cov_ranef = list(unique_name = tt),
                     bayes = T)

plot_bayes(dur_m_reduce)
ggsave(filename = "Figures/duration_posteriorDistributions_GDDModels.pdf",
       width = 8, height = 9)
# WAIC(dur_m_reduce)
# 4949
# WAIC(dur_m_noAnom_reduced)
#5032 


#_________________________________onset

#model_no_anomaly_onset_reduced_vif_ok
onset_m_noAnom_reduced <- pglmm(onset ~ accum_GDD_sc + overwinteringStage
                               + (1 | unique_name__) + (1 | id_cells),
                               data = mdf_phylo,
                               cov_ranef = list(unique_name = tt),
                               bayes = T)


#model_onset_with_anomaly_reduced_vif_ok
onset_m_reduce <- pglmm(onset ~ accum_GDD_sc + UnusualColdSum + UnusualWarmSum + 
                         voltinism + overwinteringStage
                       + (1 | unique_name__) + (1 | id_cells) + 
                         UnusualColdSum:UnusualWarmSum,
                       data = mdf_phylo,
                       cov_ranef = list(unique_name = tt),
                       bayes = T)


# WAIC(onset_m_reduce)
# 4621
#WAIC(onset_m_noAnom_reduced)
#4673

plot_bayes(onset_m_reduce)
ggsave(filename = "Figures/onset_posteriorDistributions_GDDModels.pdf",
       width = 8, height = 9)

#_________________________________offset

#model_offset_best_reduced_noAnom_vif_ok
offset_m_noAnom_reduced <- pglmm(offset ~ accum_GDD_sc + annualPrec + 
                                   voltinism + dstdoy + overwinteringStage + 
                                  (1 | unique_name__) + accum_GDD_sc:voltinism + 
                                   annualPrec:voltinism,
                                 data = mdf_phylo,
                                 cov_ranef = list(unique_name = tt),
                                 bayes = T)


#model_offset_best_reduced_Anomaly_vif_ok
offset_m_reduce <- pglmm(offset ~ accum_GDD_sc + annualPrec + UnusualColdSum + 
                           UnusualWarmSum + voltinism + dstdoy +
                           overwinteringStage + 
                           (1 | unique_name__) + (1 | id_cells)
                        + accum_GDD_sc:UnusualColdSum + annualPrec:UnusualWarmSum
                        + UnusualColdSum:voltinism + UnusualWarmSum:voltinism,
                        data = mdf_phylo,
                        cov_ranef = list(unique_name = tt),
                        bayes = T)

summary(offset_m_noAnom_reduced)
#WAIC 4841

summary(offset_m_reduce)
# WAIC 4813

plot_bayes(offset_m_reduce)
ggsave(filename = "Figures/offset_posteriorDistributions_GDDModels.pdf",
       width = 8, height = 9)
