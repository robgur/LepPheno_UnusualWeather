library(lme4)
library(lmerTest)
library(car)
library(performance)
library(readr)
library(dplyr)
library(lubridate)
library(anytime)

setwd("C:/Users/robgu/downloads")
flowering_4816_DT_15N3F <- read_csv("flowering_4816_DT_15N3_Filtered.csv")
GDD_year_cell_id <- read_csv("GDD_year_cell_id.csv")
flowering_4816_DT_15N3FJ <- left_join(flowering_4816_DT_15N3F, GDD_year_cell_id, by = c("year2","id_cells"))  
flowering_4816_DT_15N3FJ$accum_GDD_sc <- scale(flowering_4816_DT_15N3FJ$accum_GDD)

#_________________________________________ duration

#prelim_model_duration_no_anomaly
dur_m_noAnom <- lmer(duration ~ (accum_GDD_sc + annualPrec + voltinism)^2 
                     + tempSeas + precSeas + dstdoy + ndstcol + overwinteringStage + + diurnality + 
                     (1 | validName) + (1 | id_cells),
                     data = flowering_4816_DT_15N3FJ, 
                     REML = FALSE, lmerControl(optimizer = "bobyqa"))

#reduce_model_backward_selection
step(dur_m_noAnom)

#model_reduced
dur_m_noAnom_reduced <- lmer(duration ~ accum_GDD_sc + annualPrec + voltinism + dstdoy + overwinteringStage
                     + (1 | validName) + accum_GDD_sc:voltinism,
                    data = flowering_4816_DT_15N3FJ, 
                    REML = FALSE, lmerControl(optimizer = "bobyqa"))
#no_issues_with_vif_so_reduced_model_is_final


#preliminary model duration with cold and warm days
dur_m <- lmer(duration ~ (accum_GDD_sc + annualPrec + UnusualColdSum + UnusualWarmSum + voltinism )^2
              + tempSeas + precSeas + dstdoy + ndstcol + overwinteringStage + + diurnality + 
                (1 | validName) + (1 | id_cells), 
              data = flowering_4816_DT_15N3FJ, REML = FALSE, lmerControl(optimizer = "bobyqa"))


#reduce_model_backward_selection
step(dur_m)
dur_m_reduce <- lmer(duration ~ accum_GDD_sc + annualPrec + UnusualColdSum + UnusualWarmSum + voltinism +
                       dstdoy + overwinteringStage + (1 | validName) + (1 | id_cells) + accum_GDD_sc:voltinism 
                     + annualPrec:UnusualWarmSum + UnusualColdSum:voltinism, 
                     data = flowering_4816_DT_15N3FJ, REML = FALSE, lmerControl(optimizer = "bobyqa"))
#no_issues_with_vif_so_reduced_model_is_final

AIC(dur_m_reduce)
# 5050.897
AIC(dur_m_noAnom_reduced)
#5141.229

#_____________________________________________________onset

#prelim_model_duration_no_anomaly
onset_m_noAnom <- lmer(onset ~ (accum_GDD_sc + annualPrec + voltinism)^2 
                     + tempSeas + precSeas + dstdoy + ndstcol + overwinteringStage + + diurnality + 
                       (1 | validName) + (1 | id_cells),
                     data = flowering_4816_DT_15N3FJ, 
                     REML = FALSE, lmerControl(optimizer = "bobyqa"))

#reduce_model_backward_selection
step(onset_m_noAnom)

#model_best_reduced_noAnom
onset_m_noAnom_reduced <- lmer(onset ~ accum_GDD_sc + overwinteringStage
                               + (1 | validName) + (1 | id_cells),
                             data = flowering_4816_DT_15N3FJ, 
                             REML = FALSE, lmerControl(optimizer = "bobyqa"))
#no_issues_with_vif_so_reduced_model_is_final


#preliminary model duration with cold and warm days
onset_m <- lmer(onset ~ (accum_GDD_sc + annualPrec + UnusualColdSum + UnusualWarmSum + voltinism )^2
              + tempSeas + precSeas + dstdoy + ndstcol + overwinteringStage + + diurnality + 
                (1 | validName) + (1 | id_cells), 
              data = flowering_4816_DT_15N3FJ, REML = FALSE, lmerControl(optimizer = "bobyqa"))


#reduce_model_backward_selection
step(onset_m)

#model_best_reduced_Anomaly
onset_m_reduce <- lmer(onset ~ accum_GDD_sc + UnusualColdSum + UnusualWarmSum + voltinism + overwinteringStage
                       + (1 | validName) + (1 | id_cells) + UnusualColdSum:UnusualWarmSum,
                     data = flowering_4816_DT_15N3FJ, REML = FALSE, lmerControl(optimizer = "bobyqa"))
#no_issues_with_vif_so_reduced_model_is_final


#_________________________________________________ offset

#prelim_model_duration_no_anomaly
offset_m_noAnom <- lmer(offset ~ (accum_GDD_sc + annualPrec + voltinism)^2 
                       + tempSeas + precSeas + dstdoy + ndstcol + overwinteringStage + + diurnality + 
                         (1 | validName) + (1 | id_cells),
                       data = flowering_4816_DT_15N3FJ, 
                       REML = FALSE, lmerControl(optimizer = "bobyqa"))

#reduce_model_backward_selection
step(offset_m_noAnom)

#model_best_reduced_noAnom
offset_m_noAnom_reduced <- lmer(offset ~ accum_GDD_sc + annualPrec + voltinism + dstdoy + overwinteringStage + 
                                  (1 | validName) + accum_GDD_sc:voltinism + annualPrec:voltinism,
                               data = flowering_4816_DT_15N3FJ, 
                               REML = FALSE, lmerControl(optimizer = "bobyqa"))
#no_issues_with_vif_so_reduced_model_is_final

#preliminary model duration with cold and warm days
offset_m <- lmer(offset ~ (accum_GDD_sc + annualPrec + UnusualColdSum + UnusualWarmSum + voltinism )^2
                + tempSeas + precSeas + dstdoy + ndstcol + overwinteringStage + + diurnality + 
                  (1 | validName) + (1 | id_cells), 
                data = flowering_4816_DT_15N3FJ, REML = FALSE, lmerControl(optimizer = "bobyqa"))


#reduce_model_backward_selection
step(offset_m)

#model_best_reduced_Anomaly
offset_m_reduce <- lmer(offset ~ accum_GDD_sc + annualPrec + UnusualColdSum + UnusualWarmSum + voltinism + dstdoy + overwinteringStage + (1 | validName) + (1 | id_cells)
                        + accum_GDD_sc:UnusualColdSum + annualPrec:UnusualWarmSum
                        + UnusualColdSum:voltinism + UnusualWarmSum:voltinism,
                       data = flowering_4816_DT_15N3FJ, REML = FALSE, lmerControl(optimizer = "bobyqa"))
#no_issues_with_vif_so_reduced_model_is_final

AIC(offset_m_reduce)
#4966.008
AIC(offset_m_noAnom_reduced)
#5000.017


