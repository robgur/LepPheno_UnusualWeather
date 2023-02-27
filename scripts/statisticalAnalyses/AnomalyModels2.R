library(lme4)
library(lmerTest)
library(car)
library(performance)
library(readr)

flowering_4816_DT_15N3 <- read_csv("flowering_4816_DT_15N3.csv")

#prelim_model_duration_no_anomaly
dur_m_noAnom <- lmer(duration ~ (annualTemp + annualPrec + voltinism + Seas)^2 
                     + tempSeas + precSeas + dstdoy + ndstcol + overwinteringStage + + diurnality + 
                       (1 | validName) + (1 | id_cells),
                     data = flowering_4816_DT_15N3, 
                     REML = FALSE, lmerControl(optimizer = "bobyqa"))

dur_m_noAnom_reduce <- step(dur_m_noAnom)

#run_reduced_model(not_shown)_and_check_vifs
#car::vif(dur_m_noAnom_reduce)
#Seas:voltinism has very high VIFs --- this interaction is not crucial to include so re-run without them, and VIFs all below 5

#final_duration_model_no_anomaly_vif_under_5
dur_m_noAnom <- lmer(duration ~ annualTemp + annualPrec + voltinism + Seas + dstdoy + overwinteringStage + (1 | validName), 
                     data = flowering_4816_DT_15N3, REML = FALSE, lmerControl(optimizer = "bobyqa"))

#____

#prelim_model_duration
dur_m <- lmer(duration ~ (annualTemp + annualPrec + UnusualColdSum + UnusualWarmSum + voltinism + Seas)^2
              + tempSeas + precSeas + dstdoy + ndstcol + overwinteringStage + + diurnality + 
                (1 | validName) + (1 | id_cells), 
              data = flowering_4816_DT_15N3, REML = FALSE, lmerControl(optimizer = "bobyqa"))

dur_m_reduce <- step(dur_m)

#run_reduced_model(not_shown)_and_check_vifs
#Seas:voltinism has very high VIFs ---  this interaction is not crucial to include so re-run without them, and VIFs all below 5

#final_duration_model_vif_under_5
dur_m <- lmer(duration ~ annualTemp + annualPrec + UnusualColdSum + UnusualWarmSum +
                voltinism + Seas + dstdoy + ndstcol + overwinteringStage + annualTemp:UnusualColdSum + annualTemp:UnusualWarmSum
              + annualPrec:UnusualWarmSum + UnusualColdSum:UnusualWarmSum + UnusualColdSum:voltinism +  (1 | validName),
              data = flowering_4816_DT_15N3, REML = FALSE, lmerControl(optimizer = "bobyqa"))

AIC(dur_m_noAnom)
#7804.393
AIC(dur_m)
#7634.365

#____onset nodels

onset_m_noAnom <- lmer(onset ~ (annualTemp + annualPrec + voltinism + Seas)^2 
+ tempSeas + precSeas + dstdoy + ndstcol + overwinteringStage + diurnality + 
  (1 | validName) + (1 | id_cells),
data = flowering_4816_DT_15N3, 
REML = FALSE, lmerControl(optimizer = "bobyqa"))

onset_m_noAmom_reduce <- step(onset_m_noAnom)


#onset_model_noAnom_final_model
onset_m_noAnom_final <- lmer(onset ~ annualTemp + annualPrec + voltinism + Seas + tempSeas + dstdoy + overwinteringStage + 
                          (1 | validName) + (1 | id_cells) + annualTemp:annualPrec + annualTemp:voltinism,
                          data = flowering_4816_DT_15N3, 
                          REML = FALSE, lmerControl(optimizer = "bobyqa"))

#check vifs_are_good_final_model




#prelim_onset_model_with anomaly
onset_m <- lmer(onset ~ (annualTemp + annualPrec + UnusualColdDaysOnset  + UnusualWarmDaysOnset + voltinism + Seas)^2
                + tempSeas + precSeas + dstdoy + ndstcol + overwinteringStage + diurnality +
                  (1 | validName) + (1 | id_cells), 
                data = flowering_4816_DT_15N3, REML = FALSE, lmerControl(optimizer = "bobyqa"))

onset_m_reduce <- step(onset_m)

#run_reduced_model(not_shown)_and_check_vifs

#final_onset_model

onset_m_final <- lmer(onset ~ annualTemp + annualPrec + UnusualColdDaysOnset + UnusualWarmDaysOnset + voltinism + Seas + tempSeas + dstdoy +
                        overwinteringStage + (1 | validName) + (1 | id_cells) + annualTemp:annualPrec + annualTemp:UnusualColdDaysOnset
                      + annualTemp:UnusualWarmDaysOnset + annualTemp:voltinism + annualPrec:UnusualColdDaysOnset, 
                      data = flowering_4816_DT_15N3, REML = FALSE, lmerControl(optimizer = "bobyqa"))

#chec VIFs and all are below 5 so model is final

AIC(onset_m_noAnom_final)
#7564.268
AIC(onset_m_final)
#7553.726

#____offset models

offset_m_noAnom <- lmer(offset ~ (annualTemp + annualPrec + voltinism + Seas)^2 
                       + tempSeas + precSeas + dstdoy + ndstcol + overwinteringStage + diurnality + 
                         (1 | validName) + (1 | id_cells),
                       data = flowering_4816_DT_15N3, 
                       REML = FALSE, lmerControl(optimizer = "bobyqa"))

offset_m_noAmom_reduce <- step(offset_m_noAnom)

offset_m_noAnom_novif <- lmer(offset ~ annualTemp + annualPrec + voltinism + Seas + precSeas + dstdoy + overwinteringStage
          + (1 | validName) + (1 | id_cells) + annualTemp:annualPrec + annualTemp:voltinism + voltinism:Seas,
          data = flowering_4816_DT_15N3, 
          REML = FALSE, lmerControl(optimizer = "bobyqa"))

#check vif and voltinism:Seas is high, remove and re-run
offset_m_noAnom_final <- lmer(offset ~ annualTemp + annualPrec + voltinism + Seas + precSeas + dstdoy + overwinteringStage
                              + (1 | validName) + (1 | id_cells) + annualTemp:annualPrec + annualTemp:voltinism ,
                              data = flowering_4816_DT_15N3, 
                              REML = FALSE, lmerControl(optimizer = "bobyqa"))



#prelim_offset_model_anomaly
offset_m <- lmer(offset ~ (annualTemp + annualPrec + UnusualColdOnOff  + UnusualWarmOnOff + voltinism + Seas)^2
                 + tempSeas + precSeas + dstdoy + ndstcol + overwinteringStage + diurnality + 
                   (1 | validName) + (1 | id_cells), 
                 data = flowering_4816_DT_15N3, REML = FALSE, lmerControl(optimizer = "bobyqa"))

offset_m_reduce <- step(offset_m)

#run_reduced_model(not_shown)_and_check_vifs
#car::vif(offset_m_reduce_novif as above has Seas:voltinism has very high VIFs

#final_offset_model
offset_m_final <- lmer(offset ~ annualTemp + annualPrec + UnusualColdOnOff + UnusualWarmOnOff + voltinism + Seas + precSeas + dstdoy
                 + overwinteringStage + annualTemp:annualPrec + annualTemp:UnusualWarmOnOff + annualTemp:voltinism + annualPrec:UnusualWarmOnOff  
                 + UnusualColdOnOff:voltinism  + (1 | validName) + (1 | id_cells), 
                 data = flowering_4816_DT_15N3, REML = FALSE, lmerControl(optimizer = "bobyqa"))

AIC(offset_m_final)
#7423.725
AIC(offset_m_noAnom_final)
#7460.635
