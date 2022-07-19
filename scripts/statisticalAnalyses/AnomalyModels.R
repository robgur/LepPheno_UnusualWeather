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
car::vif(dur_m_noAnom_reduce)
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
car::vif(dur_m_reduce_model)
#Seas:voltinism has very high VIFs ---  this interaction is not crucial to include so re-run without them, and VIFs all below 5

#final_duration_model_vif_under_5
dur_m <- lmer(duration ~ annualTemp + annualPrec + UnusualColdSum + UnusualWarmSum +
voltinism + Seas + dstdoy + ndstcol + overwinteringStage + annualTemp:UnusualColdSum + annualTemp:UnusualWarmSum
+ annualPrec:UnusualWarmSum + UnusualColdSum:UnusualWarmSum + UnusualColdSum:voltinism +  (1 | validName),
data = flowering_4816_DT_15N3, REML = FALSE, lmerControl(optimizer = "bobyqa"))

#____

#prelim_onset_model
onset_m <- lmer(onset ~ (annualTemp + annualPrec + UnusualColdDaysOnset  + UnusualWarmDaysOnset + voltinism + Seas)^2
+ tempSeas + precSeas + dstdoy + ndstcol + overwinteringStage + diurnality +
(1 | validName) + (1 | id_cells), 
data = flowering_4816_DT_15N3, REML = FALSE, lmerControl(optimizer = "bobyqa"))

onset_m_reduce <- step(onset_m)

#run_reduced_model(not_shown)_and_check_vifs
car::vif(onset_m_reduce_model)
#All_VIFs_under_5!

#final_onset_model
offset_m <- lmer(offset ~ annualTemp + annualPrec + UnusualColdOnOff + UnusualWarmOnOff + voltinism + Seas + precSeas + dstdoy +
+ overwinteringStage  + annualTemp:annualPrec + annualTemp:UnusualWarmOnOff + annualTemp:voltinism + 
annualPrec:UnusualWarmOnOff  + UnusualColdOnOff:voltinism + (1 | validName) + (1 | id_cells),
data = flowering_4816_DT_15N3, REML = FALSE, lmerControl(optimizer = "bobyqa"))

#____

#prelim_offset_model
offset_m <- lmer(offset ~ (annualTemp + annualPrec + UnusualColdOnOff  + UnusualWarmOnOff + voltinism + Seas)^2
+ tempSeas + precSeas + dstdoy + ndstcol + overwinteringStage + diurnality + 
(1 | validName) + (1 | id_cells), 
data = flowering_4816_DT_15N3, REML = FALSE, lmerControl(optimizer = "bobyqa"))

offset_m_reduce <- step(offset_m)

#run_reduced_model(not_shown)_and_check_vifs
car::vif(offset_m_reduce_model)
#Seas:voltinism has very high VIFs and annualPrec:Seas has VIF just above 5 --  these interactions are not crucial to include so re-run without them, and VIFs all below 5

#final_offset_model
offset_m <- lmer(offset ~ annualTemp + annualPrec + UnusualColdOnOff + UnusualWarmOnOff + voltinism + Seas + precSeas + dstdoy
+ overwinteringStage + annualTemp:annualPrec + annualTemp:UnusualWarmOnOff + annualTemp:voltinism + annualPrec:UnusualWarmOnOff  
+ UnusualColdOnOff:voltinism  + (1 | validName) + (1 | id_cells), 
data = flowering_4816_DT_15N3, REML = FALSE, lmerControl(optimizer = "bobyqa"))
