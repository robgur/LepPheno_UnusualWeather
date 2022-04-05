library(lme4)
library(lmerTest)

mdf <- read.csv("data/LMM_Data/mdf.csv") %>% 
  dplyr::filter(overwinteringStage != "None")

mdf <- mdf %>% 
  mutate(annualTemp = scale(annualTemp))

m <- lmer(formula = q5 ~ annualTemp + 
            (1|id_cells) + (1|validName) + 
            (0 + annualTemp|validName), 
          data = mdf)
summary(m)

r <- residuals(m)

mdf_res <- mdf %>% 
  mutate(residual = r)

hist(mdf_res$residual)
3*sd(mdf_res$residual)
