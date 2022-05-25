library(dplyr)
library(lme4)
library(lmerTest)
library(ggplot2)

# read in data
mdf <- read.csv("data/LMM_Data/mdf.csv") %>% 
  dplyr::filter(overwinteringStage != "None") %>% 
  mutate(voltinism = if_else(condition = voltinism == "F",
                             true = "M",
                             false = voltinism))
mdf <- mdf %>% 
  mutate(volt2 = if_else(voltinism == "M",
                         true = "Multivoltine",
                         false = "Univoltine"))

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

mdf_removeOutliers <- mdf_res %>% 
  mutate(Outlier = if_else(condition = residual < 60.933 & residual > -60.933,
                           true = "No",
                           false = "Yes"))

mdf_removeOutliers$Outlier <- factor(mdf_removeOutliers$Outlier, levels = c("Yes", "No"))

mdf_removeOutliers <- mdf_removeOutliers %>% 
  mutate(yi = if_else(volt2 == "Multivoltine", 
                      true = 8,
                      false = 4))

# plot
ggplot(mdf_removeOutliers) + 
  geom_hline(aes(yintercept = yi), linetype = "dashed") +
  geom_jitter(aes(x=year, y=dstdoy, 
                  size = Outlier, color = Outlier), alpha = 0.4) + 
  scale_size_manual(values = c(4,1)) +
  scale_color_manual(values = c("gold", "dodgerblue3")) +
  labs(x = "Year", y = "Distinct observation days") +
  facet_wrap(~volt2) +
  theme_classic()

ggsave("figures/Outlier_Phenometrics.png", dpi = 600,
       width = 8, height = 4)

