library(data.table)
library(tidyverse)
library(splitstackshape)
library(pbapply)

## read in data
leps <- fread("data/munged/LepsByGrid.csv") %>% 
  mutate(binomial = paste(genus, specificEpithet, sep = " "))

## only adult data -- to the best of our knowledge
notAdults <- c("Nymph", "Immature", "Pupa", "Larva", "Egg",
               "Caterpillar", "Juvenile", "juvenile; larva", 
               "1st instar", "Larval", "larval", "pupal case",
               "Larvae", "pupa", "Pupae", "egg", "juvenile;larva",
               "Inmaduro", "juvenile", "Exuviae", "exuviae", "Egg sac",
               "Cocoon", "cocoon", "larva", "1 exuvia w/egg", "Pupal case + larva exuvia",
               "Case", "Larval?", "larval case", "larvae", "larval in sac", "Pupa case")
leps <- leps %>% 
  filter(!lifeStage %in% notAdults)

ggplot(leps, mapping = aes(x = doy)) +
  geom_histogram(bins = 365) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw()

# lines were used to visualize what the data looked like for first of years
#doy1 <- filter(leps, doy == 1)
#doy32 <- filter(leps, doy == 32) 

spp_list <- read.csv("data/traits/spp_list_traits_mike_withValidNames.csv")
spp_list <- spp_list$scientific_name

names_df <- spp_list %>% 
  dplyr::select(scientific_name, Syn, validName)

## filter lep obs to spp_list
leps <- leps %>% 
  #  filter(binomial %in% spp_list) %>% 
  mutate(gcode = paste(binomial, year, id_cells, sep = ".")) %>% 
  filter(!is.na(doy)) %>% 
  filter(day > 1) 

ggplot(leps, mapping = aes(x = doy)) +
  geom_histogram(bins = 365) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw()

# join to names by scientific_name
leps_nameJoin <- left_join(leps, names_df, 
                           by = c("binomial" = "scientific_name")) %>% 
  dplyr::distinct(id, .keep_all = T)

leps_nameJoin2 <- left_join(leps_nameJoin, names_df, 
                           by = c("binomial" = "Syn")) %>% 
  dplyr::distinct(id, .keep_all = T)

leps_validNames <- leps_nameJoin2 %>% 
  mutate(validName = if_else(condition = is.na(validName.x),
                             true = validName.y,
                             false = validName.x))

leps_validNames <- leps_validNames %>% 
  filter(!is.na(validName)) %>% 
  dplyr::select(-validName.x, -validName.y)

## find species, cell, year combination with enough data
leps_sum <- leps_validNames %>% 
  group_by(year, validName, id_cells) %>% 
  summarize(ndstcol = n_distinct(recordedBy), 
            dstdoy=n_distinct(doy), 
            nObs=n())

leps_sum <- leps_sum %>% 
  dplyr::filter(!is.na(id_cells))

## join with trait data
traits <- read.csv("data/traits/spp_list_traits_mike_withValidNames.csv") %>% 
  dplyr::distinct(validName, .keep_all = T)
all_traits <- traits %>% 
  dplyr::filter(!is.na(overwinteringStage),
                !is.na(voltinism),
                !is.na(diurnality),
                overwinteringStage != "",
                voltinism != "",
                diurnality != "")

leps_sum_traits <- left_join(leps_sum, all_traits)

enough_data <- leps_sum_traits %>% 
  filter(voltinism == "U" & ndstcol >= 3  & dstdoy >= 4 & nObs >= 5 |
         voltinism != "U" & ndstcol >= 3 & dstdoy >= 8 & nObs >= 10) %>% 
  mutate(gcode = paste(validName, year, id_cells, sep = "."))

enough_data <- enough_data %>% 
  dplyr::select(year, validName, id_cells, ndstcol, dstdoy, nObs)

leps_effort <- left_join(leps_validNames, enough_data) %>% 
  filter(!is.na(nObs),
         !is.na(id_cells),
         !is.na(binomial))

## move to bear
write.csv(leps_effort, file = "outputs/occurrence_records_forPhenometrics.csv", row.names = F)

## split into list of dataframes per species
data_list <- split(leps_effort, 
                   f = list(leps_effort$year,
                            leps_effort$id_cells,
                            leps_effort$binomial), 
                   drop = TRUE)

# phenometric function
pheno_fun <- function(x){
  
  df_dup <- data_list[[x]]
  
  q5  <- tryCatch(phenesse::quantile_ci(observations = df_dup$doy, percentile = 0.05),
                  error = function(e) data.frame(estimate = NA, low_ci = NA, high_ci = NA))
  q50 <- tryCatch(phenesse::quantile_ci(observations = df_dup$doy, percentile = 0.5),
                  error = function(e) data.frame(estimate = NA, low_ci = NA, high_ci = NA))
  q95 <- tryCatch(phenesse::quantile_ci(observations = df_dup$doy, percentile = 0.95),
                  error = function(e) data.frame(estimate = NA, low_ci = NA, high_ci = NA))
  
  rdf <- data.frame(q5 = q5$estimate, q5_low = q5$low_ci, q5_high = q5$high_ci,
                    q50 = q50$estimate, q50_low = q50$low_ci, q50_high = q50$high_ci, 
                    q95 = q95$estimate, q95_low = q95$low_ci, q95_high = q95$high_ci) %>% 
    mutate(year = df_dup$year[1],
           HEXcell = df_dup$id_cells[1],
           code = df_dup$binomial[1],
           ndstcol = df_dup$ndstcol[1],
           dstdoy = df_dup$dstdoy[1],
           nObs = df_dup$nObs[1])
  
  return(rdf)
  
}

## generate phenometrics
pheno_out_list <- pblapply(1:length(data_list), FUN = pheno_fun)
pheno_df_output <- pheno_out_list %>% 
  lapply(., function(x){if(class(x) == "data.frame") return(x) }) %>% 
  bind_rows()

write.csv(x = pheno_df_output, "outputs/lep_phenometrics.csv", row.names = F)

# deal with NAs?

errors <- pheno_df_output %>% 
  filter(is.na(q5) | is.na(q50) | is.na(q95)) %>% 
  mutate(gcode = paste(code, year, HEXcell, sep = "."))  

e5 <- errors %>% 
  filter(is.na(q5))

e50 <- errors %>% 
  filter(is.na(q50))

e95 <- errors %>% 
  filter(is.na(q95))

## filter to error data
opp_data_e5 <- leps_effort %>% 
  filter(gcode %in% e5$gcode)

opp_data_e50 <- leps_effort %>% 
  filter(gcode %in% e50$gcode)

opp_data_e95 <- leps_effort %>% 
  filter(gcode %in% e95$gcode)


# start with e5
data_list_e5 <- split(opp_data_e5, 
                      f = list(opp_data_e5$year,
                               opp_data_e5$id_cells,
                               opp_data_e5$binomial), 
                      drop = TRUE)


# estimate phenometrics for terms that errored out earlier
pheno_fun_e5 <- function(x){
  
  df_dup <- data_list_e5[[x]]
  
  q5  <- quantile(x = df_dup$doy, probs = 0.05)
  
  rdf <- data.frame(q5 = q5, q5_low = NA, q5_high = NA) %>% 
    mutate(year = df_dup$year[1],
           HEXcell = df_dup$id_cells[1],
           code = df_dup$binomial[1],
           ndstcol = df_dup$ndstcol[1],
           dstdoy = df_dup$dstdoy[1],
           nObs = df_dup$nObs[1])
  
  return(rdf)
  
}

q5_out_list <- lapply(X = 1:length(data_list_e5), FUN = pheno_fun_e5)
q5_df_output <- q5_out_list %>% 
  lapply(., function(x){if(class(x) == "data.frame") return(x) }) %>% 
  bind_rows() %>% 
  mutate(gcode = paste(code, HEXcell, year, sep = "."))

pheno_df_output_q5 <- left_join(pheno_df_output, q5_df_output,
                                by = c("q5_low", "q5_high", "year", "HEXcell", "code", 
                                       "ndstcol", "dstdoy", "nObs"))

pheno_df_output_q5 <- pheno_df_output_q5 %>%  
  mutate(q5 = if_else(condition = is.na(q5.x), 
                      true= q5.y,
                      false = q5.x))


## now let's deal with the 50% errors
data_list_e50 <- split(opp_data_e50, 
                       f = list(opp_data_e50$year,
                                opp_data_e50$id_cells,
                                opp_data_e50$binomial), 
                       drop = TRUE)


# estimate phenometrics for terms that errored out earlier for q50
pheno_fun_e50 <- function(x){
  
  df_dup <- data_list_e50[[x]]
  
  q50  <- quantile(x = df_dup$doy, probs = 0.50)
  
  rdf <- data.frame(q50 = q50, q50_low = NA, q50_high = NA) %>% 
    mutate(year = df_dup$year[1],
           HEXcell = df_dup$id_cells[1],
           code = df_dup$binomial[1],
           ndstcol = df_dup$ndstcol[1],
           dstdoy = df_dup$dstdoy[1],
           nObs = df_dup$nObs[1])
  
  return(rdf)
  
}

q50_out_list <- lapply(X = 1:length(data_list_e50), FUN = pheno_fun_e50)
q50_df_output <- q50_out_list %>% 
  lapply(., function(x){if(class(x) == "data.frame") return(x) }) %>% 
  bind_rows() %>% 
  mutate(gcode = paste(code, HEXcell, year, sep = "."))

pheno_df_output_q5q50 <- left_join(pheno_df_output_q5, q50_df_output,
                                   by = c("q50_low", "q50_high", "year", "HEXcell", "code", 
                                          "ndstcol", "dstdoy", "nObs"))

pheno_df_output_q5q50 <- pheno_df_output_q5q50 %>%  
  mutate(q50 = if_else(condition = is.na(q50.x), 
                       true= q50.y,
                       false = q50.x))

## last round for the q95

## now let's deal with the 95% errors
data_list_e95 <- split(opp_data_e95, 
                       f = list(opp_data_e95$year,
                                opp_data_e95$id_cells,
                                opp_data_e95$binomial), 
                       drop = TRUE)


# estimate phenometrics for terms that errored out earlier for q95
pheno_fun_e95 <- function(x){
  
  df_dup <- data_list_e95[[x]]
  
  q95  <- quantile(x = df_dup$doy, probs = 0.95)
  
  rdf <- data.frame(q95 = q95, q95_low = NA, q95_high = NA) %>% 
    mutate(year = df_dup$year[1],
           HEXcell = df_dup$id_cells[1],
           code = df_dup$binomial[1],
           ndstcol = df_dup$ndstcol[1],
           dstdoy = df_dup$dstdoy[1],
           nObs = df_dup$nObs[1])
  
  return(rdf)
  
}

q95_out_list <- lapply(X = 1:length(data_list_e95), FUN = pheno_fun_e95)
q95_df_output <- q95_out_list %>% 
  lapply(., function(x){if(class(x) == "data.frame") return(x) }) %>% 
  bind_rows() %>% 
  mutate(gcode = paste(code, HEXcell, year, sep = "."))

pheno_df_output_q5q50q95 <- left_join(pheno_df_output_q5q50, q95_df_output,
                                      by = c("q95_low", "q95_high", "year", "HEXcell", "code", 
                                             "ndstcol", "dstdoy", "nObs"))

pheno_df_output_q5q50q95 <- pheno_df_output_q5q50q95 %>%  
  mutate(q95 = if_else(condition = is.na(q95.x), 
                       true= q95.y,
                       false = q95.x))

## select only necesary column names
cleaned_pheno_df <- pheno_df_output_q5q50q95 %>% 
  select(q5, q5_low, q5_high, q50, q50_low, q50_high, q95, q95_low, q95_high,
         year, HEXcell, code, ndstcol, dstdoy, nObs) %>% 
  mutate(q5_CI = q5_high - q5_low, 
         q50_CI = q50_high - q50_low, 
         q95_CI = q95_high - q95_low)

## replace NA CIs with max of CIs
cleaned_pheno_df <- cleaned_pheno_df %>% 
  mutate(q5_CI = if_else(condition = is.na(q5_CI), 
                         true = max(q5_CI, na.rm = T), false = q5_CI),
         
         q50_CI = if_else(condition = is.na(q50_CI), 
                          true = max(q50_CI, na.rm = T), false = q50_CI),
         
         q95_CI = if_else(condition = is.na(q95_CI), 
                          true = max(q95_CI, na.rm = T), false = q95_CI),
         
  )
write.csv(x = cleaned_pheno_df, 
          file = "outputs/lep_phenometrics_noNAs.csv",
          row.names = F)
