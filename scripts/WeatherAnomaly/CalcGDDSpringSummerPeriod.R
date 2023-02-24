setwd("C:\\Users\\linds\\rob_gdd\\redo")


dir()

temp<-read.csv("temp_allyear_byidcells2.csv",stringsAsFactors = F)
View(temp)

library(tidyverse)
library(dplyr)

names(temp)

## March 21 to June 30, perpetual 80 and 181, leap year 81 and 182
## 1948 to 2016 
## leap years 1948, 1952, 1956, 1960, 1964, 1968, 1972, 1976, 1980, 1984, 1988, 1992, 1996, 2000, 2004, 2008, 2012, 2016

yr2<-as.list(seq(from = 1948, to=2016, by = 4))

lp<- temp %>% filter(year %in% yr2)
dim(lp)
unique(lp$year)

out<-lp %>% filter(DOY > 80 & DOY < 183)
dim(out)
View(out)
range(out$DOY)
write.csv(out,"leap_yr_filter_DOY.csv", row.names = F)


## now filter for perpetual years, so not leap years

yr2<-as.list(seq(from = 1948, to=2016, by = 4))

perpetual<- temp %>% filter(!year %in% yr2)

unique(perpetual$year)

out2<-perpetual %>% filter(DOY > 79 & DOY < 182)
dim(out2)
View(out2)
range(out2$DOY)
write.csv(out2,"perpetual_yr_filter_DOY.csv", row.names = F)


### rbind and maybe sort by year, id_cell, DOY ###

out3<-rbind(out,out2)
dim(out3)
out4<-arrange(out3,year, id_cell, DOY)

unique(out4$year)
range(out4$DOY)

write.csv(out4,"filtered_yr_DOY_final.csv",row.names=F)

### GDD prep ###

## set tbase to 5, multiply tmin and tmax by 100

out4$tbase <- 5
out4$tmin2<-out4$tmin*100
out4$tmax2<-out4$tmax*100

## then if tmin2 < 5 (tbase), set to 5 and if tmax2 > 38 (tmax_base) set to 38

out4$tmin3<-ifelse(out4$tmin2 < 5.1, 5, out4$tmin2)
out4$tmax3<-ifelse(out4$tmax2 > 38.1, 38, out4$tmax2)

## calc tave using new tmin and tmax

out4$tave3<-(out4$tmax3 + out4$tmin3)/2
View(out4)

## GDD = tmean - tbase, if tmean > tbase
## GDD = 0, if tmean is < tbase

out4$GDD<-ifelse(out4$tave3 > out4$tbase, out4$tave3 - out4$tbase, 0)
View(out4)
write.csv(out20,"GDD_tbase5_tmax_base38.csv",row.names = F)


## calculate accumulated GDD ##

out5<-out4 %>% group_by(id_cell, year) %>%
  summarise(accum_GDD = sum(GDD))

View(out5)

write.csv(out5,"GDD_year_cell_id.csv",row.names = F)

### cross-check summarize with single cell and year(my filter is only letting me filter one colum at a time ?! I know it's me but oof)

gdd_check<-out4 %>% filter(id_cell == 60)
gdd_check2<- gdd_check %>% filter(year == 1960)
View(gdd_check2)

gdd_sum<-sum(gdd_check2$GDD)
gdd_sum
# 1770.032

gdd_sum_out5<-out5 %>% filter(year == 1960)
gdd_sum_out5B<-gdd_sum_out5 %>% filter(id_cell == 60)
gdd_sum_out5B
View(gdd_sum_out5B)
# 1770.032