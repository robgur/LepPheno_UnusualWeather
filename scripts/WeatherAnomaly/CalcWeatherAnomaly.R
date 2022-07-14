library(terra)
library(sf)
library(dplyr)
library(rnaturalearth)
library(rnaturalearthdata)
library(tidyr)
library(rts)
library(stars)
library(data.table)
library(readr)
library(lubridate)
library(magrittr)
library(anytime)

#set some overal path vars
write("TMP = /srv/duo/robgur", file=file.path('~/.Renviron'))
setwd("/srv/duo/robgur")

#set the def of USA and create gridding
sf::sf_use_s2(FALSE)
test2 <- rnaturalearth::ne_countries(country = c("United States of America", "Mexico", "Canada"),returnclass = "sf")
test2<-st_transform(test2, crs="+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
grids = st_make_grid(test2, cellsize = c(250000, 250000))
grids2 = grids[test2]
grids2 = mutate(st_sf(geometry = grids2), id_cells = 1:n())


#load concatenated file of grid/day temp values and slim it down
p_allyears_allcells_1_ <- read_csv("/srv/duo/robgur/p_allyears_allcells.csv", ",", col_names = TRUE)
p_allyears_allcells_1_$ind <- as.Date(p_allyears_allcells_1_$ind)
p_allyears_new_done2 <- read_csv("p_allyears_new_done2.csv", col_types = cols(ind = col_date(format = "%m/%d/%Y")))
temp_allyear_byidcells  <- rbind(p_allyears_allcells_1_,p_allyears_new_done2)

#make DOY and year columns for the temperature data. 
temp_allyear_byidcells2  <- temp_allyear_byidcells %>%
  rename(date=ind, value = t_ave)
temp_allyear_byidcells2$date = anydate(temp_allyear_byidcells2$date)
temp_allyear_byidcells2$DOY = yday(temp_allyear_byidcells2$date)
temp_allyear_byidcells2$year = year(temp_allyear_byidcells2$date)

#upload_pheno_data
flowering_4816 <- read_csv("/home/robgur/flowering_4816_DT15N3_cleaned.csv")
flowering_4816$onset <- as.integer(flowering_4816$onset)
flowering_4816$offset <- as.integer(flowering_4816$offset)
flowering_4816$origin <- as.Date(paste0(flowering_4816$Year, "-01-01"),tz = "UTC") - days(1)
flowering_4816$onsetdate <- as.Date(flowering_4816$onset, origin = flowering_4816$origin, tz = "UTC")
flowering_4816$offsetdate <- as.Date(flowering_4816$offset, origin = flowering_4816$origin, tz = "UTC") 

#make sure everything is a data table
flowering_4816_DT <- setDT(flowering_4816)
temp_allyear_byidcells2_DT <- setDT(temp_allyear_byidcells2)

#function to calculate sum of the unusually warms days X length before onset
yrs <- unique(temp_allyear_byidcells2_DT$year)
GetUnusualWarmPreOnset <- function(dat, lngth)  {
  vec <- vector()
  for (i in 1:nrow(dat)) {
    MonthDay = format(as.Date(dat$onsetdate[i]),"%m-%d")
    targets <- as.Date(paste(yrs,MonthDay, sep = "-"))-lngth
    dates <- (c(as.Date(sapply(targets, function(x) seq.Date(from = as.Date(x),to = as.Date((x)+lngth), by = "day")))))
    winter_temp_val2 <- temp_allyear_byidcells2_DT %>% 
      filter(id_cell==dat$id_cells[i]) %>%
      filter(date %in% dates) %>%
      mutate(MonthDay=format(as.Date(date),"%m-%d")) %>%
      group_by(MonthDay) %>%
      summarize(ave_value=mean(value)+2*sd(value))
    #all_temp_val <- temp_all_byidcells31 %>% group_by(id_cells)
    winter_temp_val_year <- temp_allyear_byidcells2_DT %>% 
      filter(id_cell==dat$id_cells[i]) %>%
      filter(between(date, as.Date(dat$onsetdate[i])-lngth, as.Date(dat$onsetdate[i]))) %>%
      mutate(MonthDay=format(as.Date(date),"%m-%d"))
    winter_temp_val_merge <- winter_temp_val_year %>% left_join(winter_temp_val2, by="MonthDay")
    winter_temp_sum <- winter_temp_val_merge %>% 
      summarise(across(value, ~ sum(. > ave_value)))  %>% pull()
    #winter_temp_val_count <-  winter_temp_val_merge %>% summarize(n_hot = sum(value > ave_value )) %>% pull
    vec <- c(vec, winter_temp_sum)
  }
  return(vec)
}

#function to calculate sum of the unusually colds days X length before onset
yrs <- unique(temp_allyear_byidcells2_DT$year)
GetUnusualColdPreOnset <- function(dat, lngth)  {
  vec <- vector()
  for (i in 1:nrow(dat)) {
    MonthDay = format(as.Date(dat$onsetdate[i]),"%m-%d")
    targets <- as.Date(paste(yrs,MonthDay, sep = "-"))-lngth
    dates <- (c(as.Date(sapply(targets, function(x) seq.Date(from = as.Date(x),to = as.Date((x)+lngth), by = "day")))))
    winter_temp_val2 <- temp_allyear_byidcells2_DT %>% 
      filter(id_cell==dat$id_cells[i]) %>%
      filter(date %in% dates) %>%
      mutate(MonthDay=format(as.Date(date),"%m-%d")) %>%
      group_by(MonthDay) %>%
      summarize(ave_value=mean(value)-(2*sd(value)))
    #all_temp_val <- temp_all_byidcells31 %>% group_by(id_cells)
    winter_temp_val_year <- temp_allyear_byidcells2_DT %>% 
      filter(id_cell==dat$id_cells[i]) %>%
      filter(between(date, as.Date(dat$onsetdate[i])-lngth, as.Date(dat$onsetdate[i]))) %>%
      mutate(MonthDay=format(as.Date(date),"%m-%d"))
    winter_temp_val_merge <- winter_temp_val_year %>% left_join(winter_temp_val2, by="MonthDay")
    winter_temp_sum <- winter_temp_val_merge %>% 
      summarise(across(value, ~ sum(. < ave_value)))  %>% pull()
    #winter_temp_val_count <-  winter_temp_val_merge %>% summarize(n_hot = sum(value > ave_value )) %>% pull
    vec <- c(vec, winter_temp_sum)
  }
  return(vec)
}



#function to calculate sum of the unusually warm days between onset and offset
yrs <- unique(temp_allyear_byidcells2_DT$year)
GetUnusualWarmOnOff <- function(dat)  {
  vec <- vector()
  for (i in 1:nrow(dat)) {
    MonthDay = format(as.Date(dat$onsetdate[i]),"%m-%d")
    targets <- as.Date(paste(yrs,MonthDay, sep = "-"))
    duration <- as.numeric(as.Date(dat$offsetdate[i])-as.Date(dat$onsetdate[i]))
    dates <- (c(as.Date(sapply(targets, function(x) seq.Date(from = as.Date(x),to = as.Date((x)+duration), by = "day")))))
    winter_temp_val2 <- temp_allyear_byidcells2_DT %>% 
      filter(id_cell==dat$id_cells[i]) %>%
      filter(date %in% dates) %>%
      group_by(DOY) %>%
      summarize(ave_value=mean(value)+ (2*sd(value)))
    #all_temp_val <- temp_all_byidcells31 %>% group_by(id_cells)
    winter_temp_val_year <- temp_allyear_byidcells2_DT %>% 
      filter(id_cell==dat$id_cells[i]) %>%
      filter(between(date, as.Date(dat$onsetdate[i]), as.Date(dat$offsetdate[i])))
    winter_temp_val_merge <- winter_temp_val_year %>% left_join(winter_temp_val2, by="DOY")
    winter_temp_sum <- winter_temp_val_merge %>% 
      summarise(across(value, ~ sum(. > ave_value)))  %>% pull()
    #winter_temp_val_count <-  winter_temp_val_merge %>% summarize(n_hot = sum(value > ave_value )) %>% pull
    vec <- c(vec, winter_temp_sum)
  }
  return(vec)
}

#function to calculate sum of the unusually cold days between onset and offset
yrs <- unique(temp_allyear_byidcells2_DT$year)
GetUnusualColdOnOff <- function(dat)  {
  vec <- vector()
  for (i in 1:nrow(dat)) {
    MonthDay = format(as.Date(dat$onsetdate[i]),"%m-%d")
    targets <- as.Date(paste(yrs,MonthDay, sep = "-"))
    duration <- as.numeric(as.Date(dat$offsetdate[i])-as.Date(dat$onsetdate[i]))
    dates <- (c(as.Date(sapply(targets, function(x) seq.Date(from = as.Date(x),to = as.Date((x)+duration), by = "day")))))
    winter_temp_val2 <- temp_allyear_byidcells2_DT %>% 
      filter(id_cell==dat$id_cells[i]) %>%
      filter(date %in% dates) %>%
      group_by(DOY) %>%
      summarize(ave_value=mean(value) - (2*sd(value)))
    #all_temp_val <- temp_all_byidcells31 %>% group_by(id_cells)
    winter_temp_val_year <- temp_allyear_byidcells2_DT %>% 
      filter(id_cell==dat$id_cells[i]) %>%
      filter(between(date, as.Date(dat$onsetdate[i]), as.Date(dat$offsetdate[i])))
    winter_temp_val_merge <- winter_temp_val_year %>% left_join(winter_temp_val2, by="DOY")
    winter_temp_sum <- winter_temp_val_merge %>% 
      summarise(across(value, ~ sum(. < ave_value)))  %>% pull()
    #winter_temp_val_count <-  winter_temp_val_merge %>% summarize(n_hot = sum(value > ave_value )) %>% pull
    vec <- c(vec, winter_temp_sum)
  }
  return(vec)
}

#function to cl
GetGDDonset <- function(dat,GDD_min, lngth)  {
  vec <- vector()
  for (i in 1:nrow(dat)) {
    onset_temp_val <- temp_allyear_byidcells2_DT %>% 
      filter(id_cell== dat$id_cells[i]) %>% 
      filter(between(date, as.Date(dat$onsetdate[i])-lngth, as.Date(dat$onsetdate[i])))
    onset_temp_val <- as.data.table(onset_temp_val) 
    onset_temp_val[, GDD := value*100 - GDD_min]
    onset_temp_val[onset_temp_val$GDD < 0, "GDD"] <- 0 
    onset_temp_sum <- onset_temp_val[,sum(GDD)]
    vec <- c(vec, onset_temp_sum)
    }
  return(vec)
 }
 

GetGDDoffset <- function(dat,GDD_min)  {
  vec <- vector()
  for (i in 1:nrow(dat)) {
    offset_temp_val <- temp_allyear_byidcells2_DT %>% 
      filter(id_cell==dat$id_cells[i]) %>% 	
      filter(date >= dat$onsetdate[i], date <=dat$offsetdate[i])
    offset_temp_val <- as.data.table(offset_temp_val) 
    offset_temp_val[, GDD := value*100 - GDD_min]
    offset_temp_val[offset_temp_val$GDD < 0, "GDD"] <- 0     
    offset_temp_sum <-   offset_temp_val[,sum(GDD)]
    vec <- c(vec, offset_temp_sum)
    }
 return(vec)
}

#GetWinterCD <- function(dat)  {
#  vec <- vector()
#  for (i in 1:nrow(dat)) {
#    winter_temp_val <- temp_allyear_byidcells2_DT %>%
#      filter(id_cell==dat$id_cells[i]) %>%
#      filter(between(date, as.Date(paste0((dat$Year[i]-1),"-11-01")), as.Date(paste0(dat$Year[i],"-02-15"))))
#    winter_temp_val <- as.data.table(winter_temp_val)
#    winter_temp_sum <- winter_temp_val %>% filter(value >0 & value <0.072) %>% count()
#    vec <- c(vec, winter_temp_sum)
#  }
#  return(vec)
#}
#WinterDay <- as.data.frame(unlist(WinterCD))
#WinterCD <- WinterDay[,1]



#GetCVonset <- function(dat,lngth)  {
#  vec <- vector()
#  for (i in 1:nrow(dat)) {
#    onset_temp_val <- temp_all_byidcells31 %>% 
#      filter(id_cells== dat$dat_to_use.id_cells[i]) %>% 
#      filter(between(date, as.Date(dat$onsetdate[i])-lngth, as.Date(dat$onsetdate[i])))
#    onset_temp_val <- as.data.table(onset_temp_val)
#    onset_CV <- onset_temp_val[,sd(valueinK)]  / onset_temp_val[,mean(valueinK)] 
#    vec <- c(vec, onset_CV)
#  }
#  return(vec)
#}


#GetCVoffset <- function(dat)  {
#  vec <- vector()
#  for (i in 1:nrow(dat)) {
#    offset_temp_val <- temp_all_byidcells31 %>% 
#      filter(id_cells==dat$dat_to_use.id_cells[i]) %>% 	
#      filter(date >= dat$onsetdate[i], date <=dat$offsetdate[i])
#    offset_temp_val <- as.data.table(offset_temp_val)
#    offset_CV <- offset_temp_val[,sd(valueinK)]  / offset_temp_val[,mean(valueinK)] 
#    vec <- c(vec, offset_CV)
#  }
#  return(vec)
#}
