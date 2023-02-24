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


#set the def of USA and create gridding
sf::sf_use_s2(FALSE)
test2 <- rnaturalearth::ne_countries(country = c("United States of America", "Mexico", "Canada"),returnclass = "sf")
test2<-st_transform(test2, crs="+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
grids = st_make_grid(test2, cellsize = c(250000, 250000))
grids2 = grids[test2]
grids2 = mutate(st_sf(geometry = grids2), id_cells = 1:n())


#load concatenated file of grid/day temp values and slim it down
p_allyears_allcells_1_ <- read_csv("p_allyears_allcells.csv", col_names = TRUE)
p_allyears_allcells_1_$ind <- as.Date(p_allyears_allcells_1_$ind)
p_allyears_new_done2 <- read_csv("p_allyears_new_done2.csv", col_types = cols(ind = col_date(format = "%m/%d/%Y")))
temp_allyear_byidcells  <- rbind(p_allyears_allcells_1_,p_allyears_new_done2)

#make DOY and year columns for the temperature data. 
temp_allyear_byidcells2  <- temp_allyear_byidcells %>%
  rename(date=ind, value = t_ave)
temp_allyear_byidcells2$date = anydate(temp_allyear_byidcells2$date)
temp_allyear_byidcells2$DOY = yday(temp_allyear_byidcells2$date)
temp_allyear_byidcells2$year = year(temp_allyear_byidcells2$date)

#set values for GDD tmin_base=5, tmax_base=38
temp_allyear_byidcells2$tmin = temp_allyear_byidcells2$values.x*100
temp_allyear_byidcells2$tmax = temp_allyear_byidcells2$values.y*100
temp_allyear_byidcells2$tmin2 <- ifelse(temp_allyear_byidcells2$tmin < 5.1, 5, temp_allyear_byidcells2$tmin)
temp_allyear_byidcells2$tmax2 <- ifelse(temp_allyear_byidcells2$tmax > 38.1, 38, temp_allyear_byidcells2$tmax)
temp_allyear_byidcells2$tave2 <-(temp_allyear_byidcells2$tmax2 + temp_allyear_byidcells2$tmin2)/2

#upload_pheno_data
flowering_4816 <- read_csv("flowering_4816_DT_15N3.csv")
flowering_4816$onset <- as.integer(flowering_4816$onset)
flowering_4816$offset <- as.integer(flowering_4816$offset)
flowering_4816$origin <- as.Date(paste0(flowering_4816$Year, "-01-01"),tz = "UTC") - days(1)
flowering_4816$onsetdate <- as.Date(flowering_4816$onset, origin = flowering_4816$origin, tz = "UTC")
flowering_4816$offsetdate <- as.Date(flowering_4816$offset, origin = flowering_4816$origin, tz = "UTC") 

#make sure everything is a data table
flowering_4816_DT <- setDT(flowering_4816)
temp_allyear_byidcells2_DT <- setDT(temp_allyear_byidcells2)

#function to calculate onset GDD 
GetGDDonset <- function(dat, lngth)  {
  vec <- vector()
  for (i in 1:nrow(dat)) {
    onset_temp_val <- temp_allyear_byidcells2_DT %>% 
      filter(id_cell== dat$id_cells[i]) %>% 
      filter(between(date, as.Date(dat$onsetdate[i])-lngth, as.Date(dat$onsetdate[i])))
    onset_temp_val <- as.data.table(onset_temp_val) 
    onset_temp_val[, GDD := tave2 - 5]
    onset_temp_val[onset_temp_val$GDD < 5, "GDD"] <- 0 
    onset_temp_sum <- onset_temp_val[,sum(GDD)]
    vec <- c(vec, onset_temp_sum)
  }
  return(vec)
}

#function to calculate offset GDD 
GetGDDoffset <- function(dat)  {
  vec <- vector()
  for (i in 1:nrow(dat)) {
    offset_temp_val <- temp_allyear_byidcells2_DT %>% 
      filter(id_cell==dat$id_cells[i]) %>% 	
      filter(date >= dat$onsetdate[i], date <=dat$offsetdate[i])
    offset_temp_val <- as.data.table(offset_temp_val) 
    offset_temp_val[, GDD := tave2 - 5]
    offset_temp_val[offset_temp_val$GDD < 5, "GDD"] <- 0 
    offset_temp_sum <- offset_temp_val[,sum(GDD)]
    vec <- c(vec, offset_temp_sum)
  }
  return(vec)
}
