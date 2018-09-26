#==========
#========== Preliminaries
#==========

# load packages
library(tidyverse)
library(lubridate)

path = "data/weather_raw/raw_data/"

# import weather data: 2010-2014
weather_10_14 = read_csv(paste0(path,"weather_2010_2014.csv"))
  
# import weather data: 2015-2017
weather_15_17 = read_tsv(paste0(path,"weather_2015_2017.txt"),locale = locale(decimal_mark=','))

# import weather data: winter 2017
weather_17 = read_tsv(paste0(path,"weather_winter_2017.txt"),locale = locale(decimal_mark=','))

# import weather data: 2018
weather_18 = read_tsv(paste0(path,"weather_2018.txt"),locale = locale(decimal_mark='.'))

# import radiation data: 2015-2017
rad_15_17 = read_tsv(paste0(path,"radiation_2015_2017.txt"),
                     locale = locale(decimal_mark=','))





#==========
#========== Clean and concatenate data
#==========

# combine 2015-2018 data
weather_15_18 = list(weather_15_17,weather_17,weather_18) %>%
  lapply(function(x){
    x %>%
      rename(temp = "T", wdir = D, wspeed = "F", 
             wspeedbw = "FX", gust = "FG") %>%
      select(TIMI, temp,  wdir, wspeed, wspeedbw,  gust)
  }) %>%
  bind_rows() %>%
  # add radiation data
  left_join(rad_15_17) %>%
  select(-STOD)



# combine 2010-2018 data
weather_10_18 = weather_10_14 %>% 
  # define date time object
  mutate(TIMI = ymd_hms(paste(mdy(date), hms(paste0(0,weather_10_14$time,":00:00"))))) %>%
  select(TIMI, temp, wdir, wspeed, wspeedbw, gust, RADGL) %>%
  # combine data frames
  bind_rows(weather_15_18) %>%
  # shift TIMI by one hour so that the first hour of the day is '0'
  mutate(shift = hms(("01:00:00")),
         Date_Time = TIMI - shift)

#examine
weather_10_18

# check date_time shfit
weather_10_18 %>%
  ggplot(aes(TIMI, Date_Time))+
  geom_point()+
  theme_bw()

weather_10_18 %>%
  select(Date_Time, TIMI, temp) %>%
  gather(var, value, Date_Time, TIMI) %>%
  mutate(year = year(value),
         yday = yday(value),
         hour = hour(value)) %>%
  filter(year == 2012, yday == 160) %>%
  ggplot(aes(hour, temp, linetype = var))+
  geom_line()+
  theme_bw()

# export
weather_10_18 %>%
  select(Date_Time, temp, wdir, wspeed, wspeedbw,  gust, RADGL) %>%
  write_csv("data/weather_raw/weather_2010_2018.csv")

