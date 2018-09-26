#==========
#========== Preliminaries
#==========

# load packages
library(tidyverse)
library(lubridate)

# import sonde
sonde = list.files("data/sonde_raw/clean_data", full.names = T) %>%
  lapply(read_csv) %>%
  bind_rows()

# import hobo
hobo = list.files("data/hobo_raw/clean_data", full.names = T) %>%
  lapply(read_csv) %>%
  bind_rows()

# weather
weather = read_csv("data/weather_raw/weather_2010_2018.csv")





#==========
#========== Sonde data
#==========

# omit sept and oct data
# remove first and last 3 days (likely not in the water)
sonde_trim = sonde %>%
  mutate(year = year(Date_Time),
         yday = yday(Date_Time)) %>%
  filter(!(year == 2017 & Date_Time > "2017-08-19"),
         is.na(year) == F) %>%
  group_by(year) %>%
  filter(yday > min(yday) + 3,
         yday < max(yday) - 3) %>%
  ungroup()

# plot DO time series by year
sonde_trim %>%
  ggplot(aes(Date_Time, LDO))+
  facet_wrap(~year, scales = "free_x")+
  geom_line(size = 0.3)+
  geom_point(size = 0.5)+
  theme_bw()

# plot temp time series by year
sonde_trim %>%
  mutate(year = year(Date_Time)) %>%
  filter(is.na(year) == F) %>%
  ggplot(aes(Date_Time, Temp))+
  facet_wrap(~year, scales = "free_x")+
  geom_line(size = 0.3)+
  geom_point(size = 0.5)+
  theme_bw()

# omit LDO < 7
# omit Temp > 18
# remove values that deviate from daily means by greater than 2 standard deviations
sonde_clean = sonde_trim %>%
  filter(LDO > 7, Temp < 18) %>%
  group_by(year, yday) %>%
  filter(abs(LDO - mean(LDO, na.rm=T)) < 2*sd(LDO, na.rm=T),
         abs(Temp - mean(Temp, na.rm=T)) < 2*sd(Temp, na.rm=T))

# plot
sonde_clean %>%
  ggplot(aes(Date_Time, LDO))+
  facet_wrap(~year, scales = "free_x")+
  geom_line(size = 0.3)+
  geom_point(size = 0.5)+
  theme_bw()

# aggregate by hour
sonde_hourly = sonde_clean %>%
  mutate(hour = hour(Date_Time)) %>%
  group_by(year, yday, hour) %>%
  summarize(date_time = min(Date_Time),
            temp = mean(Temp, na.rm=T),
            conductivity = mean(SpCond, na.rm=T),
            turb_ntu = mean(TurbSC, na.rm=T),
            phyc = mean(PCYV, na.rm=T),
            o2_sat = mean(`LDO%`, na.rm=T),
            do = mean(LDO, na.rm=T)) %>%
  ungroup()

# plot
sonde_hourly  %>%
  ggplot(aes(date_time, do))+
  facet_wrap(~year, scales = "free_x")+
  geom_line(size = 0.3)+
  geom_point(size = 0.5)+
  theme_bw()





#==========
#========== HOBO & Weather data
#==========

# aggregate HOBO by hour
hobo_hourly = hobo %>%
  mutate(year = year(date_time),
         yday = yday(date_time),
         hour = hour(date_time)) %>%
  group_by(year, yday, hour) %>%
  summarize(date_time = min(date_time),
            lux = mean(lux, na.rm=T)) %>%
  ungroup()

# plot
hobo_hourly %>%
  ggplot(aes(date_time, lux))+
  facet_wrap(~year, scales = "free_x")+
  geom_line(size = 0.3)+
  geom_point(size = 0.5)+
  theme_bw() 

# aggregate weather by hour
weather_hourly = weather %>%
  mutate(year = year(Date_Time),
         yday = yday(Date_Time),
         hour = hour(Date_Time)) %>%
  group_by(year, yday, hour) %>%
  summarize(date_time = min(Date_Time),
            temp = mean(temp, na.rm=T),
            wdir = mean(wdir, na.rm=T),
            temp = mean(temp, na.rm=T),
            wspeed = mean(wspeed, na.rm=T),
            gust = mean(gust, na.rm=T),
            RADGL = mean(RADGL, na.rm=T)) %>%
  ungroup()

# plot
weather_hourly %>%
  ggplot(aes(date_time, temp))+
  facet_wrap(~year, scales = "free_x")+
  geom_line(size = 0.3)+
  geom_point(size = 0.5)+
  theme_bw() 





