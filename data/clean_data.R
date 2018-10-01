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

# profiles
profile = read_csv("data/myvatn_profile.csv")





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

# plot phycoyanin time series by year
sonde_trim %>%
  mutate(year = year(Date_Time)) %>%
  filter(is.na(year) == F) %>%
  ggplot(aes(Date_Time, PCYV))+
  facet_wrap(~year, scales = "free_x")+
  geom_line(size = 0.3)+
  geom_point(size = 0.5)+
  theme_bw()

# plot turb time series by year
sonde_trim %>%
  mutate(year = year(Date_Time)) %>%
  filter(is.na(year) == F) %>%
  ggplot(aes(Date_Time, TurbSC))+
  facet_wrap(~year, scales = "free_x")+
  geom_line(size = 0.3)+
  geom_point(size = 0.5)+
  theme_bw()

# omit LDO < 7
# omit Temp > 18
# omit PCYV > 1
# omit Turb > 200
# remove values that deviate from daily means by greater than 2 standard deviations
# results in omission of 5.6% of the data
sonde_clean = sonde_trim %>%
  filter(LDO > 5, Temp < 18, TurbSC < 200) %>%
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
            turb = mean(TurbSC, na.rm=T),
            pcyv = mean(PCYV, na.rm=T),
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

# combine HOBO and weather data (only keeping necessary variables)
# divide lux by 54 to roughly convert to PAR
# multiply RADGL by 4.57 to roughly convert ot PAR
# *Sager and McFarlane--Radiation in NCERA: Controlled Environments
hobo_weather = hobo_hourly %>%
  mutate(par_hobo = lux/54) %>%
  select(year, yday, hour, par_hobo) %>%
  full_join(weather_hourly %>%
              filter(!(year %in% c(2010, 2011, 2014))) %>%
              mutate(par_weath = 4.57*RADGL) %>%
              select(year, yday, hour, wspeed, par_weath))

# plot par_hobo vs. par_weath
hobo_weather %>%
  ggplot(aes(par_weath, par_hobo))+
  geom_point(size = 0.5, alpha = 0.5)+
  geom_line(data = data_frame(par_weath = 0:6000, par_hobo = 0:6000), size = 0.7)+
  scale_y_continuous(trans="log1p")+
  scale_x_continuous(trans="log1p")+
  theme_bw()

# plot PARs by hour
hobo_weather %>%
  filter(year == 2013, yday %in% c(170:180)) %>%
  gather(var, value, par_hobo, par_weath) %>%
  mutate(time = yday + hour/24) %>%
  ggplot(aes(time, value, color = var))+
  geom_line()+
  theme_bw()

# PAR from the weather station has slightly higher peaks than from the HOBO
# This makes sense, as RADGL probably intergrates over a wider range of wavelengths
# Lux from the HOBO (using visible light) corresponds more closely to real PAR
# So, it makes sense to use it for the analysis
# However, HOBO data are missing from 2012, so use RADGL for 2012
# Previously I tried doing a regression to roughly convert between the two, but this seems to overcorrect

hobo_weather = hobo_weather %>%
  mutate(par = ifelse(is.na(par_hobo)==F, par_hobo, par_weath))






#==========
#========== Combine sonde, HOBO, and weather data
#==========

# join sonde, hobo, and weather data
# fill in NAs 
sonde_full = sonde_hourly %>%
  select(year, yday, hour, temp, turb, pcyv, o2_sat, do) %>%
  left_join(hobo_weather %>%
               select(year, yday, hour, par, wspeed)) %>%
  full_join(sonde_hourly %>% 
              expand(year, yday, hour)) %>%
  arrange(year, yday, hour)

# plot
sonde_full %>%
  mutate(time = yday + hour/24) %>%
  {ggplot(.,aes(time, do))+
      facet_wrap(~year)+
      geom_line(data = . %>% filter(is.na(do)==F), size = 0.3)+
      geom_point(size = 0.5)+
      theme_bw()}





#==========
#========== Light 
#==========

# extract relevant data
# replace 0's with the lowest observed value (for upcoming log transformation)
profile_clean = profile %>%
  filter(sta %in% c(3,33)) %>%
  mutate(year = year(sampledate), 
         yday = yday(sampledate),
         hour = hour(sampletime),
         light = ifelse(light>0,light,{profile %>% filter(light >0)}$light %>% min))

# define function for calculating light extinction
ext_fun = function(x){
  if (nrow(x) < 2) {
    z <- NA
  } else {
    m = lm(log(light) ~ sampledepth, data = x)
    z1 = m$coef[['(Intercept)']]
    z2 = m$coef[['sampledepth']]
  }
  data_frame(sampledate = x$sampledate[1],
             year = x$year[1], 
             yday = x$yday[1], 
             hour = x$hour[1], 
             par0 = exp(z1),
             ext = -z2)
}

# calculate extinction
ext_data = profile_clean %>% 
  select(sampledate, year, yday, hour, sampledepth, light) %>%
  na.omit() %>%
  split(.$sampledate, .$hour) %>% 
  map_df(~ ext_fun(.x)) %>%
  # filter negative values of ext
  filter(ext > 0) %>%
  # add PAR (including for times missing in sonde_full)
  left_join(hobo_weather %>%
              select(year, yday, hour, par) %>%
              na.omit()) %>%
  # add daily averages of turbidity
  left_join(sonde_full %>%
              group_by(year, yday) %>%
              summarize(pcyv = mean(pcyv, na.rm=T),
                        turb = mean(turb, na.rm=T)))

# plot water extinction coefficient vs turbidity
ext_data %>%
  na.omit() %>%
  ggplot(aes(turb, ext))+
  geom_point()+
  theme_bw()

# fit linear model to predict extinction from turbidity
ext_mod1 = lm(ext ~ turb, data = ext_data)
ext_mod1 %>% summary

# plot par0 vs par
ext_data %>%
  na.omit() %>%
  ggplot(aes(par, par0))+
  geom_point()+
  geom_abline(intercept = 0, slope = 1)+
  theme_bw()

# predict surface light from HOBO light
# for intercept through 0 (since surface PAR in water must be 0 if HOBO PAR is 0)
par0_mod = lm(par0 ~ par - 1, data = ext_data)
summary(par0_mod)

# add model predictions for light extinction and surface light
sonde_full_b = sonde_full %>%
  mutate(ext = predict(ext_mod2, newdata = .),
         par0 = predict(par0_mod, newdata = .),
         # calculate light at center of water column
         # assume light varies with depth as PAR(z) = PAR_0*exp(-ext*z)
         # intergrating from the top (z=0) to the bottom (z=zmax) of the water column:
         # PAR_integrated = (PAR_0 - PAR_0*exp(-ext*zmax))/(ext*zmax)
         # Myvatn at ST33 has zmax of 3.3m
         par_int = (par0 - par0*exp(-ext*3.3))/(ext*3.3))

# check ext prediction
sonde_full_b %>%
  ggplot(aes(turb, ext))+
  geom_line(size = 0.3)+
  geom_point(size = 0.5)+
  theme_bw()

# check par prediction
sonde_full_b %>%
  ggplot(aes(par, par0))+
  geom_line(size = 0.3)+
  geom_abline(intercept = 0, slope = 1)+
  geom_point(size = 0.5)+
  theme_bw()

# check par_int
sonde_full_b %>%
  ggplot(aes(par0, par_int))+
  geom_line(size = 0.3)+
  geom_point(size = 0.5)+
  theme_bw()





#==========
#========== Gas calculations 
#==========

sonde_full_c = sonde_full_b %>%
  mutate(
    # tempreature in Kelvin
    temp_k = temp + 273.15,
    # Schmidt number and associated conversion from CO2 to O2
    # based on Wanninkhoff 1992; Holtgrieve et al 2010; Staehr
    sch_o2 = 1800.6 + 120.10*temp + 3.7818*temp^2 - 0.047608*temp^3,
    sch_conv = (sch_o2/600)^(-0.5),
    # O2 solubility in mL/L
    # based on Weiss 1970
    do_sol = exp(-173.4292 + 249.6339*(100/temp_k) + 143.3483*log(temp_k/100) - 21.8492*(temp_k/100)),
    # convert to mg/L of DO
    # use ideal gas law, solved for the number of moles n = P*V/(R*T)
    # use pressure in kPA corrected for Myvatn's elevation of ~300m = 98 kPA
    # use R = 8.3144598 L*kPa/(K*mol)  
    # use O2 molar mass of 2*32; multiply by 1000 to convert from g to mg 
    do_eq = 32*(98*do_sol)/(8.3144598*temp_k),
    # calculate equilibrium DO from back-calculation from sonde o2_sat for comparison
    do_eq2 = 100*do/o2_sat
  )

# compare do_eq vs do_eq2
sonde_full_c %>% 
  select(year,yday,hour,do_eq,do_eq2) %>%
  gather(type, do_eq, do_eq, do_eq2) %>%
  mutate(time = yday + hour/24) %>%
  ggplot(aes(time, do_eq, color=type))+
  facet_wrap(~year)+
  geom_line()+
  geom_hline(yintercept = 10,alpha=0.5,size=0.5)+
  scale_color_manual(values=c("black","firebrick"))+
  theme_bw()

# strongly correlated and generally quite close
# except for 2015, where my calculation is much higher
sonde_full_c %>% 
  mutate(time = yday + hour/24) %>%
  ggplot(aes(time, temp))+
  facet_wrap(~year)+
  geom_line()+
  geom_hline(yintercept = 10,alpha=0.5,size=0.5)+
  theme_bw()

# 2015 was very cold according to the sonde

# this is corroborated by the weather station data
weather_hourly %>%
  filter(year %in% c(2012:2013, 2015:2018),
         yday %in% c(150:230)) %>%
mutate(time = yday + hour/24) %>%
  ggplot(aes(time, temp))+
  facet_wrap(~year)+
  geom_line()+
  geom_hline(yintercept = 10,alpha=0.5,size=0.5)+
  theme_bw()

# use my calculation, given that I know what went into it and it makes sense given temperature





#==========
#========== Export data
#==========

# sonde_full_c %>%
#   select(-temp_k, -do_eq2) %>%
#   write_csv("data/sonde_final.csv")














