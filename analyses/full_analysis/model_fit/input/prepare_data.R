#==========
#========== Preliminaries
#==========

# load packages
library(tidyverse)
library(rstan)

# select years desired for analysis
years = c(2012:2013,2015:2018)

# read data
sonde = read_csv("data/sonde_final.csv") %>% 
  # extract desired years
  filter(year %in% years) %>%
  # convert hour from 0:23 to 1:24 (this makes indexing easier later)
  mutate(hour = hour + 1) 

# examine
sonde %>%
  mutate(time = yday + hour/24) %>%
  {ggplot(.,aes(time, do))+
      facet_wrap(~year)+
      geom_line(data = . %>% filter(is.na(do)==F), size = 0.3)+
      geom_point(size = 0.5)+
      theme_bw()}

# omit final days of 2015 due to large gap
sonde %>%
  filter(year == 2015, yday < 190 & yday > 180) %>%
  mutate(time = yday + hour/24) %>%
  {ggplot(.,aes(time, do))+
      geom_line(data = . %>% filter(is.na(do)==F), size = 0.3)+
      geom_point(size = 0.5)+
      theme_bw()}

sonde_2 = sonde %>% 
  filter(!(year == 2015 & yday > 185))




#==========
#========== Prepare for data analysis
#==========

# prepare data
sonde_prep = sonde_2 %>%
  arrange(year, yday, hour) %>%
  # for each year, create identifier for uninterrupted stretches of observations
  group_by(year) %>%
  mutate(i = ifelse(is.na(do)==T, 1, 0), 
         j = c(1,abs(diff(i)))) %>% 
  filter(is.na(do)==F, is.na(wspeed)==F) %>%
  mutate(series = cumsum(j)) %>% 
  ungroup() %>%
  # create unique index for each series
  # remove series with fewer than 24 observations
  mutate(unique_series = year + series/length(unique(series))) %>%
  group_by(unique_series) %>%
  mutate(series_length = length(unique_series)) %>%
  filter(series_length > 23) %>%
  ungroup() %>%
  # recreate series index and make unique index for days
  # create index for observations (for joining later)
  # replace 0 par_int with smallest non-zero value
  mutate(unique_series = as.factor(unique_series) %>% as.numeric(),
         unique_day = paste(year, yday) %>% as.factor() %>% as.numeric(),
         index = 1:length(do),
         par_int = ifelse(par_int==0, min(par[which(par_int>0)]), par_int)
         ) %>%
  select(-i, -j) 

# return missing observations for check
sonde_check = sonde_2 %>% 
  expand(year,yday,hour) %>%
  full_join(sonde_prep) %>%
  arrange(year,yday)

# check unique_series
sonde_check %>%
  mutate(time = yday + hour/24) %>%
  ggplot(aes(time, do, color=factor(unique_series)))+
  facet_wrap(~year)+
  geom_line()+
  scale_color_discrete(guide = F)+
  theme_bw()

# check unique_days
sonde_check %>%
  mutate(time = yday + hour/24) %>%
  filter(unique_day %in% (50 + c(1:10))) %>%
  ggplot(aes(time, do, color=factor(unique_day)))+
  facet_wrap(~year, scale = "free_x", nrow=2)+
  geom_line()+
  theme_bw()

# export prepared data
# sonde_check %>%
#   write_csv("analyses/full_analysis/model_fit/input/sonde_prep.csv")





#==========
#========== Package data 
#==========

# define variables in evnironment 
o2_obs = 1000*sonde_prep$do # convert to mg m^-3
o2_eq = 1000*sonde_prep$do_eq # convert to mg m^-3
light = sonde_prep$par_int
temp = sonde_prep$temp
wspeed = sonde_prep$wspeed
sch_conv = sonde_prep$sch_conv
map_days = sonde_prep$unique_day
days_per_year = c({sonde_prep %>%
    group_by(year) %>%
    summarize(value = length(unique(unique_day)))}$value)
obs_per_series = c({sonde_prep %>%
    group_by(unique_series) %>%
    summarize(value = length(unique_series))}$value) 
obs_per_day = c({sonde_prep %>%
    group_by(unique_day) %>%
    summarize(value = length(unique_day))}$value) 
z = 3.3
k0 = 2.07
k1 = 0.215
k2 = 1.7
n_obs = length(o2_obs)
n_series = length(obs_per_series) 
n_days = sum(days_per_year)
n_years = length(days_per_year)

# export as .R
# stan_rdump(c("o2_obs","o2_eq","light","temp","wspeed","sch_conv","map_days","obs_per_series","days_per_year",
#              "obs_per_day", "z","k0","k1","k2","n_obs","n_series","n_days","n_years"),
#            file="analyses/full_analysis/model_fit/input/sonde_list.R")


