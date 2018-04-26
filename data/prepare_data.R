#==========
#========== Preliminaries
#==========

# load packages
library(tidyverse)
library(rstan)

# select years desired for analysis
years = c(2013, 2015, 2016, 2017)

# read data
sonde = read_csv("data/sonde_final.csv") %>% 
  # extract desired years
  filter(year %in% years) %>%
  # convert hour from 0:23 to 1:24 (this makes indexing easier later)
  mutate(hour = hour + 1)





#==========
#========== Clean and filter data
#==========

# initial clean
sonde_prep = sonde %>% 
  # drop and rename/define pars
  select(-par) %>%
  rename(temp = wat_temp, par = par_int) %>%
  mutate(sch_conv = (sch_o2/600)^(-0.5)) %>%
  select(year,month,yday,hour,par,wspeed,temp,do,do_eq,sch_conv)  %>%
  # omit na's and select days with at least 24 hours
  na.omit %>%
  group_by(year, yday) %>%
  mutate(nhour = length(hour)) %>% 
  filter(nhour>23) %>%
  # join with year x yday x hour to fill in na's
  full_join(sonde %>% expand(year,yday,hour)) %>%
  mutate(do = 1000*do, do_eq = 1000*do_eq)

#  examine plot to identify sparse time periods
sonde_prep %>% 
  group_by(year,yday) %>% 
  summarize(do = mean(do, na.rm=T)) %>%
  ggplot(aes(yday, do))+
  facet_wrap(~year)+
  geom_line()+
  theme_classic()

# identify day with min values in 2013
sonde_prep %>% 
  group_by(year,yday) %>% 
  summarize(do = mean(do, na.rm=T), temp=mean(temp)) %>% 
  filter(year==2013) %>%
  arrange(do)

# plot hourly curves for those days
sonde_prep %>%
  filter(year==2013, yday %in% c(211:215)) %>%
  group_by(yday) %>%
  ggplot(aes(hour, do))+
  facet_wrap(~yday, scale="free_y")+
  geom_point()+
  theme_classic()

# these days show reasonable diel curves
# they are very warm, so that could explain low value

# late 2015 looks too sparse; filter for days < 184
# truncate at day 230
sonde_prep2 = sonde_prep %>%
  filter(!(year==2015 & yday > 184), yday < 230) %>%
  arrange(year, yday, hour)





#==========
#========== Prepare for export and check
#==========

# calculate variable T_S to map observations to time series
# calculate variable D_M to map observations to days
# omit NA's 
# replace 0 PAR with minimum non-0 PAr
# convert to data frame
sonde_prep3 = sonde_prep2 %>%
  group_by(year) %>%
  mutate(j = ifelse(is.na(do)==T, 1, 0), 
         k = c(1,abs(diff(j)))) %>% 
  filter(is.na(do)==F) %>%
  mutate(T_S = cumsum(k)) %>%
  ungroup() %>%
  mutate(yearT_S = paste(year, T_S),
         T_S = as.numeric(as.factor(yearT_S)),
         year_day = paste(year, yday),
         D_M = as.numeric(as.factor(year_day)),
         par = ifelse(par==0, min(par[which(par>0)]), par)) %>%
  select(-j, -k ,-yearT_S, -year_day) %>%
  as.data.frame()

# check T_S
sonde_prep2 %>% 
  expand(year,month,yday,hour) %>%
  full_join(sonde_prep4) %>%
  arrange(year,yday) %>%
  mutate(time = yday + hour/24) %>%
  ggplot(aes(time, do, color=factor(T_S)))+
  facet_wrap(~year)+
  geom_line()+
  theme_bw()

# check D_M
sonde_prep2 %>% 
  expand(year,month,yday,hour) %>%
  full_join(sonde_prep4) %>%
  arrange(year,yday) %>%
  mutate(time = yday + hour/24) %>%
  filter(D_M %in% (50 + c(1:10))) %>%
  ggplot(aes(time, do, color=factor(D_M)))+
  facet_wrap(~year, scale = "free_x", nrow=2)+
  geom_line()+
  theme_bw()

# check D_M day extraction
sonde_prep3 %>% 
  ungroup() %>%
  filter(D_M == sonde_prep3$D_M[1520])

# export prepared data
# write_csv(sonde_prep2, "data/sonde_prep.csv")




#==========
#========== Package data 
#==========

# define variables in evnironment 
o2_obs = sonde_prep3$do
o2_eq = sonde_prep3$do_eq
light = sonde_prep3$par
temp = sonde_prep3$temp
wspeed = sonde_prep3$wspeed
sch_conv = sonde_prep3$sch_conv
D_M = sonde_prep3$D_M
S = {sonde_prep3 %>%
    group_by(T_S) %>%
    summarize(value = length(T_S))}$value 
K = {sonde_prep3 %>%
    group_by(year) %>%
    summarize(value = length(unique(D_M)))}$value
temp_ref = mean(sonde_prep3$temp)
z = 3.3
sig_obs = 10
k0 = 2.07
k1 = 0.215
k3 = 1.7
N = length(o2_obs)
T_S = length(S) 
D = sum(K)
Y = length(K)
o2_st = c(1, cumsum((S)[1:(T_S-1)]) + 1)
dy_st = c(1, cumsum((K)[1:(Y-1)]) + 1) 

# export as .R
# stan_rdump(c("D_M","S","K","N","D","Y","T_S","o2_st","dy_st",
#              "o2_obs","o2_eq","light","temp","temp_ref", "wspeed",
#              "sch_conv","z","sig_obs","k0","k1","k3"), file="data/sonde_list.R")


