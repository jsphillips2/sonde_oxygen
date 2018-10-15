#==========
#========== Preliminaries
#==========

# load packages
library(tidyverse)
library(rstan)

# import data and model fit
sonde_data = read_csv("analyses/model_fit/input/sonde_prep.csv")
model_fit = read_csv("analyses/model_fit/output/summary_clean.csv")
data = read_rdump("analyses/model_fit/input/sonde_list.R")

# base theme
theme_base = theme_bw()+
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.margin = margin(0,0,0,0),
        text = element_text(size=12),
        strip.text = element_text(size=10),
        legend.text = element_text(size=10),
        axis.text=element_text(size=10, color="black"),
        axis.title.y=element_text(margin=margin(0,15,0,0)),
        axis.title.x=element_text(margin=margin(15,0,0,0)))




#==========
#========== Prepare for simulation
#==========

# prepare data
data_prep = 
  # select estimated and observed oxygen
  model_fit %>%
  filter(name %in% c("o2_pred","o2")) %>%
  select(-lower16, -upper84) %>%
  spread(name, middle) %>%
  # combine with time-varying parameters
  left_join(model_fit %>%
              filter(name %in% c("beta0","alpha","rho")) %>%
              select(-index, -lower16, -upper84) %>%
              spread(name, middle)) %>%
  # combine with sonde data
  left_join(sonde_data %>%
              # filter(is.na(unique_day)==F) %>%
              select(year, yday, index, hour,  unique_series, unique_day, do, do_eq, 
                     par_int, temp, wspeed, sch_conv)) %>%
  mutate(do = 1000*do, do_eq = 1000*do_eq) 

# define fixed parameters
fixed_vals = data_frame(gamma_1 = {model_fit %>% filter(name == "gamma_1")}$middle,
                        gamma_2 = {model_fit %>% filter(name == "gamma_2")}$middle,
                        k0 = data$k0,
                        k1 = data$k1,
                        k2 = data$k2,
                        z = data$z,
                        temp_ref = 12)

# define function to simulate
sim_fn = function(x, y){
  x %>% 
  {
    sim_o2 = c(.$o2[1], rep(NA, nrow(x) - 1))
    beta = rep(NA, nrow(x))
    gpp = rep(NA, nrow(x))
    er = rep(NA, nrow(x))
    nep = rep(NA, nrow(x))
    air = rep(NA, nrow(x))
    o2_pred = rep(NA, nrow(x))
    for(t in 2:nrow(.)){
      beta[t-1] = .$beta0[t-1]*y$gamma_1^(.$temp[t-1] - y$temp_ref);
      gpp[t-1] = beta[t-1]*tanh((.$alpha[t-1]/beta[t-1])*.$par_int[t-1]);
      er[t-1] = .$rho[t-1]*y$gamma_2^(.$temp[t-1] - y$temp_ref);
      nep[t-1] = gpp[t-1] - er[t-1];
      air[t-1] = ((y$k0 + y$k1*.$wspeed[t-1]^y$k2)/100)*.$sch_conv[t-1]*(.$do_eq[t-1] - sim_o2[t-1]);
      o2_pred[t-1] = sim_o2[t-1] + (nep[t-1] + air[t-1])/y$z;
      sim_o2[t] = o2_pred[t-1] + .$proc_err[t-1];
    }
    df = data_frame(unique_series = .$unique_series,
                    year = .$year,
                    day = .$day,
                    hour = .$hour,
                    sim_o2 = sim_o2)
    return(df)
  }
}





#==========
#========== Simulate data
#==========

# specify which variables should be 'fixed' to their mean
data_prep2 = data_prep %>% 
  group_by(year) %>%
  mutate(
    fix_beta0 = T,
    fix_alpha = T,
    fix_rho = T,
    beta0 = ifelse(fix_beta0 == T, mean(beta0, na.rm=T), beta0),
    alpha = ifelse(fix_alpha == T, mean(alpha, na.rm=T), alpha),
    rho = ifelse(fix_rho == T, mean(rho, na.rm=T), rho)) %>%
  arrange(year, yday, hour) %>%
  group_by(unique_series) %>%
  # calculate process and observation errors
  mutate(proc_err = c(0,o2[2:length(o2)] - o2_pred[1:(length(o2_pred)-1)])) %>%
  ungroup() 

# extract simulation type from specifications
type_data = data_prep2 %>% 
  expand(fix_beta0, fix_alpha, fix_rho) %>%
  gather(var, value, fix_beta0, fix_alpha, fix_rho) %>%
  mutate(name = str_split(var,"_") %>% map_chr(~as.character(.x[2]))) 

# create name for type
type = {type_data %>%
    filter(value == T)}$name %>%
  paste0(collapse = "", sep = "_") %>%
  paste0("fixed")

# simulate
data_sim = data_prep2  %>%
  # split(.$unique_series) %>%
  split(.$year) %>%
  lapply(function(x) {sim_fn(x, fixed_vals)}) %>%
  bind_rows() %>%
  mutate(sim_o2 = sim_o2/1000) %>%
  full_join(sonde_data %>%
              rename(day = unique_day))


# plot and compare to actual data
data_sim %>%
  mutate(time = yday + hour/24) %>%
  select(year, time, do, sim_o2) %>%
  gather(var, value, do, sim_o2) %>%
  ggplot(aes(time, value, color = var))+
  facet_wrap(~year)+
  geom_line(alpha = 0.7, size = 0.7)+
  scale_color_manual(values=c("firebrick","dodgerblue"))+
  scale_y_continuous(limits=c(5,17))+
  ggtitle(type)+
  theme_base







#==========
#========== Prepare for data analysis
#==========

# prepare data
sonde_prep = data_sim %>% 
  select(-do) %>%
  rename(do = sim_o2) %>%
  arrange(year, yday, hour) %>%
  # for each year, create identifier for uninterrupted stretches of observations
  group_by(year) %>%
  mutate(i = ifelse(is.na(do)==T, 1, 0), 
         j = c(1,abs(diff(i)))) %>% 
  filter(is.na(do)==F, is.na(wspeed)==F) %>%
  mutate(series = cumsum(j)) %>% 
  ungroup() %>%
  # create unique index for each series
  # remove series with fewer than 10 observations
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
sonde_check = sonde_prep %>% 
  expand(year,yday,hour) %>%
  full_join(sonde_prep) %>%
  arrange(year,yday)

# check unique_series
# sonde_check %>%
#   mutate(time = yday + hour/24) %>%
#   ggplot(aes(time, do, color=factor(unique_series)))+
#   facet_wrap(~year)+
#   geom_line()+
#   scale_color_discrete(guide = F)+
#   theme_bw()

# check unique_days
# sonde_check %>%
#   mutate(time = yday + hour/24) %>%
#   filter(unique_day %in% (35 + c(1:10))) %>%
#   ggplot(aes(time, do, color=factor(unique_day)))+
#   facet_wrap(~year, scale = "free_x", nrow=2)+
#   geom_line()+
#   theme_bw()


# export  data
export_file = type
sonde_check %>%
  write_csv(paste0("simulation/input/",
                   export_file,"/data_export.csv"))
type_data %>%
  write_csv(paste0("simulation/input/",
                   export_file,"/type_data.csv"))



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
stan_rdump(c("o2_obs","o2_eq","light","temp","wspeed","sch_conv","map_days","obs_per_series","days_per_year",
             "obs_per_day", "z","k0","k1","k2","n_obs","n_series","n_days","n_years"),
           file=paste0("simulation/input/",export_file,"/data_list.R"))


