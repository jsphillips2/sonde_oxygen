#==========
#========== Preliminaries
#==========

# load packages
library(tidyverse)
library(rstan)

# import data and model fit
sonde_data = read_csv("data/sonde_prep.csv")
model_fit = read_csv("main_analysis/model_output/summary_clean.csv")

# selects years
years = c(2012)

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
data_prep = model_fit %>%
  filter(name %in% c("o2_pred","o2")) %>%
  select(-lower16, -upper84) %>%
  spread(name, middle) %>%
  left_join(model_fit %>%
              filter(name %in% c("beta0","alpha","rho")) %>%
              select(-index, -lower16, -upper84) %>%
              spread(name, middle)) %>%
  group_by(day) %>%
  mutate(hour = 1:24) %>%
  ungroup() %>%
  select(day, hour, beta0, alpha, rho, o2, o2_pred) %>%
  full_join(sonde_data %>% rename(day = D_M)) %>%
  group_by(T_S) %>%
  mutate(proc_err = c(NA,o2[2:length(o2)] - o2_pred[1:(length(o2_pred)-1)])) %>%
  filter(year %in% years, 
         yday >= min({sonde_data %>% 
             filter(year %in% years, 
                    is.na(par)==F)}$yday)) %>%
  mutate(gamma_1 = {model_fit %>% filter(name == "gamma_1")}$middle,
         gamma_2 = {model_fit %>% filter(name == "gamma_2")}$middle,
         k0 = {model_fit %>% filter(name == "k0")}$middle,
         k1 = {model_fit %>% filter(name == "k1")}$middle,
         k2 = 1.7,
         sig_obs = 10,
         temp_ref = 12,
         z = 3.3) %>%
  ungroup()

# define function to simulate
sim_fn = function(x){
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
      beta[t-1] = .$beta0[t-1]*.$gamma_1^(.$temp[t-1] - .$temp_ref);
      gpp[t-1] = beta[t-1]*tanh((.$alpha[t-1]/beta[t-1])*.$par[t-1]);
      er[t-1] = .$rho[t-1]*.$gamma_2^(.$temp[t-1] - .$temp_ref);
      nep[t-1] = gpp[t-1] - er[t-1];
      air[t-1] = ((.$k0 + .$k1*.$wspeed[t-1]^.$k2)/100)*.$sch_conv[t-1]*(.$do_eq[t-1] - sim_o2[t-1]);
      o2_pred[t-1] = sim_o2[t-1] + (nep[t-1] + air[t-1])/.$z;
      sim_o2[t] = o2_pred[t-1] + .$proc_err[t];
    }
    df = data_frame(T_S = .$T_S,
                    day = .$day,
                    hour = .$hour,
                    sim_o2 = sim_o2)
    return(df)
  }
}


#==========
#========== Simulate data
#==========

# set simulation type and seed
type = "fixed_all"
seed = 5

# simulate data 
set.seed(seed)
data_prep2 = data_prep %>% 
  mutate(
    # comment out lines to change fixed values of simulation
    beta0 = mean(beta0, na.rm=T),
    alpha = mean(alpha, na.rm=T),
    rho = mean(rho, na.rm=T),
    # generate process errors (sample from 'real')
    proc_err = sample(na.omit(proc_err), 
                      size = length(proc_err),
                      replace = T)
  )
data_sim = data_prep2 %>%
  full_join(data_prep2  %>%
              split(.$T_S) %>%
              lapply(function(x) 
              {sim_fn(x)}) %>%
              bind_rows()) %>%
  mutate(type = type, seed = seed)

# plot and compare to originally predicted
data_sim %>%
  mutate(time = yday + hour/24) %>%
  select(time, o2_pred, sim_o2) %>%
  gather(var, value, o2_pred, sim_o2) %>%
  ggplot(aes(time, value/1000, color = var))+
  geom_line(alpha = 0.7, size = 0.7)+
  scale_color_manual(values=c("firebrick","dodgerblue"))+
  theme_base





#==========
#========== Export Data
#==========

# select vars to export
data_exp = data_sim %>% 
  select(-do) %>%
  rename(do = sim_o2)

# calculate variable T_S to map observations to time series
# calculate variable D_M to map observations to days
# omit NA's 
# replace 0 PAR with minimum non-0 PAR
# convert to data frame
data_exp2 = data_exp %>%
  arrange(year, yday, hour) %>%
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
data_exp %>% 
  expand(year,month,yday,hour) %>%
  full_join(data_exp2) %>%
  arrange(year,yday) %>%
  mutate(time = yday + hour/24) %>%
  ggplot(aes(time, do, color=factor(T_S)))+
  geom_line()+
  theme_bw()

# check D_M
data_exp %>% 
  expand(year,month,yday,hour) %>%
  full_join(data_exp2) %>%
  arrange(year,yday) %>%
  mutate(time = yday + hour/24) %>%
  filter(D_M %in% (10 + c(1:10))) %>%
  ggplot(aes(time, do, color=factor(D_M)))+
  geom_line()+
  theme_bw()

# check D_M day extraction
data_exp2 %>% 
  ungroup() %>%
  filter(D_M == data_exp2$D_M[1000])

# export prepared data
export_file = paste0(type,"/rep_",seed)
data_exp %>%
  left_join(data_exp2 %>% select(year, month, yday, hour, T_S, D_M)) %>%
write_csv(paste0("simulation/simulated_data/",export_file,"/data_export.csv"))





#==========
#========== Package data 
#==========

# define variables in evnironment 
# add noise to o2
sig_obs = 10
o2_obs = data_exp2$do + rnorm(data_exp2$do , mean = 0, sd = sig_obs)
o2_eq = data_exp2$do_eq
light = data_exp2$par
temp = data_exp2$temp
wspeed = data_exp2$wspeed
sch_conv = data_exp2$sch_conv
D_M = data_exp2$D_M
S = c({data_exp2 %>%
    group_by(T_S) %>%
    summarize(value = length(T_S))}$value,1) 
K = c({data_exp2 %>%
    group_by(year) %>%
    summarize(value = length(unique(D_M)))}$value,1)
temp_ref = 11.99139
z = 3.3
k2 = 1.7
N = length(o2_obs)
T_S = length(S)-1 
D = sum(K)-1
Y = length(K)-1
o2_st = c(1, if(T_S < 2) 1 else c(cumsum((S)[1:(T_S-1)]) + 1, 1))
dy_st = c(1, if(Y < 2) 1 else c(cumsum((K)[1:(Y-1)]) + 1, 1))

# export as .R
# stan_rdump(c("D_M","S","K","N","D","Y","T_S","o2_st","dy_st",
#              "o2_obs","o2_eq","light","temp","temp_ref", "wspeed",
#              "sch_conv","z","sig_obs","k2"), file=paste0("simulation/simulated_data/",export_file,"/data_list.R"))

