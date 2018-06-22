#==========
#========== Preliminaries
#==========

# load packages
library(tidyverse)
library(rstan)

# read data
data = read_csv("data/sonde_prep.csv")
model_fit = read_csv("main_analysis/model_output/summary_clean.csv")




#==========
#========== Subset & Prepare Data
#==========

# subset for desired years
year = 2012
data_sub = data %>%
  filter(year == 2012) %>%
  na.omit() %>%
  mutate(day = D_M) %>%
  # merge with model fit
  left_join(model_fit %>%
              filter(name %in% c("beta")) %>%
              arrange(index) %>%
              group_by(day) %>%
              mutate(hour = 1:24) %>%
              ungroup %>%
              select(name, hour, day, middle) %>%
              spread(name, middle)) %>%
  left_join(model_fit %>%
              select(name,day, middle) %>%
              filter(name %in% c("alpha","rho")) %>%
              spread(name, middle)) %>%
  cbind(model_fit %>%
          select(name, middle) %>%
          filter(name %in% c("gamma_1","gamma_2","sig_proc","k0","k1")) %>%
          spread(name, middle)) %>%
  tbl_df 

# define fixed values and indexing variables
o2_obs = data_sub$do
o2_eq = data_sub$do_eq
light = data_sub$par
temp = data_sub$temp
wspeed = data_sub$wspeed
sch_conv = data_sub$sch_conv
D_M = data_sub$D_M
S = c({data_sub %>%
    group_by(T_S) %>%
    summarize(value = length(T_S))}$value,1) 
K = c({data_sub %>%
    group_by(year) %>%
    summarize(value = length(unique(D_M)))}$value,1)
temp_ref = 11.99139
z = 3.3
sig_obs = 10
k2 = 1.7
N = length(o2_obs)
T_S = length(S)-1 
D = sum(K)-1
Y = length(K)-1
o2_st = c(1,c(cumsum((S)[1:(T_S-1)]) + 1))
dy_st = c(1, c(cumsum((K)[1:(Y-1)]) + 1))





#==========
#========== Simulate Data
#==========

# declare empty variables
o2_s = rep(NA, N)
gpp_s = rep(NA, N)
er_s = rep(NA, N)
nep_s = rep(NA, N)
air_s = rep(NA, N)
o2_pred_s = rep(NA, N)
o2_obs_s = rep(NA, N)

# function for simulations
sim_fun = function(data){
  as.list(data) %>%
    with(
      for (t in 1:T_S){
        o2_s[o2_st[t]] = o2_obs[o2_st[t]];
        o2_obs_s[o2_st[t]] = rnorm(1, o2_s[o2_st[t]], sd = sig_obs);
        for (n in (o2_st[t]+1):(o2_st[t]+S[t])){
          gpp_s[n-1] = beta[n-1]*tanh((alpha[D_M[n-1]]/beta[n-1])*light[n-1]);
          er_s[n-1] = rho[D_M[n-1]]*gamma_2^(temp[n-1] - temp_ref);
          nep_s[n-1] = gpp_s[n-1] - er_s[n-1];
          air_s[n-1] = ((k0 + k1*wspeed[n-1]^k2)/100)*sch_conv[n-1]*(o2_eq[n-1] 
                                                                         - o2_s[n-1]);
          o2_pred_s[n-1] = o2_s[n-1] + (nep_s[n-1] + air_s[n-1])/z;      
          o2_s[n] = rnorm(1, mean = o2_pred_s[n-1], sd = sig_proc);
          o2_obs_s[n] = rnorm(1, o2_s[n], sd = sig_obs)
        }
        return(o2_obs_s);
      })
}

# simulate
o2_obs_s = sim_fun(data_sub)

plot(na.omit(o2_obs_s))

