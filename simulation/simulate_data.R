#==========
#========== Preliminaries
#==========

# load packages
library(tidyverse)

# import data and model fit
sonde_data = read_csv("data/sonde_prep.csv")
model_fit = read_csv("main_analysis/model_output/summary_clean.csv")

# selects years
years = c(2012)





#==========
#========== Prepare data for simulation
#==========

test = model_fit %>%
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
         z = 3.3)




#==========
#========== Simulate data
#==========

test2 = test %>%
  split(.$T_S) %>%
  lapply(function(x){
    x %>% {
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
  })

# add simulated data to real data
test3 = test %>% 
  left_join(test2 %>% bind_rows())

# plot

