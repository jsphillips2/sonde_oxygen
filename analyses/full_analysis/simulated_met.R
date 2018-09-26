#==========
#========== Preliminaries
#==========

# load packages
library(tidyverse)
library(rstan)

# import data and model fit
sonde_data = read_csv("analyses/full_analysis/model_fit/input/sonde_prep.csv")
model_fit = read_csv("analyses/full_analysis/model_fit/output/sig_obs10/summary_clean.csv")

# select years
years = sonde_data$year %>% unique

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
  # add variable for hour
  group_by(day) %>%
  arrange(day, index) %>%
  mutate(hour = 1:24) %>%
  ungroup() %>%
  select(day, hour, beta0, alpha, rho, o2, o2_pred) %>%
  # combine with sonde data
  full_join(sonde_data %>% rename(day = D_M)) %>%
  group_by(T_S) %>%
  # calculate process and observation errors
  mutate(proc_err = c(0,o2[2:length(o2)] - o2_pred[1:(length(o2_pred)-1)]),
         obs_err = do - o2) %>%
  filter(year %in% years) %>%
  # add remaining parameters
  mutate(gamma_1 = {model_fit %>% filter(name == "gamma_1")}$middle,
         gamma_2 = {model_fit %>% filter(name == "gamma_2")}$middle,
         k0 = {model_fit %>% filter(name == "k0")}$middle,
         k1 = {model_fit %>% filter(name == "k1")}$middle,
         k2 = 1.7,
         sig_obs = 10,
         temp_ref = 11.99139,
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
      sim_o2[t] = o2_pred[t-1] + .$proc_err[t] + .$obs_err[t-1];
    }
    df = data_frame(T_S = .$T_S,
                    day = .$day,
                    hour = .$hour,
                    sim_o2 = sim_o2,
                    gpp = gpp,
                    er = er,
                    nep = nep)
    return(df)
  }
}





#==========
#========== Simulate data
#==========

# simulate
met_sim = data_prep %>% 
  mutate(
    fix_beta0 = T,
    fix_alpha = F,
    fix_rho = T,
    fix_temp = T,
    fix_par = T,
    beta0 = ifelse(fix_beta0 == T, mean(beta0, na.rm=T), beta0),
    alpha = ifelse(fix_alpha == T, mean(alpha, na.rm=T), alpha),
    rho = ifelse(fix_rho == T, mean(rho, na.rm=T), rho),
    temp = ifelse(fix_temp == T, mean(temp, na.rm=T), temp),
    par = ifelse(fix_par == T, mean(par, na.rm=T), par)
  ) %>% 
  split(.$T_S) %>%
  lapply(function(x){sim_fn(x)}) %>%
  bind_rows() %>%
  full_join(data_prep) %>%
  select(year, yday, hour, gpp, er, nep) %>%
  arrange(year, yday, hour) %>%
  mutate(time = yday + hour/24) %>%
  mutate(type = "alpha") %>%
  bind_rows(data_prep %>% 
              mutate(
                fix_beta0 = F,
                fix_alpha = T,
                fix_rho = T,
                fix_temp = T,
                fix_par = T,
                beta0 = ifelse(fix_beta0 == T, mean(beta0, na.rm=T), beta0),
                alpha = ifelse(fix_alpha == T, mean(alpha, na.rm=T), alpha),
                rho = ifelse(fix_rho == T, mean(rho, na.rm=T), rho),
                temp = ifelse(fix_temp == T, mean(temp, na.rm=T), temp),
                par = ifelse(fix_par == T, mean(par, na.rm=T), par)
              ) %>% 
              split(.$T_S) %>%
              lapply(function(x){sim_fn(x)}) %>%
              bind_rows() %>%
              full_join(data_prep) %>%
              select(year, yday, hour, gpp, er, nep) %>%
              arrange(year, yday, hour) %>%
              mutate(time = yday + hour/24) %>%
              mutate(type = "beta0")
  ) %>%
  group_by(type, year, yday) %>%
  summarize(nep = sum(nep)/1000) %>%
  ungroup %>%
  group_by(type, year) %>%
  mutate(nep = nep - mean(nep, na.rm=T))

met_sim  %>%
  {ggplot(.,aes(yday, nep, color = type))+
      facet_wrap(~year)+
      geom_line()+
      geom_line(data = . %>% filter(is.na(nep) == F), linetype = 2)+
      scale_color_manual(values = c("firebrick","dodgerblue"))+
      # scale_y_continuous(limits=c(4,))+
      theme_base
  }

met_sim %>% 
  group_by(type) %>%
  summarize(sd = sd(gpp, na.rm=T))



