#==========
#========== Preliminaries
#==========

# load packages
library(tidyverse)

# import specifications
type = c("fixed_all","fixed_none")
rep = paste0("rep_",c(1:5))

# import data
# observed
sonde_data = read_csv("data/sonde_prep.csv")

# simulated data
sim_data = type %>%
  lapply(function(x){
    rep %>%
      lapply(function(y){
        read_csv(paste0("simulation/simulated_data/",x,"/",y,"/data_export.csv")) %>%
          mutate(type = x,
                 rep = y)
      }) %>% bind_rows()
  }) %>%
  bind_rows() 

# model fit
model_fit = type %>%
  lapply(function(x){
    rep %>%
      lapply(function(y){
        read_csv(paste0("simulation/simulated_output/",x,"/",y,"/summary_clean.csv")) %>%
          mutate(type = x,
                 rep = y)
      }) %>% bind_rows()
  }) %>%
  bind_rows() 

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
#========== Plot
#==========

sim_mean = sim_data %>%
  na.omit() %>%
  group_by(type, year, yday) %>%
  summarize(beta0 = mean(beta0),
            alpha = mean(alpha),
            rho = mean(rho)) 

# beta0
model_fit %>%
  filter(name %in% c("beta0")) %>%
  left_join(sonde_data %>%
              group_by(year, yday) %>%
              summarize(day = unique(D_M))) %>%
  mutate(middle = middle/1000,
         lower16 = lower16/1000,
         upper84 = upper84/1000) %>%
  arrange(type, rep, year,yday) %>%
  ggplot(aes(yday, middle))+
  facet_wrap(~type)+
  geom_line(aes(group = rep), alpha = 0.5)+
  geom_line(data = sim_mean, aes(x = yday, y = beta0/1000), size = 0.7)+
  theme_base

# alpha
model_fit %>%
  filter(name %in% c("alpha")) %>%
  left_join(sonde_data %>%
              group_by(year, yday) %>%
              summarize(day = unique(D_M))) %>%
  mutate(middle = middle/1000,
         lower16 = lower16/1000,
         upper84 = upper84/1000) %>%
  arrange(type, rep, year,yday) %>%
  ggplot(aes(yday, middle))+
  facet_wrap(~type)+
  geom_line(aes(group = rep), alpha = 0.5)+
  geom_line(data = sim_mean, aes(x = yday, y = alpha/1000), size = 0.7)+
  theme_base

# rho
model_fit %>%
  filter(name %in% c("rho")) %>%
  left_join(sonde_data %>%
              group_by(year, yday) %>%
              summarize(day = unique(D_M))) %>%
  mutate(middle = middle/1000,
         lower16 = lower16/1000,
         upper84 = upper84/1000) %>%
  arrange(type, rep, year,yday) %>%
  ggplot(aes(yday, middle))+
  facet_wrap(~type)+
  geom_line(aes(group = rep), alpha = 0.5)+
  geom_line(data = sim_mean, aes(x = yday, y = rho/1000), size = 0.7)+
  theme_base




