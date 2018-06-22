#==========
#========== Preliminaries
#==========

# load packages
library(tidyverse)

# import "real" fits
model_fit = read_csv("simulations/fit_16/model_output_16/summary_clean.csv")

# import simulated fits
input_folder = "simulations/sim_16/sim_fits/sim_a"
sim_fits = lapply(list.files(input_folder, full.names = T, pattern = "*.csv"), 
               read_csv) %>%
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
#========== Plot: beta0, alpha, rho
#==========

plot_vars = c("beta0","alpha","rho")
model_fit %>%
  filter(name %in% plot_vars) %>%
  select(name, index, day, middle) %>%
  rename(middle_true = middle)  %>%
  ggplot(aes(day, middle_true))+
  facet_wrap(~name, scales="free_y",nrow=3)+
  geom_line(data = sim_fits %>%
              filter(name %in% plot_vars),
            aes(day, middle, group=rep), size = 0.5, alpha=0.5)+
  geom_line(size = 0.7)+
  theme_bw()





#==========
#========== Plot: GPP, ER, NEP
#==========

plot_vars = c("GPP","ER","NEP")
model_fit %>%
  filter(name %in% plot_vars) %>%
  select(name, index, day, middle) %>%
  rename(middle_true = middle)  %>%
  ggplot(aes(day, middle_true))+
  facet_wrap(~name, scales="free_y",nrow=3)+
  geom_line(data = sim_fits %>%
              filter(name %in% plot_vars),
            aes(day, middle, group=rep), size = 0.5, alpha=0.5)+
  geom_line(size = 0.7)+
  theme_bw()
