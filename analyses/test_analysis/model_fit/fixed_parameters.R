#==========
#========== Preliminaries
#==========

# load packages
library(tidyverse)
library(GGally)
library(truncnorm)

# import data and model fit
input_dir = "analyses/test_analysis/model_fit/"
priors = read_csv(paste0(input_dir,"input/priors.csv"))
params_full = read_csv(paste0(input_dir,
                              "output/sig_obs10/fixed_pars_full.csv"))
model_fit = read_csv(paste0(input_dir,
                            "output/sig_obs10/summary_clean.csv"))

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
#========== Chains
#==========

params_full %>%
  gather(par, value, -chain, -step) %>%
  filter(par != "lp__") %>%
  ggplot(aes(step, value, linetype=factor(chain)))+
  facet_wrap(~par, scales="free_y")+
  geom_line(alpha=0.5, show.legend = F)+
  scale_y_continuous("Value")+
  scale_x_continuous("Iteration", breaks = c(0, 500, 1000))+
  theme_base+
  theme(legend.position = c(0.8,0.15))





#==========
#========== Pairs Plot
#==========

params_full %>% 
  filter(chain == 1) %>%
  select(-chain, -step, -lp__) %>% 
  ggpairs()+
  theme_base






#==========
#========== Gas Exchange
#==========

# set pars to plot
pars = c("k0", "k1")

# generate prior densities
set.seed(1)
prior_dens = priors %>%
  filter(name %in% paste0(pars, c("_prior"))) %>%
  split(.$name) %>%
  {lapply(1:length(.), function(y){
    data_frame(name = str_split(names(.)[y], "_prior") %>% 
                 map_chr(~as.character(.x[1])),
               value = rtruncnorm(n = 10000, a = 0, b = Inf,
                                  mean = .[[y]]$mean, sd = .[[y]]$sd),
               type = "Prior"
    )
  }
  )} %>%
  bind_rows()

# posterior densities
post_dens = params_full %>%
  select("chain", "step", pars) %>%
  gather(name, value, pars) %>%
  select(-chain, - step) %>%
  mutate(type = "Posterior")

# combine densities
comb_dens = bind_rows(prior_dens, post_dens)

# plot
comb_dens %>%
  ggplot(aes())+
  facet_wrap(~name, nrow = 2, scales = "free")+
  stat_density(aes(value, linetype = type),
               geom = "line", size = 0.4, position="identity")+
  geom_rect(data = model_fit %>%
              filter(name %in% pars),
            aes(xmin = lower16, xmax = upper84, ymin = 0, ymax = Inf),
            inherit.aes = F, alpha = 0.2)+
  geom_vline(data = model_fit %>%
               filter(name %in% pars),
             aes(xintercept = middle),
             size = 0.8)+
  scale_y_continuous("Probability Density")+
  scale_x_continuous("Value")+
  scale_linetype_manual("",values=c(1,2))+
  theme_base+
  theme(legend.position = c(0.8,0.916))

# compare sd of prior and posterior
comb_dens %>%
  group_by(name, type) %>%
  summarize(sd = sd(value)) %>%
  spread(type, sd) %>%
  mutate(ratio = Posterior/Prior)





#==========
#========== Figure 4: Process Error
#==========

# set pars to plot
pars = c("sig_proc")

# generate prior densities
set.seed(1)
prior_dens = priors %>%
  filter(name %in% paste0(pars, c("_prior"))) %>%
  split(.$name) %>%
  {lapply(1:length(.), function(y){
    data_frame(name = str_split(names(.)[y], "_prior") %>% 
                 map_chr(~as.character(.x[1])),
               value = rtruncnorm(n = 10000, a = 0, b = Inf,
                                  mean = .[[y]]$mean, sd = .[[y]]$sd),
               type = "Prior"
    )
  }
  )} %>%
  bind_rows()

# posterior densities
post_dens = params_full %>%
  select("chain", "step", pars) %>%
  gather(name, value, pars) %>%
  select(-chain, - step) %>%
  mutate(type = "Posterior")

# combine densities
comb_dens = bind_rows(prior_dens, post_dens)

# plot
comb_dens %>%
  ggplot(aes())+
  facet_wrap(~name, nrow = 2)+
  stat_density(aes(value, linetype = type),
               geom = "line", size = 0.4, position="identity")+
  geom_rect(data = model_fit %>%
              filter(name %in% pars),
            aes(xmin = lower16, xmax = upper84, ymin = 0, ymax = Inf),
            inherit.aes = F, alpha = 0.2)+
  geom_vline(data = model_fit %>%
               filter(name %in% pars),
             aes(xintercept = middle),
             size = 0.8)+
  scale_y_continuous("Probability Density")+
  scale_x_continuous("Value",
                     limits=c(75,125))+
  scale_linetype_manual("",values=c(1,2))+
  theme_base+
  theme(legend.position = c(0.8,0.916))

# compare sd of prior and posterior
comb_dens %>%
  group_by(name, type) %>%
  summarize(sd = sd(value)) %>%
  spread(type, sd) %>%
  mutate(ratio = Posterior/Prior)





#==========
#========== Figure 5: Gammas
#==========

# set pars to plot
pars = c("gamma_1", "gamma_2")

# generate prior densities
set.seed(1)
prior_dens = priors %>%
  filter(name %in% paste0(pars, c("_prior"))) %>%
  split(.$name) %>%
  {lapply(1:length(.), function(y){
    data_frame(name = str_split(names(.)[y], "_prior") %>% 
                 map_chr(~as.character(.x[1])),
               value = rtruncnorm(n = 10000, a = 1, b = Inf,
                                  mean = .[[y]]$mean, sd = .[[y]]$sd),
               type = "Prior"
    )
  }
  )} %>%
  bind_rows()

# posterior densities
post_dens = params_full %>%
  select("chain", "step", pars) %>%
  gather(name, value, pars) %>%
  select(-chain, - step) %>%
  mutate(type = "Posterior")

# combine densities
comb_dens = bind_rows(prior_dens, post_dens)

# plot
comb_dens %>%
  ggplot(aes())+
  facet_wrap(~name, nrow = 2)+
  stat_density(aes(value, linetype = type),
               geom = "line", size = 0.4, position="identity")+
  geom_rect(data = model_fit %>%
              filter(name %in% pars),
            aes(xmin = lower16, xmax = upper84, ymin = 0, ymax = Inf),
            inherit.aes = F, alpha = 0.2)+
  geom_vline(data = model_fit %>%
               filter(name %in% pars),
             aes(xintercept = middle),
             size = 0.8)+
  scale_y_continuous("Probability Density")+
  scale_x_continuous("Value", 
                     limits=c(1,1.5))+
  scale_linetype_manual("",values=c(1,2))+
  theme_base+
  theme(legend.position = c(0.8,0.916))

# compare sd of prior and posterior
comb_dens %>%
  group_by(name, type) %>%
  summarize(sd = sd(value)) %>%
  spread(type, sd) %>%
  mutate(ratio = Posterior/Prior)





#==========
#========== Figure 6: sig_beta0 and sig_rho
#==========

# set pars to plot
pars = c("sig_beta0", "sig_rho")

# generate prior densities
set.seed(1)
prior_dens = priors %>%
  filter(name %in% paste0(pars, c("_prior"))) %>%
  split(.$name) %>%
  {lapply(1:length(.), function(y){
    data_frame(name = str_split(names(.)[y], "_prior") %>% 
                 map_chr(~as.character(.x[1])),
               value = rtruncnorm(n = 10000, a = 0, b = Inf,
                                  mean = .[[y]]$mean, sd = .[[y]]$sd),
               type = "Prior"
    )
  }
  )} %>%
  bind_rows()

# posterior densities
post_dens = params_full %>%
  select("chain", "step", pars) %>%
  gather(name, value, pars) %>%
  select(-chain, - step) %>%
  mutate(type = "Posterior")

# combine densities
comb_dens = bind_rows(prior_dens, post_dens)

# plot
comb_dens %>%
  ggplot(aes())+
  facet_wrap(~name, nrow = 2)+
  stat_density(aes(value, linetype = type),
               geom = "line", size = 0.4, position="identity")+
  geom_rect(data = model_fit %>%
              filter(name %in% pars),
            aes(xmin = lower16, xmax = upper84, ymin = 0, ymax = Inf),
            inherit.aes = F, alpha = 0.2)+
  geom_vline(data = model_fit %>%
               filter(name %in% pars),
             aes(xintercept = middle),
             size = 0.8)+
  scale_y_continuous("Probability Density")+
  scale_x_continuous("Value",
                     limits=c(0,0.3))+
  scale_linetype_manual("",values=c(1,2))+
  theme_base+
  theme(legend.position = c(0.8,0.916))

# compare sd of prior and posterior
comb_dens %>%
  group_by(name, type) %>%
  summarize(sd = sd(value)) %>%
  spread(type, sd) %>%
  mutate(ratio = Posterior/Prior)





#==========
#========== Figure 7: sig_alpha
#==========

# set pars to plot
pars = c("sig_alpha")

# generate prior densities
set.seed(1)
prior_dens = priors %>%
  filter(name %in% paste0(pars, c("_prior"))) %>%
  split(.$name) %>%
  {lapply(1:length(.), function(y){
    data_frame(name = str_split(names(.)[y], "_prior") %>% 
                 map_chr(~as.character(.x[1])),
               value = rtruncnorm(n = 10000, a = 0, b = Inf,
                                  mean = .[[y]]$mean, sd = .[[y]]$sd),
               type = "Prior"
    )
  }
  )} %>%
  bind_rows()

# posterior densities
post_dens = params_full %>%
  select("chain", "step", pars) %>%
  gather(name, value, pars) %>%
  select(-chain, - step) %>%
  mutate(type = "Posterior")

# combine densities
comb_dens = bind_rows(prior_dens, post_dens)

# plot
comb_dens %>%
  ggplot(aes())+
  facet_wrap(~name, nrow = 2)+
  stat_density(aes(value, linetype = type),
               geom = "line", size = 0.4, position="identity")+
  geom_rect(data = model_fit %>%
              filter(name %in% pars),
            aes(xmin = lower16, xmax = upper84, ymin = 0, ymax = Inf),
            inherit.aes = F, alpha = 0.2)+
  geom_vline(data = model_fit %>%
               filter(name %in% pars),
             aes(xintercept = middle),
             size = 0.8)+
  scale_y_continuous("Probability Density")+
  scale_x_continuous("Value",
                     limits=c(0,1))+
  scale_linetype_manual("",values=c(1,2))+
  theme_base+
  theme(legend.position = c(0.8,0.916))

# compare sd of prior and posterior
comb_dens %>%
  group_by(name, type) %>%
  summarize(sd = sd(value)) %>%
  spread(type, sd) %>%
  mutate(ratio = Posterior/Prior)
