#==========
#========== Preliminaries
#==========

# load packages
library(tidyverse)
library(truncnorm)

# import data and model fit
input_dir = "analyses/test_analysis/model_fit/"
sonde_data = read_csv(paste0(input_dir,"input/sonde_prep.csv"))
priors = read_csv(paste0(input_dir,"input/priors.csv"))
params_full = read_csv(paste0(input_dir,"output/fixed_pars_full.csv"))
post_pred = read_csv(paste0(input_dir,"output/post_pred_full.csv"))
model_fit = read_csv(paste0(input_dir,"output/summary_clean.csv"))

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
#========== Figure 1: Posterior predictive check
#==========





#==========
#========== Figure 2: Net daily fluxes
#==========

f2 = model_fit %>%
  filter(name=="Flux") %>%
  left_join(sonde_data %>%
              group_by(year, yday) %>%
              arrange(year,yday,hour) %>%
              mutate(d_do = c(NA,diff(do))) %>%
              summarize(day = unique(D_M),
                        d_do = 24*mean(d_do,na.rm=T))) %>%
  full_join(sonde_data %>% expand(year,yday)) %>%
  select(year,yday,middle,d_do) %>%
  gather(var, value, c(middle,d_do)) %>%
  mutate(var = ifelse(var=="middle","NEP + AIR","Observed"),
         var = factor(var, levels=c("Observed","NEP + AIR")),
         value = 3.3*value/1000) %>%
  ggplot(aes(yday, value, color=var, size=var))+
  facet_wrap(~year)+
  geom_hline(yintercept = 0, alpha=0.5, size=0.5)+
  geom_line()+
  scale_color_manual("",values=c("gray50","black"))+
  scale_size_manual("",values=c(0.5,0.7))+
  scale_y_continuous(expression(Net~Flux~"("*g~O[2]~m^{-2}~day^{-1}*")"),
                     limits=c(-6,8), breaks=c(-4,0,4))+
  scale_x_continuous("Day of Year", 
                     limits=c(147,205), breaks=c(157,176,195))+
  theme_base+
  theme(legend.position = c(0.65,0.916))

# examine and export
f2
# ggsave("main_analysis/analysis/figures/fig_2.pdf", f2, dpi = 300,
#        height = 5, width = 5, units = "in")





#==========
#========== Figure 3: Gas Exchange
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

# plot
bind_rows(prior_dens, post_dens) %>%
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

# plot
bind_rows(prior_dens, post_dens) %>%
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

# plot
bind_rows(prior_dens, post_dens) %>%
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

# plot
bind_rows(prior_dens, post_dens) %>%
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

# plot
bind_rows(prior_dens, post_dens) %>%
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


