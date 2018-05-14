#==========
#========== Preliminaries
#==========

# load packages
library(tidyverse)

# import data and model fit
sonde_data = read_csv("data/sonde_prep.csv")
params_full = read_csv("model_output/o2_model/fixed_pars_full.csv")
post_pred = read_csv("model_output/o2_model/post_pred_full.csv")
model_fit = read_csv("model_output/o2_model/summary_clean.csv")

# base theme
theme_base = theme_bw()+
  theme(panel.grid=element_blank(),
        strip.background=element_blank(),
        text = element_text(size=12),
        strip.text = element_text(size=12),
        legend.text = element_text(size=12),
        axis.text=element_text(size=10, color="black"),
        axis.title.y=element_text(margin=margin(0,15,0,0)),
        axis.title.x=element_text(margin=margin(15,0,0,0)))





#==========
#========== Figure 1: Posterior Predictive Check
#==========

# split by type (simualted vs. observe)
post_pred_split = post_pred %>%
  gather(var, chi_squared) %>%
  mutate(stoch = strsplit(var,"\\_") %>% map_chr(~.x[2]),
         type = strsplit(var,"\\_") %>% map_chr(~.x[3])) %>%
  split(.$type) 

# combine with separete columns for simulated and observed
# plot
lapply(1:length(post_pred_split), function(x){
  y = post_pred_split[[x]] %>% select(chi_squared, stoch)
  names(y) = c(names(post_pred_split)[x], "stoch")
  return(y)}) %>% 
  bind_cols %>% 
  select(-stoch1) %>%
  mutate(Stoch = ifelse(stoch=="obs","Observation Error","Process Error"),
         Stoch = factor(Stoch, levels=c("Process Error","Observation Error"))) %>%
  ggplot(aes(real,sim))+
  facet_wrap(~Stoch, nrow=2)+
  geom_point(size=2, alpha=0.5)+
  geom_abline(intercept=0, slope=1)+
  scale_y_continuous(expression(chi^2~"Simulated"), limits=c(4750,5600), 
                     breaks=c(4900,5200,5500))+
  scale_x_continuous(expression(chi^2~"Real"),limits=c(4750,5600), 
                     breaks=c(4900,5200,5500))+
  coord_equal()+
  theme_base




#==========
#========== Figure 2: Net Daily Fluxes
#==========

model_fit %>%
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
         value = value/1000) %>%
  ggplot(aes(yday, value, color=var, size=var))+
  facet_wrap(~year)+
  geom_hline(yintercept = 0, alpha=0.5, size=0.5)+
  geom_line()+
  scale_color_manual("",values=c("gray50","black"))+
  scale_size_manual("",values=c(0.5,0.7))+
  scale_y_continuous(expression("Net Flux (g "*O[2]~m^{-2}~day^{-1}*")"),
                     limits=c(-2,2), breaks=c(-1.5,0,1.5))+
  scale_x_continuous("Day of Year", 
                     limits=c(150,240), breaks=c(160,180,200,220))+
  theme_base+
  theme(legend.position = c(0.87,0.924))




#==========
#========== Figure 3: Daily Metabolism
#==========

model_fit %>%
  filter(name %in% c("GPP","ER","NEP")) %>%
  left_join(sonde_data %>%
              group_by(year, yday) %>%
              summarize(day = unique(D_M))) %>%
  full_join(sonde_data %>% expand(year,yday,name=c("GPP","ER","NEP"))) %>%
  mutate(middle = ifelse(name=="ER",-middle,middle)/1000,
         lower16 = ifelse(name=="ER",-lower16,lower16)/1000,
         upper84 = ifelse(name=="ER",-upper84,upper84)/1000,
         name = factor(name, levels=c("GPP","NEP","ER"))) %>%
  {ggplot(., aes(yday, middle, color=name))+
      facet_wrap(~year)+
      geom_hline(yintercept = 0, alpha=0.5, size=0.5)+
      geom_ribbon(aes(ymin=lower16, ymax=upper84, fill=name),
                  linetype=0, alpha=0.35)+
      geom_line(size=0.6)+
      geom_line(data=.%>% filter(is.na(middle)==F), size=0.3, linetype=2)+
      scale_color_manual("",values=c("dodgerblue","gray40","firebrick"))+
      scale_fill_manual("",values=c("dodgerblue","gray40","firebrick"))+
      scale_y_continuous(expression("Metabolism (g "*O[2]~m^{-2}~day^{-1}*")"))+
      scale_x_continuous("Day of Year")+
      theme_base+
      theme(legend.position = c(0.9,0.9))
  }
  




#==========
#========== Figure 4: RESP vs. GPP
#==========

# prep data
met_d = model_fit %>%
  filter(name %in% c("GPP","ER","NEP")) %>%
  left_join(sonde_data %>%
              group_by(year, yday) %>%
              summarize(day = unique(D_M))) %>%
  select(-lower16, -upper84) %>%
  mutate(middle = middle/1000) %>%
  spread(name, middle) 

# plot
met_d %>%
  ggplot(aes(GPP, ER, color=factor(year)))+
  geom_point(size=2)+
  geom_abline(intercept=0, slope=1)+
  scale_color_manual("",values=c("gray60","gray40","gray20","black"))+
  scale_y_continuous(expression("ER (g "*O[2]~m^{-2}~day^{-1}*")"), limits=c(2,9))+
  scale_x_continuous(expression("GPP (g "*O[2]~m^{-2}~day^{-1}*")"), limits=c(2,9))+
  coord_equal()+
  theme_base+
  theme(legend.position = c(0.1,0.85))

# calculate correlation
met_d %>% 
with(cor(GPP, ER))

# calculate correlation by year
met_d %>%
  split(.$year) %>%
  lapply(function(x){with(x,cor(GPP, ER))})

# means
model_fit %>%
  filter(name %in% c("GPP_mean","ER_mean","NEP_mean")) %>%
  mutate(middle = middle/1000,
         lower16 = lower16/1000,
         upper84 = upper84/1000)





#==========
#========== Figure 5: Daily Responses
#==========

# prep data
resp_d = model_fit %>%
  filter(name %in% c("beta0","rho")) %>%
  left_join(sonde_data %>%
              group_by(year, yday) %>%
              summarize(day = unique(D_M))) %>%
  mutate(middle = middle/1000,
         lower16 = lower16/1000,
         upper84 = upper84/1000) %>%
  arrange(year,yday)

# plots
resp_d %>%
  ggplot(aes(yday, middle, color=name))+
  facet_wrap(~year, labeller=label_parsed)+
  geom_hline(yintercept = 0.3, alpha=0.5, size=0.5)+
  geom_ribbon(aes(ymin=lower16, ymax=upper84, fill=name),
                linetype=0, alpha=0.35)+
  geom_line(size=0.6)+
  scale_color_manual("",values=c("dodgerblue","firebrick"), 
                     labels = c(expression(beta^0), expression(rho)))+
  scale_fill_manual("",values=c("dodgerblue","firebrick"), 
                    labels = c(expression(beta^0), expression(rho)))+
  scale_y_continuous(expression("Metabolism Parameter (g "*O[2]~m^{-2}~h^{-1}*")"),
                     breaks=c(0.15,0.3,0.45))+
  scale_x_continuous("Day of Year")+
  theme_base+
  theme(legend.text.align = 0,
        legend.position = c(0.9,0.9))

# calculate correlation
resp_d %>% 
  select(-lower16, -upper84) %>%
  spread(name, middle) %>%
with(cor(β, ρ))

# calculate correlation (without 2017)
resp_d %>% 
  select(-lower16, -upper84) %>%
  spread(name, middle) %>%
  filter(year!=2017) %>%
  with(cor(β, ρ))

# calculate correlation for each year
resp_d %>% 
  select(-lower16, -upper84) %>%
  spread(name, middle) %>%
  split(.$year) %>%
  lapply(function(x){with(x,cor(β, ρ))})
  





#==========
#========== Figure 6: beta0 vs. rho
#==========

model_fit %>%
  filter(name %in% c("beta0","rho")) %>%
  left_join(sonde_data %>%
              group_by(year, yday) %>%
              summarize(day = unique(D_M))) %>%
  select(-lower16, -upper84) %>%
  mutate(middle = middle/1000) %>%
  spread(name, middle) %>%
  {ggplot(.,aes(beta0, rho))+
      geom_path(aes(color=factor(year)), size=0.8, alpha=0.8)+
      geom_point(data=. %>% group_by(year) %>% 
                   summarize(beta0 = beta0[1], rho = rho[1]),
                 aes(color=factor(year)),
                 size = 3)+
      geom_point(data=. %>% group_by(year) %>% 
                   summarize(beta0 = beta0[length(beta0)], rho = rho[length(rho)]),
                 aes(color=factor(year)),
                 size = 4, shape=17)+
      geom_text(data=. %>% group_by(year) %>% 
                  summarize(beta0 = beta0[1] - 0.015, rho = rho[1] + 0.003),
                aes(label=year))+
      scale_color_manual("",values=c("gray70","gray40","gray10","black"), guide=F)+
      scale_x_continuous(expression(beta^0~"("*g~O[2]~m^{-2}~day^{-1}*")"))+
      scale_y_continuous(expression(rho~"("*g~O[2]~m^{-2}~day^{-1}*")"))+
      theme_base}





#==========
#========== Appendix: Fixed Parameters
#==========

# pairs plot
ggpairs(params_full %>% select(-step, -chain, -lp__))

# values
model_fit %>%
  filter(name %in% 
           c("alpha","gamma_1","gamma_2","sig_beta0","sig_rho","sig_proc")) %>%
  select(-index, -day)
  

