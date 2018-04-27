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
        strip.text = element_text(size=10),
        axis.text=element_text(size=10),
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
  mutate(Stoch = ifelse(stoch=="obs","Observation Error","Process Error")) %>%
  ggplot(aes(real,sim))+
  facet_wrap(~Stoch, nrow=2)+
  geom_point(size=2, alpha=0.5)+
  geom_abline(intercept=0, slope=1)+
  scale_y_continuous("Chi-squared Simulated",limits=c(4750,5600), 
                     breaks=c(4900,5200,5500))+
  scale_x_continuous("Chi-squared Real",limits=c(4750,5600), 
                     breaks=c(4900,5200,5500))+
  coord_equal()+
  theme_base




#==========
#========== Figure 2: Total Daily Fluxes
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
  scale_y_continuous(expression("Net Flux (g "*O[2]~m^{-2}~day^{-1}*")"))+
  scale_x_continuous("Day of Year")+
  theme_base




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
      theme_base
  }
  



#==========
#========== Figure 4: Daily Responses
#==========

model_fit %>%
  filter(name %in% c("beta0","rho")) %>%
  left_join(sonde_data %>%
              group_by(year, yday) %>%
              summarize(day = unique(D_M))) %>%
  mutate(middle = middle/1000,
         lower16 = lower16/1000,
         upper84 = upper84/1000) %>%
  arrange(year,yday) %>%
  ggplot(aes(yday, middle, color=name))+
  facet_wrap(~year)+
  geom_ribbon(aes(ymin=lower16, ymax=upper84, fill=name),
                linetype=0, alpha=0.35)+
  geom_line(size=0.6)+
  scale_color_manual("",values=c("dodgerblue","firebrick"))+
  scale_fill_manual("",values=c("dodgerblue","firebrick"))+
  scale_y_continuous(expression("Metabolism Parameter (g "*O[2]~m^{-2}~day^{-1}*")"))+
  scale_x_continuous("Day of Year")+
  theme_base







