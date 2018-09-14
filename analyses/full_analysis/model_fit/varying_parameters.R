#==========
#========== Preliminaries
#==========

# load packages
library(tidyverse)
library(GGally)
library(truncnorm)

# import data and model fit
input_dir = "analyses/full_analysis/model_fit/"
sonde_data = read_csv(paste0(input_dir,"input/sonde_prep.csv"))
priors = read_csv(paste0(input_dir,"input/priors.csv"))
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
#========== Daily Metabolism
#==========

# plot
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
         {ggplot(., aes(yday, middle))+
             facet_wrap(~year)+
             geom_hline(yintercept = 0, alpha=0.5, size=0.2)+
             geom_ribbon(data = . %>% filter(is.na(middle) == F), aes(ymin=lower16, ymax=upper84, group = name),
                         linetype=0, alpha=0.35)+
             geom_line(aes(linetype=name), size=0.6)+
             scale_linetype_manual("",
                                   values=c(1,3,2),
                                   guide = guide_legend(direction = "horizontal",
                                     label.position = "top"))+
             scale_y_continuous(expression(Metabolism~
                                             "("*g~O[2]~m^{-2}~day^{-1}*")"))+
             # scale_x_continuous("Day of Year", 
             #                    limits=c(151,201), breaks=c(160,176,192))+
             theme_base+
             theme(legend.position = "top",
                   legend.key.width = unit(3, "line"),
                   legend.key.height = unit(0.8, "line"),
                   legend.text.align = 0.5)
         }





#==========
#========== beta0 and rho
#==========

# prepare data
met_p = model_fit %>%
  filter(name %in% c("beta0","rho")) %>%
  left_join(sonde_data %>%
              group_by(year, yday) %>%
              summarize(day = unique(D_M))) %>%
  mutate(middle = middle/1000,
         lower16 = lower16/1000,
         upper84 = upper84/1000) %>%
  arrange(year,yday)

# plot
met_p %>%
  ggplot(aes(yday, middle, linetype=name))+
  facet_wrap(~year)+
  geom_ribbon(aes(ymin=lower16, ymax=upper84, group=name),
              linetype=0, alpha=0.35)+
  geom_line(size=0.6)+
  geom_hline(yintercept = 0.3, size = 0.3, alpha = 0.5)+
  scale_linetype_manual("",values=c(1,2),
                        guide = guide_legend(direction = "horizontal",
                          label.position = "top"),
                        labels = c(expression(beta^0), expression(rho)))+
  scale_y_continuous(expression("Metabolism Parameter (g "*O[2]~m^{-2}~h^{-1}*")"))+
  # scale_x_continuous("Day of Year", 
  #                    limits=c(151,201), breaks=c(160,176,192))+
  theme_base+
  theme(legend.position = "top",
        legend.key.width = unit(3, "line"),
        legend.key.height = unit(0.8, "line"),
        legend.text.align = 0.5)





#==========
#========== alpha
#==========

model_fit %>%
  filter(name %in% c("alpha")) %>%
  left_join(sonde_data %>%
              group_by(year, yday) %>%
              summarize(day = unique(D_M))) %>%
  mutate(middle = middle,
         lower16 = lower16,
         upper84 = upper84) %>%
  arrange(year,yday) %>%
  ggplot(aes(yday, middle))+
  facet_wrap(~year)+
  geom_hline(yintercept = 3, size = 0.3, alpha = 0.5)+
  geom_ribbon(aes(ymin=lower16, ymax=upper84),
              linetype=0, alpha=0.35)+
  geom_line(size=0.6)+
  scale_y_continuous(expression(alpha~"("*mg~O[2]~s~mu*mol-photons^{-1}~h^{-1}*")"))+
  scale_x_continuous("Day of Year")+
  theme_base





#==========
#========== GPP vs ER
#==========

model_fit %>%
  filter(name %in% c("GPP","ER","NEP")) %>%
  left_join(sonde_data %>%
              group_by(year, yday) %>%
              summarize(day = unique(D_M))) %>%
  select(-lower16, -upper84) %>%
  mutate(middle = middle/1000) %>%
  spread(name, middle) %>%
  ggplot(aes(GPP, ER))+
  geom_point(size=2)+
  geom_abline(intercept=0, slope=1)+
  scale_y_continuous(expression(ER~"("*g~O[2]~m^{-2}~day^{-1}*")"), 
                     limits=c(2,9), breaks=c(4,6,8))+
  scale_x_continuous(expression(GPP~"("*g~O[2]~m^{-2}~day^{-1}*")"), 
                     limits=c(2,9), breaks=c(4,6,8))+
  coord_equal()+
  theme_base+
  theme(legend.position = c(0.25,0.9), legend.direction = "horizontal")





#==========
#========== beta0 vs rho
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
      geom_path(aes(group=factor(year)), size=0.5, alpha = 0.7)+
      geom_point(data=. %>% group_by(year) %>% 
                   summarize(beta0 = beta0[1], rho = rho[1]),
                 size = 2.5, shape = 15)+
      geom_point(data=. %>% group_by(year) %>% 
                   summarize(beta0 = beta0[length(beta0)], rho = rho[length(rho)]),
                 size = 2.5, shape = 17)+
      geom_text(data=. %>% group_by(year) %>% 
                  summarize(beta0 = beta0[1] - 0.02, rho = rho[1] + 0.003), 
                aes(label=year), size = 3.5)+
      scale_y_continuous(expression(rho~"("*g~O[2]~m^{-2}~h^{-1}*")"))+
      scale_x_continuous(expression(beta^0~"("*g~O[2]~m^{-2}~h^{-1}*")"))+
      theme_base}

