#==========
#========== Preliminaries
#==========

# load packages
library(tidyverse)
library(GGally)
library(truncnorm)

# import data and model fit
input_dir = "analyses/test_analysis/model_fit/"
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
#========== beta and rho
#==========

# prep data
beta_rho = model_fit %>%
  filter(name %in% c("beta0","rho")) %>%
  left_join(sonde_data %>%
              group_by(year, yday) %>%
              summarize(day = unique(D_M))) %>%
  mutate(middle = middle/1000,
         lower16 = lower16/1000,
         upper84 = upper84/1000) %>%
  arrange(year,yday)

# plot
beta_rho %>%
  ggplot(aes(yday, middle, color=name))+
  facet_wrap(~year, labeller=label_parsed)+
  geom_hline(yintercept = 0.3, alpha=0.5, size=0.2)+
  geom_ribbon(aes(ymin=lower16, ymax=upper84, fill=name),
              linetype=0, alpha=0.35)+
  geom_line(size=0.6)+
  scale_color_manual("",values=c("dodgerblue","firebrick"), 
                     labels = c(expression(beta^0), expression(rho)))+
  scale_fill_manual("",values=c("dodgerblue","firebrick"), 
                    labels = c(expression(beta^0), expression(rho)))+
  scale_y_continuous(expression("Metabolism Parameter (g "*O[2]~m^{-2}~h^{-1}*")"),
                     breaks=c(0.15,0.3,0.45))+
  scale_x_continuous("Day of Year", 
                     limits=c(151,201), breaks=c(160,176,192))+
  theme_base+
  theme(legend.text.align = 0,
        legend.position = c(0.2,0.85))





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
  geom_ribbon(aes(ymin=lower16, ymax=upper84),
              linetype=0, alpha=0.35)+
  geom_line(size=0.6)+
  scale_y_continuous(expression(alpha~"("*mg~O[2]~s~mu*mol-photons^{-1}~h^{-1}*")"))+
  scale_x_continuous("Day of Year", 
                     limits=c(151,201), breaks=c(160,176,192))+
  theme_base







