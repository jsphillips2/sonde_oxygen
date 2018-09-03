#==========
#========== Preliminaries
#==========

# load packages
library(tidyverse)

# import data and model fit
import_file = "fixed_beta_rho"
sonde_data = read_csv(paste0("simulation/simulated_data/",import_file,"/data_export.csv"))
params_full = read_csv(paste0(output_path,import_file,"/post_pred_full.csv"))
post_pred = read_csv(paste0(output_path,import_file,"/daily_full.csv"))
model_fit = read_csv(paste0(output_path,import_file,"/summary_clean.csv"))

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
#========== Figure 5: beta and rho
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
  ggplot(aes(day, middle, color=name))+
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
  # scale_x_continuous("Day of Year", 
  #                    limits=c(150,240), breaks=c(165,195,225))+
  theme_base+
  theme(legend.text.align = 0,
        legend.position = c(0.9,0.9))

