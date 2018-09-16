#==========
#========== Preliminaries
#==========

# load packages
library(tidyverse)
library(GGally)
library(truncnorm)

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
#========== Import data
#==========

# set pahts for main and simulation analyses
main_path = "analyses/full_analysis/model_fit/"

# posterior predictive checik for both sig_obs = 10 and sig_obs = 100
# parameters for sig_obs = 100
post_pred10 = read_csv(paste0(main_path,
                              "output/sig_obs10/post_pred_full.csv"))
post_pred100 = read_csv(paste0(main_path,
                               "output/sig_obs100/post_pred_full.csv"))
model_fit100 = read_csv(paste0(main_path,
                               "output/sig_obs100/summary_clean.csv"))





#==========
#========== Posterior Predictive Check
#==========

# create function to process data
post_pred_fn = function(data, sig_obs){
  data %>% 
    gather(var, chi_squared) %>%
    mutate(stoch = strsplit(var,"\\_") %>% map_chr(~.x[2]),
           type = strsplit(var,"\\_") %>% map_chr(~.x[3])) %>%
    split(.$type) %>%
    {lapply(1:length(.), function(x){
      y = .[[x]] %>% select(chi_squared, stoch)
      names(y) = c(names(.)[x], "stoch")
      return(y)}) %>% 
        bind_cols %>% 
        select(-stoch1) %>%
        mutate(Stoch = ifelse(stoch=="obs","Observation Error","Process Error"),
               Stoch = factor(Stoch, levels=c("Process Error","Observation Error")),
               sig_obs = sig_obs)} 
}

# process data for sig_obs = 10 and sig_obs = 100
post_pred_split = post_pred_fn(post_pred10, 10) %>%
  bind_rows(post_pred_fn(post_pred100, 100))

# plot
post_pred_split %>%
  ggplot(aes(real,sim, color = factor(sig_obs)))+
  facet_wrap(~Stoch, nrow=2)+
  geom_point(size=1.5, alpha=0.5)+
  geom_abline(intercept=0, slope=1)+
  # scale_y_continuous(expression(chi^2~"Simulated"), limits=c(900,1500))+
  # scale_x_continuous(expression(chi^2~"Real"), limits=c(650,1500))+
  scale_color_manual(expression(sigma[obs]), values = c("forestgreen","black"))+
  # coord_equal()+
  theme_base+
  theme(legend.position = c(0.15,0.9))

# examine & export\