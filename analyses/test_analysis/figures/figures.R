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
main_path = "analyses/test_analysis/model_fit/"
sim_path = "analyses/test_analysis/simulation/"

# data
sonde_data = read_csv(paste0(main_path,"input/sonde_prep.csv"))

# main analysis
priors = read_csv(paste0(main_path,"input/priors.csv"))
params_full = read_csv(paste0(main_path,
                              "output/sig_obs10/fixed_pars_full.csv"))
model_fit = read_csv(paste0(main_path,
                            "output/sig_obs10/summary_clean.csv"))


# posterior predictive checik for both sig_obs = 10 and sig_obs = 100
# parameters for sig_obs = 100
post_pred10 = read_csv(paste0(main_path,
                            "output/sig_obs10/post_pred_full.csv"))
post_pred100 = read_csv(paste0(main_path,
                            "output/sig_obs100/post_pred_full.csv"))
model_fit100 = read_csv(paste0(main_path,
                                 "output/sig_obs100/summary_clean.csv"))


# simulation anlaysis
types = list.files(paste0(sim_path,"output/")) 
sim_params = lapply(types, function(x){
  read_csv(paste0(sim_path, "output/", x, "/fixed_pars_full.csv")) %>%
    cbind(read_csv(paste0(sim_path, "input/", x, "/type_data.csv")) %>%
            select(var, value) %>%
            spread(var, value)) %>%
    tbl_df %>%
    mutate(type = x)
}) %>%
  bind_rows

sig_sim_params = sim_params %>%
  select(sig_beta0, sig_alpha, sig_rho, fix_beta0, fix_alpha, fix_rho) 



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
p = post_pred_split %>%
  ggplot(aes(real,sim, color = factor(sig_obs)))+
  facet_wrap(~Stoch, nrow=2)+
  geom_point(size=1.5, alpha=0.5)+
  geom_abline(intercept=0, slope=1)+
  scale_y_continuous(expression(chi^2~"Simulated"), limits=c(900,1500))+
  scale_x_continuous(expression(chi^2~"Real"), limits=c(650,1500))+
  scale_color_manual(expression(sigma[obs]), values = c("forestgreen","black"))+
  coord_equal()+
  theme_base+
  theme(legend.position = c(0.15,0.9))

# examine & export
p
# ggsave("analyses/test_analysis/figures/fig_1.pdf", p, dpi = 300,
#        height = 5, width = 4, units = "in")

# examine estiamtes for sig_proc
list("sig_obs10" =  model_fit,"sig_obs100" =   model_fit100) %>%
  lapply(function(x){
    x %>% 
      filter(name == "sig_proc") %>%
      select(-index, -day)})
  




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
p = comb_dens %>%
  mutate(name = ifelse(name == "k0", 'k[0]', "k[1]")) %>%
  ggplot(aes())+
  facet_wrap(~name, nrow = 2, scales = "free", labeller = label_parsed)+
  stat_density(aes(value, linetype = type),
               geom = "line", size = 0.4, position="identity")+
  geom_rect(data = model_fit %>%
              filter(name %in% pars) %>%
              mutate(name = ifelse(name == "k0", 'k[0]', "k[1]")),
            aes(xmin = lower16, xmax = upper84, ymin = 0, ymax = Inf),
            inherit.aes = F, alpha = 0.2)+
  geom_vline(data = model_fit %>%
               filter(name %in% pars) %>%
               mutate(name = ifelse(name == "k0", 'k[0]', "k[1]")),
             aes(xintercept = middle),
             size = 0.8)+
  scale_y_continuous("Posterior Probability", breaks = NULL)+
  scale_x_continuous("Value")+
  scale_linetype_manual("",values=c(1,2))+
  theme_base+
  theme(legend.position = c(0.8,0.916))

# examine & export
p
# ggsave("analyses/test_analysis/figures/fig_2.pdf", p, dpi = 300,
#        height = 5, width = 4, units = "in")

# compare sd of prior and posterior
comb_dens %>%
  group_by(name, type) %>%
  summarize(sd = sd(value)) %>%
  spread(type, sd) %>%
  mutate(ratio = Posterior/Prior)





#==========
#========== Daily Metabolism
#==========

# plot
p = model_fit %>%
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
             geom_hline(yintercept = 0, alpha=0.5, size=0.2)+
             geom_ribbon(aes(ymin=lower16, ymax=upper84, group = name),
                         linetype=0, alpha=0.35)+
             geom_line(aes(linetype=name), size=0.6)+
             scale_linetype_manual("",
                                   values=c(1,3,2),
                                   guide = guide_legend(
                                     direction = "horizontal",
                                     label.position = "top"))+
             scale_y_continuous(expression(Metabolism~
                                             "("*g~O[2]~m^{-2}~day^{-1}*")"),
                                limits=c(-10,10), breaks=c(-6,0,6))+
             scale_x_continuous("Day of Year", 
                                limits=c(151,201), breaks=c(160,176,192))+
             theme_base+
             theme(legend.position = c(0.5,0.93),
                   legend.key.width = unit(3, "line"),
                   legend.key.height = unit(0.8, "line"),
                   legend.text.align = 0.5)
         }

# examine & export
p
# ggsave("analyses/test_analysis/figures/fig_3.pdf", p, dpi = 300,
#        height = 5, width = 4, units = "in")





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
p = met_p %>%
  ggplot(aes(yday, middle, linetype=name))+
  geom_ribbon(aes(ymin=lower16, ymax=upper84, group=name),
              linetype=0, alpha=0.35)+
  geom_line(size=0.6)+
  scale_linetype_manual("",values=c(1,2),
                        guide = guide_legend(
                          direction = "horizontal",
                          label.position = "top"),
                        labels = c(expression(beta^0), expression(rho)))+
  scale_y_continuous(expression("Metabolism Parameter (g "*O[2]~m^{-2}~h^{-1}*")"),
                     limits=c(0.1,0.5), breaks=c(0.15,0.3,0.45))+
  scale_x_continuous("Day of Year", 
                     limits=c(151,201), breaks=c(160,176,192))+
  theme_base+
  theme(legend.position = c(0.5,0.94),
        legend.key.width = unit(3, "line"),
        legend.key.height = unit(0.8, "line"),
        legend.text.align = 0.5)

# examine & export
p
# ggsave("analyses/test_analysis/figures/fig_4.pdf", p, dpi = 300,
#        height = 5, width = 4, units = "in")





#==========
#========== alpha
#==========

p = model_fit %>%
  filter(name %in% c("alpha")) %>%
  left_join(sonde_data %>%
              group_by(year, yday) %>%
              summarize(day = unique(D_M))) %>%
  mutate(middle = middle,
         lower16 = lower16,
         upper84 = upper84) %>%
  arrange(year,yday) %>%
  ggplot(aes(yday, middle))+
  geom_ribbon(aes(ymin=lower16, ymax=upper84),
              linetype=0, alpha=0.35)+
  geom_line(size=0.6)+
  scale_y_continuous(expression(alpha~"("*mg~O[2]~s~mu*mol-photons^{-1}~h^{-1}*")"),
                     limits=c(1,5))+
  scale_x_continuous("Day of Year", 
                     limits=c(151,201), breaks=c(160,176,192))+
  theme_base

# examine & export
p
# ggsave("analyses/test_analysis/figures/fig_5.pdf", p, dpi = 300,
#        height = 5, width = 4, units = "in")





#==========
#========== sig_beta0
#==========

p = sig_sim_params %>%
  mutate(fix_beta0 = ifelse(fix_beta0 == T, "Fixed", "Not Fixed"),
         fix_alpha = ifelse(fix_alpha == T, "alpha:~Fixed", "alpha:~Not~Fixed"),
         fix_rho = ifelse(fix_rho == T, "rho:~Fixed", "rho:~Not~Fixed")) %>%
  ggplot(aes(sig_beta0, linetype = fix_beta0))+
  facet_grid(fix_rho ~ fix_alpha, labeller = label_parsed)+
  stat_density(geom = "line", position = "identity")+
  scale_linetype_manual(expression(beta[0]), values = c(2,1))+
  scale_y_continuous("Posterior Probability", breaks=NULL)+
  scale_x_continuous(expression(sigma[beta[0]]), breaks = c(0.03, 0.09, 0.15))+
  theme_base+
  theme(legend.position = c(0.85,0.85))

# examine & export
p
# ggsave("analyses/test_analysis/figures/fig_6.pdf", p, dpi = 300,
#        height = 5, width = 5, units = "in")





#==========
#========== sig_rho
#==========

p = sig_sim_params %>%
  mutate(fix_beta0 = ifelse(fix_beta0 == T, 
                            "beta[0]:~Fixed", "beta[0]:~Not~Fixed"),
         fix_alpha = ifelse(fix_alpha == T, "alpha:~Fixed", "alpha:~Not~Fixed"),
         fix_rho = ifelse(fix_rho == T, "Fixed", "Not Fixed")) %>%
  ggplot(aes(sig_rho, linetype = fix_rho))+
  facet_grid(fix_beta0 ~ fix_alpha, labeller = label_parsed)+
  stat_density(geom = "line", position = "identity")+
  scale_linetype_manual(expression(rho), values = c(2,1))+
  scale_y_continuous("Posterior Probability", breaks=NULL)+
  scale_x_continuous(expression(sigma[rho]), breaks = c(0.03, 0.10, 0.17))+
  theme_base+
  theme(legend.position = c(0.85,0.85))


# examine & export
p
# ggsave("analyses/test_analysis/figures/fig_7.pdf", p, dpi = 300,
#        height = 5, width = 5, units = "in")





#==========
#========== sig_alpha
#==========

p = sig_sim_params %>%
  mutate(fix_beta0 = ifelse(fix_beta0 == T, 
                            "beta[0]:~Fixed", "beta[0]:~Not~Fixed"),
         fix_alpha = ifelse(fix_alpha == T, "Fixed", "Not Fixed"),
         fix_rho = ifelse(fix_rho == T, "rho:~Fixed", "rho:~Not~Fixed")) %>%
  ggplot(aes(sig_alpha, linetype = fix_alpha))+
  facet_grid(fix_rho ~ fix_beta0, labeller = label_parsed)+
  stat_density(geom = "line", position = "identity")+
  scale_linetype_manual(expression(alpha), values = c(2,1))+
  scale_y_continuous("Posterior Probability", breaks=NULL)+
  scale_x_continuous(expression(sigma[alpha]), breaks = c(0.07, 0.22, 0.37))+
  theme_base+
  theme(legend.position = c(0.85,0.85))

# examine & export
p
# ggsave("analyses/test_analysis/figures/fig_8.pdf", p, dpi = 300,
#        height = 5, width = 5, units = "in")


