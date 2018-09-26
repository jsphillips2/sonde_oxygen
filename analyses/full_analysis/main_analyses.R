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

# set paths for main and simulation analyses
main_path = "analyses/full_analysis/model_fit/"
sim_path = "analyses/full_analysis/simulation/"

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
types = c("_fixed","beta0_alpha_rho_fixed")
sim_params = lapply(types, function(x){
  read_csv(paste0(sim_path, "output/", x, "/fixed_pars_full.csv")) %>%
    cbind(read_csv(paste0(sim_path, "input/", x, "/type_data.csv")) %>%
            select(var, value) %>%
            spread(var, value)) %>%
    tbl_df %>%
    mutate(type = x)
}) %>%
  bind_rows




#==========
#========== Table I: Parameter estimates
#==========

# export parameter estimates
# model_fit %>%
#   filter(name %in% c("k0","k1","gamma_1","gamma_2","sig_beta0","sig_alpha","sig_rho","sig_prco")) %>%
#   select(name, middle, lower16, upper84) %>%
#   mutate(middle = round(middle, 3), 
#          lower16 = round(lower16, 3),
#          upper84 = round(upper84, 3)) %>%
#   write_csv("analyses/full_analysis/figures/table_i.csv")






#==========
#========== Fig 1: Posterior Predictive Check
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
post_pred_split = post_pred_fn(post_pred10, 10) 

# plot
p = post_pred_split %>%
  ggplot(aes(real,sim))+
  facet_wrap(~Stoch, nrow=2)+
  geom_point(size=1.5, alpha=0.5)+
  geom_abline(intercept=0, slope=1)+
  scale_y_continuous(expression(chi^2~"Simulated"), limits=c(5800,6800), 
                     breaks = c(6000, 6300, 6600))+
  scale_x_continuous(expression(chi^2~"Real"), limits=c(5800,6800), 
                     breaks = c(6000, 6300, 6600))+
  coord_equal()+
  theme_base+
  theme(legend.position = c(0.15,0.9))

# examine & export
p
# ggsave("analyses/full_analysis/figures/fig_1.pdf", p, dpi = 300,
#        height = 5, width = 3, units = "in")

# examine estiamtes for sig_proc
list("sig_obs10" =  model_fit,"sig_obs100" =   model_fit100) %>%
  lapply(function(x){
    x %>% 
      filter(name == "sig_proc") %>%
      select(-index, -day)})

# bayesian "p-value"
post_pred_split %>%
  bind_rows(post_pred_fn(post_pred100, 100)) %>%
  mutate(p = ifelse(sim > real, 1, 0)) %>%
  group_by(Stoch, sig_obs) %>%
  summarize(p = mean(p))





#==========
#========== Fig 2: Sigmas
#==========

# plot
p = pars_fit %>%
  select(sig_beta0, sig_alpha, sig_rho, fix_beta0, fix_alpha, fix_rho) %>%
  gather(par, value, sig_beta0, sig_alpha, sig_rho) %>%
  mutate(fixed = ifelse(par == "sig_beta0", fix_beta0, 
                        ifelse(par == "sig_alpha", fix_alpha, fix_rho)),
         fixed = ifelse(fixed == T, "Fixed", "Not Fixed"),
         par = ifelse(par == "sig_beta0", "sigma[beta[0]]", 
                      ifelse(par == "sig_alpha", "sigma[alpha]", "sigma[rho]"))) %>%
  select(par, value, fixed) %>%
  ggplot(aes(value, linetype = fixed))+
  facet_wrap(~par, labeller = label_parsed, nrow = 3, scales = "free_y")+
  stat_density(position = "identity", geom = "line")+
  scale_y_continuous("Posterior Probability Density", breaks = NULL)+
  scale_x_continuous("Value (Dimensionsless)", breaks = c(0.05, 0.13, 0.21))+
  scale_linetype_manual("", values = c(2,1))+
  theme_base+
  theme(legend.position = c(0.75,0.92))

# examine & export
p
# ggsave("analyses/full_analysis/figures/fig_2.pdf", p, dpi = 300,
#        height = 6, width = 3, units = "in")






#==========
#========== Figure 3: beta0 and rho
#==========

# prepare data
p = model_fit %>%
  filter(name %in% c("beta0","rho")) %>%
  left_join(sonde_data %>%
              group_by(year, yday) %>%
              summarize(day = unique(D_M))) %>%
  mutate(middle = middle/1000,
         lower16 = lower16/1000,
         upper84 = upper84/1000) %>%
  {ggplot(., aes(yday, middle, color=name))+
      facet_wrap(~year)+
      geom_hline(yintercept = mean(.$middle), alpha=0.5, size=0.2)+
      geom_ribbon(aes(ymin=lower16, ymax=upper84, fill = name),
              linetype=0, alpha = 0.3)+
      geom_line(aes(color=name), size=0.5)+
      geom_text(data = data_frame(year = 2015, yday = 190, name = c("beta0","rho"), middle = c(0.45, 0.15)), 
                aes(color = name), label=c(expression(beta[0]), expression(rho)), hjust = 0)+
      scale_y_continuous(expression(Metabolism~Parameter~
                                  "("*g~O[2]~m^{-2}~day^{-1}*")"),
                         breaks = c(0.1, 0.3, 0.5))+
      scale_color_manual("",values=c("dodgerblue","firebrick2"), guide = F)+
      scale_fill_manual("",values=c("dodgerblue","firebrick2"), guide = F)+
      scale_x_continuous("Day of Year", breaks = c(160, 190, 220))+
      theme_base
    }

# examine & export
p
# ggsave("analyses/full_analysis/figures/fig_3.pdf", p, dpi = 300,
#        height = 5, width = 5, units = "in")

# gray scale version
p = model_fit %>%
  filter(name %in% c("beta0","rho")) %>%
  left_join(sonde_data %>%
              group_by(year, yday) %>%
              summarize(day = unique(D_M))) %>%
  mutate(middle = middle/1000,
         lower16 = lower16/1000,
         upper84 = upper84/1000) %>%
         {ggplot(., aes(yday, middle, color=name))+
             facet_wrap(~year)+
             geom_hline(yintercept = mean(.$middle), alpha=0.5, size=0.2)+
             geom_ribbon(aes(ymin=lower16, ymax=upper84, group = name),
                         linetype=0, fill = "gray70")+
             geom_line(aes(color=name), size=0.5)+
             geom_text(data = data_frame(year = 2015, yday = 190, name = c("beta0","rho"), middle = c(0.45, 0.15)), 
                       aes(group = name), color = "black",
                       label=c(expression(beta[0]), expression(rho)), hjust = 0)+
             scale_linetype_manual("",
                                   values=c(1,5,2),
                                   guide = F)+
             scale_y_continuous(expression(Metabolism~Parameter~
                                             "("*g~O[2]~m^{-2}~day^{-1}*")"),
                                breaks = c(0.1, 0.3, 0.5))+
             scale_color_manual(values=c("black","white"), guide = F)+
             scale_x_continuous("Day of Year", breaks = c(160, 190, 220))+
             theme_base
         }

# examine & export
p
# ggsave("analyses/full_analysis/figures/fig_3_gray.pdf", p, dpi = 300,
#        height = 5, width = 5, units = "in")





#==========
#========== Figure 4: alpha
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
  {ggplot(., aes(yday, middle))+
      facet_wrap(~year)+
      geom_hline(yintercept = mean(.$middle), alpha=0.5, size=0.2)+
      geom_ribbon(aes(ymin=lower16, ymax=upper84),
              linetype=0, alpha=0.35)+
      geom_line(size=0.6)+
      scale_y_continuous(expression(alpha~"("*mg~O[2]~s~mu*mol-photons^{-1}~h^{-1}*")"))+
      scale_x_continuous("Day of Year", breaks = c(160, 190, 220))+
      theme_base
    }

# examine & export
p
# ggsave("analyses/full_analysis/figures/fig_4.pdf", p, dpi = 300,
#        height = 5, width = 5, units = "in")





#==========
#========== Figure 5: beta0 vs rho
#==========

p = model_fit %>%
  filter(name %in% c("beta0","rho")) %>%
  left_join(sonde_data %>%
              group_by(year, yday) %>%
              summarize(day = unique(D_M))) %>%
  select(-lower16, -upper84) %>%
  mutate(middle = middle/1000) %>%
  spread(name, middle) %>%
  {ggplot(.,aes(beta0, rho))+
      facet_wrap(~year)+
      geom_hline(yintercept = mean(.$rho), size = 0.3, linetype = 2)+
      geom_vline(xintercept = mean(.$beta0), size = 0.2, linetype = 2)+
      geom_path(aes(group=factor(year)), size = 0.6, alpha = 0.7)+
      geom_point(data=. %>% group_by(year) %>% 
                   summarize(beta0 = beta0[1], rho = rho[1]),
                 size = 3, shape = 15)+
      geom_point(data=. %>% group_by(year) %>% 
                   summarize(beta0 = beta0[length(beta0)], rho = rho[length(rho)]),
                 size = 3, shape = 17)+
      scale_y_continuous(expression(rho~"("*g~O[2]~m^{-2}~h^{-1}*")"),
                         breaks = c(0.12, 0.16, 0.20))+
      scale_x_continuous(expression(beta^0~"("*g~O[2]~m^{-2}~h^{-1}*")"),
                         breaks = c(0.31, 0.39, 0.47))+
      theme_base}

# examine & export
p
# ggsave("analyses/full_analysis/figures/fig_5.pdf", p, dpi = 300,
#        height = 5, width = 5, units = "in")


# correlation between beta0 and rho
model_fit %>%
  filter(name %in% c("beta0","rho")) %>%
  select(name, middle, day) %>%
  spread(name, middle) %>%
  {cor(.$beta0, .$rho)}





#==========
#========== Fig 6: Daily Metabolism
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
         {ggplot(., aes(yday, middle, color = name))+
             facet_wrap(~year)+
             geom_hline(yintercept = 0, alpha=0.5, size=0.2)+
             geom_line(data = . %>% filter(is.na(middle) == F), aes(group = name), size = 0.5, linetype = 3)+
             geom_ribbon(aes(ymin=lower16, ymax=upper84, fill = name),
                         linetype=0, alpha = 0.3)+
             geom_line(aes(color=name), size=0.5)+
             geom_text(data = data_frame(year = 2015, yday = 190, name = c("GPP","NEP","ER"), middle = c(5,2,-3)), 
                       aes(label = name), hjust = 0)+
             scale_color_manual("",values=c("dodgerblue","black","firebrick2"), guide = F)+
             scale_fill_manual("",values=c("dodgerblue","black","firebrick2"), guide = F)+
             scale_y_continuous(expression(Metabolism~
                                             "("*g~O[2]~m^{-2}~day^{-1}*")"), breaks = c(-6, 0, 6))+
             scale_x_continuous("Day of Year", breaks = c(160, 190, 220))+
             theme_base
         }


# examine & export
p
# ggsave("analyses/full_analysis/figures/fig_6.pdf", p, dpi = 300,
#        height = 5, width = 5, units = "in")

# gray scale version
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
             facet_wrap(~year)+
             geom_hline(yintercept = 0, alpha=0.5, size=0.2)+
             geom_line(data = . %>% filter(is.na(middle) == F), aes(group = name), size = 0.5, linetype = 3)+
             geom_ribbon(aes(ymin=lower16, ymax=upper84, group = name),
                         linetype=0, fill = "gray70")+
             geom_line(aes(color=name), size=0.5)+
             geom_text(data = data_frame(year = 2015, yday = 190, name = c("GPP","NEP","ER"), middle = c(5,2,-3)), 
                       aes(label = name), hjust = 0)+
             scale_linetype_manual("",
                                   values=c(1,5,2),
                                   guide = F)+
             scale_y_continuous(expression(Metabolism~
                                             "("*g~O[2]~m^{-2}~day^{-1}*")"), breaks = c(-6, 0, 6))+
             scale_color_manual(values=c("black","gray40","white"), guide = F)+
             scale_x_continuous("Day of Year", breaks = c(160, 190, 220))+
             theme_base
         }

# examine & export
p
# ggsave("analyses/full_analysis/figures/fig_6_gray.pdf", p, dpi = 300,
#        height = 5, width = 5, units = "in")

# metabolism estimates
model_fit %>% 
  filter(name %in% c("GPP_mean","ER_mean","NEP_mean")) %>%
  mutate(lower16 = lower16/1000,
         middle = middle/1000,
         upper84 = upper84/1000) %>%
  select(name, lower16, middle, upper84)





#==========
#========== Figure 7: GPP vs ER
#==========

p = model_fit %>%
  filter(name %in% c("GPP","ER","NEP")) %>%
  left_join(sonde_data %>%
              group_by(year, yday) %>%
              summarize(day = unique(D_M))) %>%
  select(-lower16, -upper84) %>%
  mutate(middle = middle/1000) %>%
  spread(name, middle) %>%
  ggplot(aes(GPP, ER))+
  geom_point(size=2, alpha = 0.5)+
  geom_abline(intercept=0, slope=1)+
  scale_y_continuous(expression(ER~"("*g~O[2]~m^{-2}~day^{-1}*")"), 
                     limits=c(2,9), breaks=c(3,5.5,8))+
  scale_x_continuous(expression(GPP~"("*g~O[2]~m^{-2}~day^{-1}*")"), 
                     limits=c(2,9), breaks=c(3,5.5,8))+
  coord_equal()+
  theme_base+
  theme(legend.position = c(0.25,0.9), legend.direction = "horizontal")

# examine & export
p
# ggsave("analyses/full_analysis/figures/fig_7.pdf", p, dpi = 300,
#        height = 4, width = 3, units = "in")


# correlation between GPP and ER
model_fit %>%
  filter(name %in% c("GPP","ER")) %>%
  select(name, middle, day) %>%
  spread(name, middle) %>%
  {cor(.$GPP, .$ER)}

