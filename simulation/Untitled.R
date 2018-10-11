#==========
#========== Preliminaries
#==========

# load packages
library(tidyverse)
library(lubridate)
library(truncnorm)
library(nlme)

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
sonde_data = read_csv(paste0(sim_path,"input/beta0_alpha_rho_fixed/data_export.csv"))

# main analysis
priors = read_csv(paste0(sim_path,"input/beta0_alpha_rho_fixed/priors.csv"))
params_full = read_csv(paste0(sim_path,
                              "output/beta0_alpha_rho_fixed/fixed_pars_full.csv"))
model_fit = read_csv(paste0(sim_path,
                            "output/beta0_alpha_rho_fixed/summary_clean.csv"))





#==========
#========== Table I: Parameter estimates
#==========

# export parameter estimates
# model_fit %>%
#   filter(name %in% c("k0","k1","gamma_1","gamma_2","sig_beta0","sig_alpha","sig_rho","sig_proc")) %>%
#   select(name, middle, lower16, upper84) %>%
#   mutate(middle = round(middle, 3), 
#          lower16 = round(lower16, 3),
#          upper84 = round(upper84, 3)) %>%
#   write_csv("analyses/full_analysis/figures/table_i.csv")





#==========
#========== Fig 1: P-I Curve
#==========
set.seed(1)
test = model_fit %>%
  filter(name %in% c("beta0","rho","alpha")) %>%
  select(-lower16, -upper84) %>%
  mutate(middle = middle/1000) %>%
  spread(name, middle) %>%
  split(.$index) %>%
  lapply(function(x){y = data_frame(index = x$index, 
                                    beta0 = x$beta0, 
                                    rho = x$rho, 
                                    alpha = x$alpha,
                                    light = 1:200, 
                                    temp = 12)
  }) %>%
  bind_rows %>%
  mutate(nep = beta0*tanh(alpha*light/beta0) - rho) 
test = test %>% filter(index %in% sample(x= min(test$index):max(test$index), size = 20))

test2 = test %>%
  summarize(beta0 = median(beta0),
            rho = median(rho),
            alpha = median(alpha)) 

test3 = data_frame(beta0 = test2$beta0, rho = test2$rho, alpha = test2$alpha,
                   temp = 12, light = 1:200) %>%
  mutate(nep = beta0*tanh(alpha*light/beta0) - rho)

alph_dat = data_frame(light = -10:50, rho = test2$rho, alpha = test2$alpha, nep  = alpha*light - rho)
resp_dat = data_frame(light = c(-10, 150), nep = -test2$rho)
bet_dat = data_frame(light = 150, nep = c(-test2$rho, test2$beta0 - test2$rho - 0.01))
text_d = data_frame(light = c(5, 75, 175),
                    nep = c(0.05, -0.145, 0.05),
                    label = c("Initial Slope","ER","Max GPP"))

p = test  %>%
  ggplot(aes(light, nep))+
  geom_line(aes(group = index), alpha = 0.2, size = 0.2)+
  geom_line(data = test3, size = 0.8)+
  geom_line(data = alph_dat, linetype = 2)+
  geom_line(data = resp_dat, linetype = 2)+
  geom_line(data = bet_dat, linetype = 2)+
  geom_text(data = text_d, label = text_d$label)+
  scale_y_continuous(expression(NEP~
                                  "("*g~O[2]~m^{-2}~day^{-1}*")"))+
  scale_x_continuous(expression("Light ("*mu*mol~photons~m^{-2}~s^{-1}*")"))+
  theme_base

# examine & export
p
# ggsave("analyses/full_analysis/figures/fig_1.pdf", p, dpi = 300,
#        height = 4, width = 5, units = "in")




#==========
#========== Fig 2: Net flux
#==========

p = model_fit %>%
  filter(name %in% c("nep","air")) %>%
  select(name, day, index, middle) %>%
  spread(name, middle) %>%
  arrange(day, index) %>%
  select(-day, -index) %>%
  bind_cols(sonde_data %>% 
              na.omit() %>%
              arrange(year, yday, hour)) %>%
  select(year, yday, hour, par_int, do, nep, air) %>%
  arrange(year, yday, hour) %>%
  group_by(year, yday) %>%
  mutate(air = air/1000,
         nep = nep/1000,
         flux = c(diff(do), NA) - air) %>%
  group_by(year, yday) %>%
  summarize(flux = 24*mean(flux, na.rm=T),
            nep = 24*mean(nep, na.rm=T)) %>%
  gather(var, value, flux, nep) %>%
  mutate(var = factor(var, levels = c("flux","nep"), 
                      labels = c("Estimated NEP","Obesrved Flux"))) %>%
                      {ggplot(., aes(yday, value, color = var))+
                          facet_wrap(~year)+
                          geom_line(alpha = 0.8)+
                          scale_color_manual("",values = c("black","forestgreen"))+
                          scale_y_continuous(expression(Daily~Flux~
                                                          "("*g~O[2]~m^{-2}~day^{-1}*")"))+
                          scale_x_continuous("Day of Year", breaks = c(160, 190, 220))+
                          theme_base+
                          theme(legend.position = c(0.85, 0.65))}

# examine & export
p
# ggsave("analyses/full_analysis/figures/fig_2.pdf", p, dpi = 300,
#        height = 5, width = 5, units = "in")

# gray scale version
p = model_fit %>%
  filter(name %in% c("nep","air")) %>%
  select(name, day, index, middle) %>%
  spread(name, middle) %>%
  arrange(day, index) %>%
  select(-day, -index) %>%
  bind_cols(sonde_data %>% 
              na.omit() %>%
              arrange(year, yday, hour)) %>%
  select(year, yday, hour, par_int, do, nep, air) %>%
  arrange(year, yday, hour) %>%
  group_by(year, yday) %>%
  mutate(air = air/1000,
         nep = nep/1000,
         flux = c(diff(do), NA) - air) %>%
  group_by(year, yday) %>%
  summarize(flux = 24*mean(flux, na.rm=T),
            nep = 24*mean(nep, na.rm=T)) %>%
  gather(var, value, flux, nep) %>%
  mutate(var = factor(var, levels = c("flux","nep"), 
                      labels = c("Estimated NEP","Obesrved Flux"))) %>%
                      {ggplot(., aes(yday, value, linetype = var))+
                          facet_wrap(~year)+
                          geom_line()+
                          scale_linetype_manual("",values = c(1,2))+
                          scale_y_continuous(expression(Daily~Flux~
                                                          "("*g~O[2]~m^{-2}~day^{-1}*")"))+
                          scale_x_continuous("Day of Year", breaks = c(160, 190, 220))+
                          theme_base+
                          theme(legend.position = c(0.85, 0.65))}

# examine & export
p
# ggsave("analyses/full_analysis/figures/fig_2_gray.pdf", p, dpi = 300,
#        height = 5, width = 5, units = "in")









#==========
#========== Figure 3: beta0 and rho
#==========

# prepare data
p = model_fit %>%
  filter(name %in% c("beta0","rho")) %>%
  left_join(sonde_data %>%
              filter(is.na(unique_day)==F) %>%
              group_by(year, yday) %>%
              summarize(day = unique(unique_day))) %>%
  mutate(middle = middle/1000,
         lower16 = lower16/1000,
         upper84 = upper84/1000) %>%
         {ggplot(., aes(yday, middle, color=name))+
             facet_wrap(~year)+
             geom_ribbon(aes(ymin=lower16, ymax=upper84, fill = name),
                         linetype=0, alpha = 0.3)+
             geom_line(aes(color=name), size=0.5)+
             geom_text(data = data_frame(year = 2015, 
                                         yday = 190, 
                                         name = c("beta0","beta0","rho","rho"), 
                                         middle = c(0.49, 0.41, 0.19, 0.11)), 
                       aes(color = name), 
                       label=c("Max GPP",expression((beta[0])), "Basline ER", expression((rho))), 
                       hjust = 0)+
             scale_y_continuous(expression(Metabolism~Parameter~
                                             "("*g~O[2]~m^{-2}~day^{-1}*")"))+
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
              filter(is.na(unique_day)==F) %>%
              group_by(year, yday) %>%
              summarize(day = unique(unique_day))) %>%
  mutate(middle = middle/1000,
         lower16 = lower16/1000,
         upper84 = upper84/1000) %>%
         {ggplot(., aes(yday, middle, color=name))+
             facet_wrap(~year)+
             geom_hline(yintercept = mean(.$middle), alpha=0.5, size=0.2)+
             geom_ribbon(aes(ymin=lower16, ymax=upper84, group = name),
                         linetype=0, fill = "gray70")+
             geom_line(aes(color=name), size=0.5)+
             geom_text(data = data_frame(year = 2015, 
                                         yday = 190, 
                                         name = c("beta0","beta0","rho","rho"), 
                                         middle = c(0.49, 0.41, 0.19, 0.11)), 
                       color = "black", 
                       label=c("Max GPP",expression((beta[0])), "Basline ER", expression((rho))), 
                       hjust = 0)+
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
              filter(is.na(unique_day)==F) %>%
              group_by(year, yday) %>%
              summarize(day = unique(unique_day))) %>%
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
              filter(is.na(unique_day)==F) %>%
              group_by(year, yday) %>%
              summarize(day = unique(unique_day))) %>%
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
              filter(is.na(unique_day)==F) %>%
              group_by(year, yday) %>%
              summarize(day = unique(unique_day))) %>%
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
              filter(is.na(unique_day)==F) %>%
              group_by(year, yday) %>%
              summarize(day = unique(unique_day))) %>%
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
              filter(is.na(unique_day)==F) %>%
              group_by(year, yday) %>%
              summarize(day = unique(unique_day))) %>%
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





#==========
#========== Figure 8: beta0 and phycocyanin
#==========

# prepare data
beta0_phyc = model_fit %>%
  filter(name %in% c("beta0")) %>%
  left_join(sonde_data %>%
              filter(is.na(unique_day)==F) %>%
              filter(pcyv < 0.3) %>%
              group_by(year, yday) %>%
              summarize(day = unique(unique_day),
                        pcyv = mean(pcyv))) %>%
  select(year, yday, middle, pcyv) %>%
  rename(beta0 = middle) 

# plot
p = beta0_phyc %>%
  gather(var, val, beta0, pcyv)%>%
  group_by(var) %>%
  mutate(val = (val - mean(val, na.rm=T))/sd(val, na.rm=T)) %>%
  ungroup() %>%
  mutate(var = factor(var, levels=c("beta0","pcyv"), labels=c("Max GPP", "Phycocyanin"))) %>%
  ggplot(aes(yday, val, color = var))+
  facet_wrap(~year)+
  geom_line()+
  scale_color_manual("",values = c("black","cyan4"))+
  scale_y_continuous("Z-Score (dimensionless)", breaks = NULL)+
  scale_x_continuous("Day of Year", breaks = c(160, 190, 220))+
  theme_base+
  theme(legend.position = c(0.85, 0.9))

# examine & export
p
# ggsave("analyses/full_analysis/figures/fig_8.pdf", p, dpi = 300,
#        height = 5, width = 5, units = "in")

# gray scale version

# plot
p = beta0_phyc %>%
  group_by(var) %>%
  mutate(val = (val - mean(val, na.rm=T))/sd(val, na.rm=T)) %>%
  ungroup() %>%
  mutate(var = factor(var, levels=c("beta0","pcyv"), labels=c("Max GPP", "Phycocyanin"))) %>%
  ggplot(aes(yday, val, linetype = var))+
  facet_wrap(~year)+
  geom_line()+
  scale_linetype_manual("",values = c(1,2))+
  scale_y_continuous("Z-Score (dimensionless)", breaks = NULL)+
  scale_x_continuous("Day of Year", breaks = c(160, 190, 220))+
  theme_base+
  theme(legend.position = c(0.85, 0.9))

# examine & export
p
# ggsave("analyses/full_analysis/figures/fig_8_gray.pdf", p, dpi = 300,
#        height = 5, width = 5, units = "in")


m = gls(log(beta0) ~ pcyv, correlation = corCAR1(form = ~ yday|year), data = beta0_phyc)
summary(m)


#==========
#========== beta0, phycocyanin, and midges
#==========

# prepare beta0 and phycocyanin data
beta0_phyc = model_fit %>%
  filter(name %in% c("beta0")) %>%
  left_join(sonde_data %>%
              filter(is.na(unique_day)==F) %>%
              filter(pcyv < 0.3) %>%
              group_by(year, yday) %>%
              summarize(day = unique(unique_day),
                        pcyv = mean(pcyv))) %>%
  select(year, yday, middle, pcyv) %>%
  rename(beta0 = middle) 

# prepare midge data
midges_summary = midges %>%
  filter(sta %in% c(3, 33)) %>%
  mutate(year = year(sampledate),
         yday = yday(sampledate)) %>%
  group_by(year, yday, coreid) %>%
  summarize(tanyt = sum(tanyt/fract_count),
            chiro = sum(chiro/fract_count),
            midges = tanyt + chiro) %>%
  group_by(year, yday) %>%
  summarize(tanyt = mean(tanyt, na.rm=T),
            chiro = mean(chiro, na.rm=T),
            midges = mean(midges, na.rm=T))

# combine 
beta0_phyc_midge = beta0_phyc %>%
  left_join(midges_summary) %>%
  na.omit()

# plot
beta0_phyc_midge %>%
  gather(var, val, beta0, pcyv, midges) %>%
  group_by(var) %>%
  mutate(val = (val - mean(val, na.rm=T))/sd(val, na.rm=T)) %>%
  ungroup() %>%
  mutate(var = factor(var, levels=c("beta0","pcyv","midges"), labels=c("Max GPP", "Phycocyanin", "Midges"))) %>%
  ggplot(aes(yday, val, color = var))+
  facet_wrap(~year)+
  geom_line()+
  scale_color_manual("",values = c("black","cyan4","magenta4"))+
  scale_y_continuous("Z-Score (dimensionless)", breaks = NULL)+
  scale_x_continuous("Day of Year", breaks = c(160, 190, 220))+
  theme_base+
  theme(legend.position = c(0.14, 0.88))

m = gls(log(beta0) ~ pcyv + midges, correlation = corCAR1(form = ~ yday|year), data = beta0_phyc_midge)
summary(m)





#==========
#========== Fig 2: Sigmas
#==========

# plot
p = sim_params %>%
  select(type, sig_beta0, sig_alpha, sig_rho) %>%
  gather(par, value, sig_beta0, sig_alpha, sig_rho) %>%
  mutate(par = ifelse(par == "sig_beta0", "sigma[beta[0]]",
                      ifelse(par == "sig_alpha", "sigma[alpha]", "sigma[rho]"))) %>%
  select(par, value, type) %>%
  ggplot(aes(value, linetype = type))+
  facet_wrap(~par, labeller = label_parsed, nrow = 3, scales = "free")+
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
