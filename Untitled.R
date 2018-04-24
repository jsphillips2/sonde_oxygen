# fit_clean = read_csv("model_output/summaries_clean.csv")

fit_clean %>%
  filter(name %in% c("beta0")) %>%
  left_join(sonde_prep3 %>%
              group_by(year, yday) %>%
              summarize(day = unique(D_M))) %>%
  ggplot(aes(yday, middle, color=name))+
  facet_wrap(~year)+
  geom_line(size=0.8)+
  geom_ribbon(aes(ymin=lower25, ymax=upper75, fill=name), size=0, alpha=0.5)+
  scale_color_manual(values=c("dodgerblue"))+
  scale_fill_manual(values=c("dodgerblue"))+
  theme_bw()

fit_clean %>%
  filter(name %in% c("rho")) %>%
  left_join(sonde_prep3 %>%
              group_by(year, yday) %>%
              summarize(day = unique(D_M))) %>%
  ggplot(aes(yday, middle, color=name))+
  facet_wrap(~year)+
  geom_line(size=0.8)+
  geom_ribbon(aes(ymin=lower25, ymax=upper75, fill=name), size=0, alpha=0.5)+
  scale_color_manual(values=c("firebrick"))+
  scale_fill_manual(values=c("firebrick"))+
  theme_bw()

library(GGally)

fixed_par_v = c("alpha","gamma_1","gamma_2","sig_beta0","sig_rho","sig_proc","lp__")
fixed_pars = rstan::extract(fit, pars=fixed_par_v) %>%
  lapply(as_data_frame) %>%
  bind_cols()
names(fixed_pars) = fixed_par_v
 
ggpairs(fixed_pars)
