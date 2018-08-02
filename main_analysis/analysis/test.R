# beta0
{model_fit %>% filter(name == "beta0")}$middle %>% summary


dy = {model_fit %>% filter(name == "beta0", middle  > 515 | middle < 302, day < 147)}$day

sonde_data %>%
  rename(day = D_M) %>%
  filter(day %in% dy) %>%
  mutate(do_flux = c(NA, 3.3*diff(do))) %>%
  select(day, par, do_flux) %>%
  bind_cols(model_fit %>% filter(name=="nep",day %in% dy) %>%
              rename(nep = middle) %>%
              select(nep)) %>%
  ggplot(aes(par, do_flux))+
  facet_wrap(~day)+
  geom_point(size=2)+
  geom_hline(yintercept = 0, size= 0.2)+
  geom_line(aes(y = nep), size=0.9)+
  scale_y_continuous(limits=c(-1000, 1500))+
  theme_base







# alpha
{model_fit %>% filter(name == "alpha")}$middle %>% summary


dy =
  {model_fit %>% filter(name == "alpha", middle  > 4.28 | middle < 1.31, day < 147)}$day

sonde_data %>%
  rename(day = D_M) %>%
  filter(day %in% dy) %>%
  mutate(do_flux = c(NA, 3.3*diff(do))) %>%
  select(day, par, do_flux) %>%
  bind_cols(model_fit %>% filter(name=="nep",day %in% dy) %>%
              rename(nep = middle) %>%
              select(nep)) %>%
  ggplot(aes(par, do_flux))+
  facet_wrap(~day)+
  geom_point(size=2)+
  geom_hline(yintercept = 0, size= 0.2)+
  geom_line(aes(y = nep), size=0.9)+
  scale_y_continuous(limits=c(-1000, 1500))+
  theme_base



#==========
#========== Figure 7: beta0 vs. rho
#==========

model_fit %>%
  filter(name %in% c("beta0","alpha")) %>%
  left_join(sonde_data %>%
              group_by(year, yday) %>%
              summarize(day = unique(D_M))) %>%
  select(-lower16, -upper84) %>%
  mutate(middle = middle/1000) %>%
  spread(name, middle) %>%
  {ggplot(.,aes(beta0, alpha))+
      geom_path(aes(color=factor(year)), size=0.4)+
      geom_point(aes(color=factor(year)), size=0.8, shape=1)+
      geom_point(data=. %>% group_by(year) %>% 
                   summarize(beta0 = beta0[1], alpha = alpha[1]),
                 aes(color=factor(year)),
                 size = 2, shape = 15)+
      geom_point(data=. %>% group_by(year) %>% 
                   summarize(beta0 = beta0[length(beta0)], alpha = alpha[length(alpha)]),
                 aes(color=factor(year)),
                 size = 2, shape = 17)+
      geom_text(data=. %>% group_by(year) %>% 
                  summarize(beta0 = beta0[1] - 0.02, alpha = alpha[1] + 0.003), 
                aes(label=year), size = 3.5)+
      scale_color_manual("",values=c("gray70","gray40","gray25","black"), guide=F)+
      theme_base}








