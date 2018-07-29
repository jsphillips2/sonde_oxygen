
{model_fit %>% filter(name == "rho")}$middle %>% summary


dy = {model_fit %>% filter(name == "rho", middle < 85)}$day

# dy = {resp_d %>% filter(name == "beta0", middle < 0.28)}$day

sonde_data %>%
  rename(day = D_M) %>%
  filter(day %in% dy) %>%
  mutate(do_flux = c(NA, 3.3*diff(do))) %>%
  select(day, par, do_flux) %>%
  bind_cols(model_fit %>% filter(name=="nep",day %in% dy) %>%
              rename(nep = middle) %>%
              select(nep)) %>%
  ggplot(aes(par, do_flux))+
  facet_wrap(~day, scales="free")+
  geom_point(size=2)+
  geom_line(aes(y = nep), size=0.9)+
  theme_base


params_full %>%
  select(-lp__, -step) %>%
  gather(variable, value, -chain) %>%
  ggplot(aes(value, color=chain))+
  facet_wrap(~variable, scales="free")+
  geom_density()+
  theme_base





