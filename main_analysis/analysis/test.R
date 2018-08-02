
{model_fit %>% filter(name == "beta0")}$middle %>% summary


dy = {model_fit %>% filter(name == "beta0", middle  > 500, day < 147)}$day

dy = c(183:191)
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
  facet_wrap(~day)+
  geom_point(size=2)+
  geom_hline(yintercept = 0, size= 0.2)+
  geom_line(aes(y = nep), size=0.9)+
  scale_y_continuous(limits=c(-1500, 1500))+
  theme_base


sonde_data %>%
  rename(day = D_M) %>%
  filter(day %in% dy) %>%
  mutate(do_flux = c(3.3*diff(do), NA)) %>%
  select(day, par, do_flux) %>%
  bind_cols(model_fit %>% filter(name=="nep",day %in% dy) %>%
              rename(nep = middle) %>%
              select(nep)) %>%
  bind_cols(model_fit %>% filter(name=="air",day %in% dy) %>%
              rename(air = middle) %>%
              select(air)) %>%
  mutate(flux = do_flux - air) %>%
  ggplot(aes(par, flux))+
  facet_wrap(~day, scales="free")+
  geom_point(size=2)+
  geom_line(aes(y = nep), size=0.9)+
  theme_base


