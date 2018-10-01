test = model_fit %>%
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
  mutate(flux = c(diff(do), NA) - air/1000,
         nep = nep/1000)

test %>%
  group_by(year, yday) %>%
  summarize(flux = 24*mean(flux, na.rm=T),
            nep = 24*mean(nep, na.rm=T)) %>%
  ggplot(aes(yday, flux))+
  facet_wrap(~year)+
  geom_line(color="blue")+
  geom_line(aes(yday, nep),color="red")+
  theme_base


test %>%
  group_by(year, yday) %>%
  summarize(flux = 24*mean(flux, na.rm=T),
            nep = 24*mean(nep, na.rm=T)) %>%
  gather(var, value, flux, nep) %>%
  ggplot(aes(yday, value, color = var))+
  facet_wrap(~year)+
  geom_line()+
  scale_color_manual(values = c("black","forestgreen"))+
  theme_base


test %>%
  group_by(year, yday) %>%
  summarize(flux = 24*mean(flux, na.rm=T),
            nep = 24*mean(nep, na.rm=T)) %>%
  ggplot(aes(nep, flux))+
  facet_wrap(~year)+
  geom_point()+
  geom_abline(intercept = 0, slope = 1)+
  theme_base

test %>%
  group_by(year, yday) %>%
  summarize(flux = 24*mean(flux, na.rm=T),
            nep = 24*mean(nep, na.rm=T)) %>%
  {cor.test(.$flux, .$nep)}
  
  
test %>%
  group_by(year, yday) %>%
  summarize(flux = 24*mean(flux, na.rm=T),
            nep = 24*mean(nep, na.rm=T)) %>%
  gather(var, value, flux, nep) %>%
  ggplot(aes(yday, value, color = var))+
  facet_wrap(~year)+
  geom_point(alpha = 0.5, size = 1)+
  geom_smooth(method="gam", formula = y~s(x, k = 10), se = F)+
  scale_color_manual(values=c("dodgerblue","firebrick"))+
  theme_base





