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

sonde_data %>%
  group_by(year, yday) %>%
  summarize(par_int = mean(par0, na.rm=T)*exp(-mean(ext*1.5, na.rm=T))) %>%
  ggplot(aes(yday, par_int))+
  facet_wrap(~year)+
  geom_line()+
  theme_base

sonde_data %>%
  group_by(year, yday) %>%
  summarize(par_int = mean(par_int, na.rm=T),
            pcyv = mean(pcyv, na.rm=T)) %>%
  filter(pcyv < 0.3) %>%
  gather(var, value, par_int, pcyv) %>%
  group_by(var) %>%
  mutate(value = (value - mean(value, na.rm=T))/sd(value, na.rm=T)) %>%
  ggplot(aes(yday, value, color=var))+
  facet_wrap(~year)+
  geom_line()+
  theme_base


sonde_data %>%
  filter(pcyv < 0.3) %>%
  group_by(year, yday) %>%
  summarize(pcyv = mean(pcyv, na.rm=T))%>%
  ggplot(aes(yday, pcyv))+
  facet_wrap(~year)+
  geom_line()+
  theme_base

model_fit %>%
  filter(name %in% c("beta0")) %>%
  left_join(sonde_data %>%
              filter(is.na(unique_day)==F) %>%
              filter(pcyv < 0.3) %>%
              group_by(year, yday) %>%
              summarize(day = unique(unique_day),
                        pcyv = mean(pcyv))) %>%
  select(year, yday, middle, pcyv) %>%
  rename(beta0 = middle) %>%
  gather(var, val, beta0, pcyv) %>%
  group_by(var) %>%
  mutate(val = (val - mean(val, na.rm=T))/sd(val, na.rm=T)) %>%
  ggplot(aes(yday, val, color = var))+
  facet_wrap(~year)+
  geom_line()+
  theme_base



model_fit %>%
  filter(name %in% c("GPP")) %>%
  left_join(sonde_data %>%
              filter(is.na(unique_day)==F) %>%
              filter(pcyv < 0.3) %>%
              group_by(year, yday) %>%
              summarize(day = unique(unique_day),
                        pcyv = mean(pcyv))) %>%
  select(year, yday, middle, pcyv) %>%
  rename(GPP = middle) %>%
  gather(var, val, GPP, pcyv) %>%
  group_by(var) %>%
  mutate(val = (val - mean(val, na.rm=T))/sd(val, na.rm=T)) %>%
  ggplot(aes(yday, val, color = var))+
  facet_wrap(~year)+
  geom_line()+
  theme_base

model_fit %>%
  filter(name %in% c("alpha")) %>%
  left_join(sonde_data %>%
              filter(is.na(unique_day)==F) %>%
              group_by(year, yday) %>%
              summarize(day = unique(unique_day),
                        par_int = mean(par_int))) %>%
  select(year, yday, middle, par_int) %>%
  rename(alpha = middle) %>%
  gather(var, val, alpha, par_int) %>%
  group_by(var) %>%
  mutate(val = (val - mean(val, na.rm=T))/sd(val, na.rm=T)) %>%
  ggplot(aes(yday, val, color = var))+
  facet_wrap(~year)+
  geom_line()+
  theme_base

