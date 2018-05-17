model_fit %>%
  filter(name=="Flux") %>%
  left_join(sonde_data %>%
              group_by(year, yday) %>%
              arrange(year,yday,hour) %>%
              mutate(d_do = c(NA,diff(do))) %>%
              summarize(day = unique(D_M),
                        d_do = 24*mean(d_do,na.rm=T))) %>%
  full_join(sonde_data %>% expand(year,yday)) %>%
  select(year,yday,middle,d_do) %>%
  gather(var, value, c(middle,d_do)) %>%
  mutate(var = ifelse(var=="middle","NEP + AIR","Observed"),
         var = factor(var, levels=c("Observed","NEP + AIR")),
         value = value/1000) %>%
  # filter(year == 2017, yday %in% c(210:217)) %>%
  ggplot(aes(yday, value, color=var, size=var))+
  facet_wrap(~year)+
  geom_hline(yintercept = 0, alpha=0.5, size=0.5)+
  geom_line()+
  scale_color_manual("",values=c("gray50","black"))+
  scale_size_manual("",values=c(0.5,0.7))+
  scale_y_continuous(expression(Net~Flux~"("*g~O[2]~m^{-2}~day^{-1}*")"),
                     limits=c(-2,2), breaks=c(-1.5,0,1.5))+
  theme_base+
  theme(legend.position = c(0.87,0.91))

model_fit %>%
  filter(name %in% c("nep","air")) %>%
  group_by(day, name) %>%
  mutate(hour = 1:24) %>%
  select(-lower16,-upper84) %>%
  spread(name, middle) %>%
  mutate(flux = nep + air) %>%
left_join(sonde_data %>%
            group_by(year, yday) %>%
            arrange(year,yday,hour) %>%
            mutate(d_do = c(NA,diff(do)),
                   day = D_M)) %>%
  select(year,yday,day,hour,par, nep, air, flux, d_do) %>%
  filter(year==2017, yday %in% c(210:217)) %>%
  mutate(air = air/3300,
         nep = nep/3300,
         flux = flux/3300,
         d_do = d_do/1000) %>%
  ggplot(aes(par, d_do))+
  facet_wrap(~yday, scales = "free")+
  geom_point(alpha = 0.5)+
  geom_line(aes(par, flux), size=0.8)+
  theme_base  

sonde_data %>%
  select(year, par, temp, wspeed) %>%
  na.omit() %>%
  group_by(year) %>%
  summarize(par_mean = mean(par),
            par_se = sd(par)/sqrt(length(par)),
            temp_mean = mean(temp),
            temp_se = sd(temp)/sqrt(length(temp)))
    
    

test = model_fit %>%
  filter(name=="Flux") %>%
  left_join(sonde_data %>%
              group_by(year, yday) %>%
              arrange(year,yday,hour) %>%
              mutate(d_do = c(NA,diff(do))) %>%
              summarize(day = unique(D_M),
                        d_do = 24*mean(d_do,na.rm=T))) %>%
  full_join(sonde_data %>%
              select(year,yday,par,temp,wspeed) %>%
              group_by(year, yday) %>%
              summarize_all(mean)) %>%
  select(year,yday,par,temp,wspeed,middle,d_do) %>%
  gather(var, value, c(middle,d_do)) %>%
  mutate(var = ifelse(var=="middle","Flux","Observed"),
         var = factor(var, levels=c("Observed","Flux")),
         value = value/1000) %>%
  spread(var,value) %>%
  mutate(Error = Observed - Flux) 

test %>%
  ggplot(aes(wspeed, Error, color=factor(year)))+
  geom_point()+
  geom_hline(yintercept=0, size=0.5)+
  theme_base

test %>%
  ggplot(aes(par, Error, color=factor(year)))+
  geom_point()+
  geom_hline(yintercept=0, size=0.5)+
  theme_base

test %>%
  ggplot(aes(temp, Error, color=factor(year)))+
  geom_point()+
  geom_hline(yintercept=0, size=0.5)+
  theme_base

test %>%
  mutate(wspeed = (wspeed - mean(wspeed,na.rm = T))/(3*sd(wspeed,na.rm = T))) %>%
  ggplot(aes(yday, Observed))+
  facet_wrap(~year)+
  geom_line(color="gray50",size=0.5)+
  geom_line(aes(yday, Flux), color="black",size=0.7)+
  geom_line(aes(yday, wspeed), color="red")+
  geom_hline(yintercept=0, size=0.5)+
  theme_base

test %>%
  mutate(wspeed = (wspeed - mean(wspeed,na.rm = T))/(3*sd(wspeed,na.rm = T))) %>%
  ggplot(aes(yday, -Error))+
  facet_wrap(~year)+
  geom_line()+
  geom_line(aes(yday, wspeed), color="red")+
  geom_hline(yintercept=0, size=0.5)+
  theme_base
