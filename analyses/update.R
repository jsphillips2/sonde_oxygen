# prepare data
flux_d <- lapply(list("main","main_revised"), 
                 function(x){read_csv(paste0("model/output/",x,"/summary_clean.csv")) %>%
                     filter(name %in% c("NEP","AIR")) %>%
                     left_join(read_csv(paste0("model/input/",x,"/sonde_prep.csv")) %>%
                                 filter(is.na(unique_day)==F) %>%
                                 group_by(year, yday) %>%
                                 summarize(day = unique(unique_day))) %>%
                     full_join(read_csv(paste0("model/input/",x,"/sonde_prep.csv")) %>%
                                 expand(year,yday,name=c("NEP","AIR"))) %>%
                     mutate(middle = middle/1000,
                            name = factor(name, levels=c("NEP","AIR"))) %>%
                     mutate(analysis = x)}) %>% 
  bind_rows() %>%
  mutate(analysis = ifelse(analysis=="main_revised","Revsied","Original"),
         analysis = factor(analysis, levels = c("Revsied","Original")),
         name = ifelse(name=="AIR","EXC","NEP"))

# prepare data
flux_d %>%
  ggplot(aes(yday, middle, color = name, linetype = analysis))+
  facet_wrap(~year)+
  geom_line(size = 0.5)+
  geom_hline(yintercept = 0, size = 0.3, alpha = 0.5)+
  scale_color_manual("",values=c("firebrick","dodgerblue"))+
  scale_linetype_discrete("")+
  scale_y_continuous(expression(Daily~DO~Flux~"("*g~O[2]~m^{-2}~d^{-1}*")"))+
  scale_x_continuous("Day of Year", breaks = c(170, 200, 230), limits = c(150, 250))+
  theme(legend.position = "top")







# prepare data
beta_rho_d = read_csv("model/output/main/summary_clean.csv") %>%
  filter(name %in% c("beta0","rho")) %>%
  left_join(read_csv("model/input/main/sonde_prep.csv") %>%
              filter(is.na(unique_day)==F) %>%
              group_by(year, yday) %>%
              summarize(day = unique(unique_day))) %>%
  mutate(type = "Original") %>%
  bind_rows(read_csv("model/output/main_revised/summary_clean.csv") %>%
              filter(name %in% c("beta0","rho")) %>%
              left_join(read_csv("model/input/main_revised/sonde_prep.csv") %>%
                          filter(is.na(unique_day)==F) %>%
                          group_by(year, yday) %>%
                          summarize(day = unique(unique_day))) %>%
              mutate(type = "Revised"))

# plot
beta_rho_d %>%
  ggplot(aes(yday, middle, color=name))+
  facet_wrap(~year)+
  geom_line(aes(linetype = type, size = type), alpha = 0.7)+
  geom_text(data = tibble(year = 2015, yday = 187, 
                          name = c("beta0","rho"), 
                          middle = c(490, 240)), 
            aes(color = name), 
            label=c("Max GPP", "Basline ER"), 
            hjust = 0)+
  geom_text(data = tibble(year = 2015, yday = 187, 
                          name = c("beta0","rho"), 
                          middle = c(410,160)), 
            aes(color = name), 
            label=c("(beta^0)", "(rho)"), parse = T, 
            hjust = 0)+
  scale_y_continuous(expression(Metabolism~Parameter~"("*mg~O[2]~m^{-2}~h^{-1}*")"))+
  scale_color_manual("",values=c("dodgerblue","firebrick2"), guide = F)+
  scale_linetype_manual("",values=c(2,1))+
  scale_size_manual(values=c(0.7, 0.4))+
  scale_x_continuous("Day of Year", breaks = c(170, 200, 230), limits = c(150, 250))+
  guides(size = F)+
  theme(legend.position = "top")






# prepare data
alpha_d = read_csv("model/output/main/summary_clean.csv") %>%
  filter(name %in% c("alpha")) %>%
  left_join(read_csv("model/input/main/sonde_prep.csv") %>%
              filter(is.na(unique_day)==F) %>%
              group_by(year, yday) %>%
              summarize(day = unique(unique_day))) %>%
  mutate(type = "Original") %>%
  bind_rows(read_csv("model/output/main_revised/summary_clean.csv") %>%
              filter(name %in% c("alpha")) %>%
              left_join(read_csv("model/input/main_revised/sonde_prep.csv") %>%
                          filter(is.na(unique_day)==F) %>%
                          group_by(year, yday) %>%
                          summarize(day = unique(unique_day))) %>%
              mutate(type = "Revised"))

# plot
alpha_d %>%
  filter(name %in% c("alpha")) %>%
  ggplot(aes(yday, middle))+
  facet_wrap(~year)+
  geom_line(aes(linetype = type, size = type), alpha = 0.7)+
  scale_linetype_manual("",values=c(2,1))+
  scale_size_manual(values=c(0.7, 0.4))+
  scale_y_continuous(expression(Initial~Slope~(alpha)~"("*mg~O[2]~s~mu*mol-photons^{-1}~h^{-1}*")"))+
  scale_x_continuous("Day of Year", breaks = c(160, 190, 220))+
  guides(size = F)+
  theme(legend.position = "top")
