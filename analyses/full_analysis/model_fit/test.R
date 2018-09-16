model_fit %>%
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
             geom_line(aes(linetype=name), size=0.6)+
             scale_linetype_manual("",
                                   values=c(1,5,2),
                                   guide = guide_legend(direction = "horizontal",
                                                        label.position = "top"))+
             scale_y_continuous(expression(Metabolism~
                                             "("*g~O[2]~m^{-2}~day^{-1}*")"))+
             # scale_x_continuous("Day of Year", 
             #                    limits=c(151,201), breaks=c(160,176,192))+
             theme_base+
             theme(legend.position = "top",
                   legend.key.width = unit(3, "line"),
                   legend.key.height = unit(0.8, "line"),
                   legend.text.align = 0.5)
         }


model_fit %>%
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
                         linetype=0, fill = "gray80")+
             geom_line(aes(color=name), size=0.6)+
             scale_linetype_manual("",
                                   values=c(1,5,2),
                                   guide = guide_legend(direction = "horizontal",
                                                        label.position = "top"))+
             scale_y_continuous(expression(Metabolism~
                                             "("*g~O[2]~m^{-2}~day^{-1}*")"))+
             scale_color_manual(values=c("black","gray50","white"))+
             # scale_x_continuous("Day of Year", 
             #                    limits=c(151,201), breaks=c(160,176,192))+
             theme_base+
             theme(legend.position = "top",
                   legend.key.width = unit(3, "line"),
                   legend.key.height = unit(0.8, "line"),
                   legend.text.align = 0.5)
         }


model_fit %>%
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
             geom_ribbon(aes(ymin=lower16, ymax=upper84, fill = name),
                         linetype=0, alpha =0.3)+
             geom_line(aes(color=name), size=0.6)+
             scale_linetype_manual("",
                                   values=c(1,5,2),
                                   guide = guide_legend(direction = "horizontal",
                                                        label.position = "top"))+
             scale_y_continuous(expression(Metabolism~
                                             "("*g~O[2]~m^{-2}~day^{-1}*")"))+
             scale_color_manual(values=c("dodgerblue","black","firebrick2"))+
             scale_fill_manual(values=c("dodgerblue","black","firebrick2"))+
             # scale_x_continuous("Day of Year", 
             #                    limits=c(151,201), breaks=c(160,176,192))+
             theme_base+
             theme(legend.position = "top",
                   legend.key.width = unit(3, "line"),
                   legend.key.height = unit(0.8, "line"),
                   legend.text.align = 0.5)
         }
