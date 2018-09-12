#==========
#========== Preliminaries
#==========

# load packages
library(tidyverse)

# import specifications
input_path = "analyses/test_analysis/simulation/"
types = list.files(paste0(input_path,"output/")) 

# import parameter estimates
pars_fit = lapply(types, function(x){
  read_csv(paste0(input_path, "output/", x, "/fixed_pars_full.csv")) %>%
    cbind(read_csv(paste0(input_path, "input/", x, "/type_data.csv")) %>%
            select(var, value) %>%
            spread(var, value)) %>%
    tbl_df %>%
    mutate(type = x)
}) %>%
  bind_rows

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
#========== Preliminaries
#==========

sig_pars = pars_fit %>%
  select(sig_beta0, sig_alpha, sig_rho, fix_beta0, fix_alpha, fix_rho) 

# sig_beta0
sig_pars %>%
  mutate(fix_beta0 = ifelse(fix_beta0 == T, "Fixed", "Not Fixed"),
         fix_alpha = ifelse(fix_alpha == T, "alpha:~Fixed", "alpha:~Not~Fixed"),
         fix_rho = ifelse(fix_rho == T, "rho:~Fixed", "rho:~Not~Fixed")) %>%
  ggplot(aes(sig_beta0, linetype = fix_beta0))+
  facet_grid(fix_rho ~ fix_alpha, labeller = label_parsed)+
  stat_density(geom = "line", position = "identity")+
  scale_linetype_manual(expression(beta[0]), values = c(2,1))+
  scale_y_continuous("Probability Density", breaks=NULL)+
  scale_x_continuous(expression(sigma[beta[0]]), breaks = c(0.03, 0.09, 0.15))+
  theme_base+
  theme(legend.position = c(0.85,0.85))

# sig_alpha
sig_pars %>%
  mutate(fix_beta0 = ifelse(fix_beta0 == T, 
                            "beta[0]:~Fixed", "beta[0]:~Not~Fixed"),
         fix_alpha = ifelse(fix_alpha == T, "Fixed", "Not Fixed"),
         fix_rho = ifelse(fix_rho == T, "rho:~Fixed", "rho:~Not~Fixed")) %>%
  ggplot(aes(sig_alpha, linetype = fix_alpha))+
  facet_grid(fix_rho ~ fix_beta0, labeller = label_parsed)+
  stat_density(geom = "line", position = "identity")+
  scale_linetype_manual(expression(alpha), values = c(2,1))+
  scale_y_continuous("Probability Density", breaks=NULL)+
  scale_x_continuous(expression(sigma[alpha]), breaks = c(0.07, 0.22, 0.37))+
  theme_base+
  theme(legend.position = c(0.85,0.85))

# sig_rho
sig_pars %>%
  mutate(fix_beta0 = ifelse(fix_beta0 == T, 
                            "beta[0]:~Fixed", "beta[0]:~Not~Fixed"),
         fix_alpha = ifelse(fix_alpha == T, "alpha:~Fixed", "alpha:~Not~Fixed"),
         fix_rho = ifelse(fix_rho == T, "Fixed", "Not Fixed")) %>%
  ggplot(aes(sig_rho, linetype = fix_rho))+
  facet_grid(fix_beta0 ~ fix_alpha, labeller = label_parsed)+
  stat_density(geom = "line", position = "identity")+
  scale_linetype_manual(expression(rho), values = c(2,1))+
  scale_y_continuous("Probability Density", breaks=NULL)+
  scale_x_continuous(expression(sigma[rho]), breaks = c(0.03, 0.10, 0.17))+
  theme_base+
  theme(legend.position = c(0.85,0.85))


