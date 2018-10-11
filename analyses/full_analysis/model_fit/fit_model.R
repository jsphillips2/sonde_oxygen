#==========
#========== Preliminaries
#==========

# load packages
library(tidyverse)
library(rstan)
library(truncnorm)
library(GGally)
source("model/stan_utility.R")

# stan settings
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()-2)

# read data
data = read_rdump("analyses/full_analysis/model_fit/input/sonde_list.R")

# set observation error 
data$sig_obs = 0.01

# set reference temperature
# data_full$temp %>% mean()
data$temp_ref = 12





#==========
#========== Fit model
#==========

# initial specifications
model = "o2_model"
model_path = paste0("model/",model,".stan")
chains = 4
iter = 2000

# fit model
fit = stan(file = model_path, data = data, seed=1, chains = chains, iter = iter)

# summary of fit
fit_summary = summary(fit, probs=c(0.16, 0.5, 0.84))$summary %>% 
{as_data_frame(.) %>%
    mutate(var = rownames(summary(fit)$summary))}

# check Rhat & n_eff
fit_summary %>% filter(Rhat > 1.05) %>% select(Rhat, n_eff, var) %>% arrange(-Rhat)
fit_summary %>% filter(n_eff < 0.5*(chains*iter/2)) %>% select(Rhat, n_eff, var) %>% arrange(n_eff)

# additional diagnostics
check_div(fit)
check_treedepth(fit)
check_energy(fit)





#==========
#========== Examine Chains
#==========

# fixed parameters by step
fixed_par_v = c("gamma_1","gamma_2","sig_b0","sig_a","sig_r","sig_proc","sig_obs","lp__")
fixed_pars = rstan::extract(fit, pars=fixed_par_v) %>%
  lapply(as_data_frame) %>%
  bind_cols() %>%
  mutate(chain = rep(1:chains, each = iter/2), step = rep(c(1:(iter/2)), chains))
names(fixed_pars) = c(fixed_par_v,"chain","step")

# examine chains for parameters
fixed_pars %>%
  gather(par, value, -chain, -step) %>%
  filter(par != "lp__") %>%
  ggplot(aes(step, value, color=factor(chain)))+
  facet_wrap(~par, scales="free_y")+
  geom_line(alpha=0.5)+
  theme_bw()

# pairs plot for parameters
ggpairs(fixed_pars %>% select(-chain, -step))

# posterior densities
fixed_pars %>%
  gather(par, value, -chain, -step) %>%
  filter(par != "lp__") %>%
  ggplot(aes(value))+
  facet_wrap(~par, scales="free")+
  stat_density(alpha=0.5, geom = "line")+
  theme_bw()





#==========
#========== Prepare output for export
#==========

# beta0 and rho full
daily_pars = c("beta0","alpha","rho")
daily = rstan::extract(fit, pars=daily_pars) %>% 
{lapply(1:length(daily_pars), function(x){y = .[[x]] %>% as_data_frame %>%
  mutate(chain = rep(1:chains, each = iter/2), step = rep(c(1:(iter/2)), chains)) %>%
  gather(var, value, -chain, -step) %>%
  mutate(day = strsplit(var, "V") %>% map_int(~as.integer(.x[2])),
         name = daily_pars[x]) %>%
  select(name, chain, step, day, value)
return(y)
})} %>%
  bind_rows()

# clean variable names in summary
fit_clean = fit_summary %>%
  rename(lower16 = `16%`, middle = `50%`, upper84 = `84%`)  %>%
  mutate(name = strsplit(var, "\\[|\\]|,") %>% map_chr(~.x[1]),
         index = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2])),
         day = ifelse(name %in% c("beta0","alpha","rho","GPP","ER","NEP","AIR","Flux"), 
                      index, 
                      data$map_days[index])) %>%
  select(name, index, day, middle, lower16, upper84) %>%
  filter(!(name %in% c("log_beta0","log_rho","lp__")))

# export path
output_path = "analyses/full_analysis/model_fit/output"
# output_path = "analyses/full_analysis/model_fit/output/sig_obs"

# export
write_csv(fixed_pars, paste0(output_path,"/fixed_pars_full.csv"))
write_csv(daily, paste0(output_path,"/daily_full.csv"))
write_csv(fit_clean, paste0(output_path,"/summary_clean.csv"))



