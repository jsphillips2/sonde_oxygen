#==========
#========== Preliminaries
#==========

# load packages
library(tidyverse)
library(rstan)
library(GGally)
source("main_analysis/fit_model/stan_utility.R")

# stan settings
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()-1)

# read data
data = read_rdump("data/sonde_list.R")

# function for initial values
init_fn = function(){
  list(
    alpha = runif(1, 1, 5),
    gamma_1 = runif(1, 1, 2),
    gamma_2 = runif(1, 1, 2),
    sig_beta = runif(1, 0, 2),
    sig_alpha = runif(1, 0, 2),
    sig_rho = runif(1, 0, 2),
    sig_proc = runif(1, 50, 150),
    log_beta0 = runif(data$D, log(0.75) + 6, log(1.25) + 6),
    log_rho = runif(data$D, log(0.75) + 6, log(1.25) + 6)
  )
}





#==========
#========== Fit model
#==========

# initial specifications
model = "o2_model"
model_path = paste0("main_analysis/fit_model/",model,".stan")
chains = 1
iter = 1000

# fit model
fit = stan(file=model_path, data=data, seed=194838, chains = chains,
           init = init_fn, iter = iter)

# summary of fit
fit_summary = summary(fit, probs=c(0.16, 0.5, 0.84))$summary %>% 
{as_data_frame(.) %>%
    mutate(var = rownames(summary(fit)$summary))}

# check Rhat
fit_summary %>% filter(Rhat > 1.05) %>% select(Rhat, n_eff, var)

# additional diagnostics
check_div(fit)
check_treedepth(fit)
check_energy(fit)





#==========
#========== Examine Chains
#==========

# fixed parameters by step
fixed_par_v = c("gamma_1","gamma_2","k0","k1","sig_beta0","sig_alpha","sig_rho","sig_proc","lp__")
fixed_pars = rstan::extract(fit, pars=fixed_par_v) %>%
  lapply(as_data_frame) %>%
  bind_cols() %>%
  mutate(chain = rep(chains, each = iter/2), step = rep(c(1:(iter/2))))
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





#==========
#========== Poterior Predictive Check
#==========

# extract relevant variables
post_pred_v = c("chi_proc_real","chi_proc_sim","chi_obs_real","chi_obs_sim")
post_pred = rstan::extract(fit, pars=post_pred_v) %>%
  lapply(as_data_frame) %>%
  bind_cols() %>%
  mutate(chain = rep(chains, each = iter/2), step = rep(c(1:(iter/2))))
names(post_pred) = c(post_pred_v,"chain","step")

# process error
post_pred %>%
  ggplot(aes(chi_proc_real,chi_proc_sim))+
  geom_point()+
  geom_abline(intercept=0, slope=1)+
  scale_y_continuous(limits=c(4800,5600))+
  scale_x_continuous(limits=c(4800,5600))+
  theme_bw()

# observation error
post_pred %>%
  ggplot(aes(chi_obs_real,chi_obs_sim))+
  geom_point()+
  geom_abline(intercept=0, slope=1)+
  scale_y_continuous(limits=c(4800,5600))+
  scale_x_continuous(limits=c(4800,5600))+
  theme_bw()
  



#==========
#========== Prepare output for export
#==========

# beta0 and rho full
daily_pars = c("beta0","alpha","rho")
daily = rstan::extract(fit, pars=daily_pars) %>% 
{lapply(1:3, function(x){y = .[[x]] %>% as_data_frame %>%
  mutate(chain = rep(chains, each = iter/2), step = rep(c(1:(iter/2)))) %>%
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
         day = ifelse(name %in% c("beta0","alpha","rho","GPP","ER","NEP","AIR","Flux"), index, data$D_M[index])) %>%
  select(name, index, day, middle, lower16, upper84) %>%
  filter(!(name %in% c("log_beta0","log_rho","lp__")))

# Export
output_path = "main_analysis/model_output"
# write_csv(fixed_pars, paste0(output_path,"/fixed_pars_full.csv"))
# write_csv(post_pred, paste0(output_path,"/post_pred_full.csv"))
# write_csv(daily, paste0(output_path,"/daily_full.csv"))
# write_csv(fit_clean, paste0(output_path,"/summary_clean.csv"))



