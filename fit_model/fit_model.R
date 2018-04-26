#==========
#========== Preliminaries
#==========

# load packages
library(tidyverse)
library(rstan)
source("fit_model/stan_utility.R")
library(GGally)

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
model = "o2_model_nc3"
model_path = paste0("fit_model/",model,".stan")
chains = 1
iter = 1000

# fit model
fit = stan(file=model_path, data=data, seed=194838, chains = chains,
           init = init_fn, iter = iter)

# summary of fit
fit_summary = summary(fit)$summary %>% 
{as_data_frame(.) %>%
    mutate(var = rownames(summary(fit)$summary))}

# check Rhat
fit_summary %>% filter(Rhat > 1.05) %>% select(Rhat, n_eff, var)

# additional diagnostics
check_div(fit)
check_treedepth(fit)
check_energy(fit)





#==========
#========== Additional diagnostics
#==========

# fixed parameters by step
fixed_par_v = c("alpha","gamma_1","gamma_2","sig_beta0","sig_rho","sig_proc","lp__")
fixed_pars = rstan::extract(fit, pars=fixed_par_v) %>%
  lapply(as_data_frame) %>%
  bind_cols()
names(fixed_pars) = fixed_par_v

# examine chains for parameters
fixed_pars %>%
  mutate(chain = rep(chains, each = iter/2), step = rep(c(1:(iter/2)), chains)) %>%
  gather(par, value, names(fixed_pars)) %>%
  filter(par != "lp__") %>%
  ggplot(aes(step, value, color=factor(chain)))+
  facet_wrap(~par, scales="free_y")+
  geom_line(alpha=0.5)+
  theme_bw()

# pairs plot for parameters
ggpairs(fixed_pars)





#==========
#========== Poterior Predictive Check
#==========

post_pred_v = c("chi_proc_real","chi_proc_sim","chi_obs_real","chi_obs_sim")
post_pred = rstan::extract(fit, pars=post_pred_v) %>%
  lapply(as_data_frame) %>%
  bind_cols()
names(post_pred) = post_pred_v

post_pred %>%
  ggplot(aes(chi_proc_real,chi_proc_sim))+
  geom_point()+
  geom_abline(intercept=0, slope=1)+
  scale_y_continuous(limits=c(4800,5600))+
  scale_x_continuous(limits=c(4800,5600))+
  theme_bw()

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

# clean variable names in summary
fit_clean = fit_summary %>%
  rename(lower2 = `2.5%`, lower25 = `25%`, middle = `50%`, upper75 = `75%`, upper97 = `97.5%`)  %>%
  mutate(name = strsplit(var, "\\[|\\]|,") %>% map_chr(~.x[1]),
         index = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2])),
         day = ifelse(name %in% c("beta0","rho"), index, data$D_M[index])) %>%
  select(name, index, day, middle, lower2, lower25, upper75, upper97) %>%
  filter(!(name %in% c("log_beta0","log_rho","lp__")))

# Export
output_path = paste0("model_output/",model)
# write_csv(fixed_pars, paste0(output_path,"/fixed_pars_full.csv"))
# write_csv(fit_clean, paste0(output_path,"/summaries_clean.csv"))



