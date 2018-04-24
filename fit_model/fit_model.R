#==========
#========== Preliminaries
#==========

# load packages
library(tidyverse)
library(rstan)
source("fit_model/stan_utility.R")

# stan settings
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

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
model = "o2_model_nc"
model_path = paste0("fit_model/",model,".stan")

# fit model
fit = stan(file='fit_model/o2_model_nc.stan', data=data, seed=194838, chains = 1,
           init = init_fn, iter = 2000)

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
#========== Process output
#==========

# fixed parameters by step
fixed_par_v = c("alpha","gamma_1","gamma_2","sig_beta0","sig_rho","sig_proc","lp__")
fixed_pars = rstan::extract(fit, pars=fixed_par_v) %>%
  lapply(as_data_frame) %>%
  bind_cols()
names(fixed_pars) = fixed_par_v

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



