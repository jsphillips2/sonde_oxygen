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
data = read_rdump("simulations/data_16/sonde_list_16.R")
model_fit = read_csv("simulations/fit_16/model_output_16/summary_clean.csv")

# add fitted values to data
vars = c("beta","alpha","rho","gamma_1","gamma_2","k0","k1","sig_proc")
vars_f =  vars %>%
  lapply(function(x){
    nm = paste0(x,"_f")
    data$nm = {model_fit %>% filter(name==x)}$middle
  })
names(vars_f) = vars %>% lapply(function(x){paste0(x,"_f")})
data_f = data %>% append(vars_f)

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
model = "o2_sim"
model_path = paste0("simulations/sim_16/",model,".stan")
chains = 1
iter = 1000

# fit model
fit = stan(file=model_path, data=data_f, seed=1, chains = chains,
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

