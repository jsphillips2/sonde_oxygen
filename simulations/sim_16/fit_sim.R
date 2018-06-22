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
    gamma_1 = runif(1, 1, 2),
    gamma_2 = runif(1, 1, 2),
    sig_beta = runif(1, 0, 2),
    sig_alpha = runif(1, 0, 2),
    sig_rho = runif(1, 0, 2),
    sig_proc = runif(1, 50, 150),
    log_beta0 = runif(data$D, log(0.5) + 5.5, log(1.5) + 5.5),
    log_alpha = runif(data$D, log(0.5) + 1, log(1.5) + 1),
    log_rho = runif(data$D, log(0.5) + 5.5, log(1.5) + 5.5)
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
seed = 6

# fit model
fit = stan(file=model_path, data=data_f, seed=seed, chains = chains,
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
#========== Export 
#==========

# clean variable names in summary
fit_clean = fit_summary %>%
  rename(lower16 = `16%`, middle = `50%`, upper84 = `84%`)  %>%
  mutate(name = strsplit(var, "\\[|\\]|,") %>% map_chr(~.x[1]),
         index = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2])),
         day = ifelse(name %in% c("beta0","alpha","rho","GPP","ER","NEP","AIR","Flux"), index, data$D_M[index])) %>%
  select(name, index, day, middle, lower16, upper84) %>%
  filter(!(name %in% c("log_beta0","log_rho","lp__"))) %>%
  mutate(seed = seed)

# Export
output_path = paste0("simulations/sim_16/sim_fits/sim_a/summary_clean","_",seed,".csv")
# write_csv(fit_clean, output_path)



  

