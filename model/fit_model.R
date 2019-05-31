#==========
#========== Preliminaries
#==========

# load packages
library(tidyverse)
library(rstan)
library(loo)
source("model/stan_utility.R")

# stan settings
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()-2)

# specify analysis
analysis <- "main"

# read data
data <- read_rdump(paste0("model/input/",analysis,"/sonde_list.R"))

# set observation error 
if(analysis!="sig_obs"){data$sig_obs <- 0.01}

# set reference temperature
# data_full$temp %>% mean()
data$temp_ref <- 12





#==========
#========== Fit model
#==========

# model
model <- if(analysis=="main"|analysis=="surface_par"|analysis=="alt_k"){"o2_model.stan"
} else if(analysis=="fixed"){"o2_model_fixed.stan"
    } else if(analysis=="sig_obs"){"o2_model_sig_obs.stan"}

# model path
model_path <- paste0("model/stan/",model)

# MCMC specificaiotns
chains <- 4
iter <- 2000
adapt_delta <- 0.8
max_treedepth <- 10

# fit model
fit <- stan(file = model_path, data = data, seed=1, chains = chains, iter = iter, 
            control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth))

# summary of fit
fit_summary <- summary(fit, probs=c(0.16, 0.5, 0.84))$summary %>% 
{as_tibble(.) %>%
    mutate(var = rownames(summary(fit)$summary))}

# check Rhat & n_eff
fit_summary %>% filter(Rhat > 1.01) %>% select(Rhat, n_eff, var) %>% arrange(-Rhat)
fit_summary %>% filter(n_eff < 0.5*(chains*iter/2)) %>% select(Rhat, n_eff, var) %>% arrange(n_eff) %>%
  mutate(eff_frac = n_eff/(chains*iter/2))

# additional diagnostics
check_div(fit)
check_treedepth(fit,max_treedepth)
check_energy(fit)

# export path
output_path <- paste0("model/output/",analysis)

# save model full output
saveRDS(fit, paste0(output_path,"/fit.rds"))





#==========
#========== Examine Chains
#==========

# function for selecting fixed parameters
fixed_par_fn <- function(x){
  if(x=="o2_model.stan"){return(c("gamma_1","gamma_2","sig_b0","sig_a","sig_r","sig_proc","lp__"))}
  if(x=="o2_model_sig_obs.stan"){return(c("gamma_1","gamma_2","sig_b0","sig_a","sig_r",
                                          "sig_proc","sig_obs","lp__"))}
  if(x=="o2_model_fixed.stan"){return(c("gamma_1","gamma_2","sig_proc","lp__"))}
}

# fixed parameters by step
fixed_pars <- rstan::extract(fit, pars=fixed_par_fn(model)) %>%
  lapply(as_tibble) %>%
  bind_cols() %>%
  set_names(fixed_par_fn(model)) %>%
  mutate(chain = rep(1:chains, each = iter/2), step = rep(c(1:(iter/2)), chains))

# examine chains for parameters
fixed_pars %>%
  gather(par, value, -chain, -step) %>%
  filter(par != "lp__") %>%
  ggplot(aes(step, value, color=factor(chain)))+
  facet_wrap(~par, scales="free_y")+
  geom_line(alpha=0.5)+
  theme_bw()

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
daily_pars <- c("beta0","alpha","rho")
daily <- rstan::extract(fit, pars=daily_pars) %>% 
{lapply(1:length(daily_pars), function(x){y = .[[x]] %>% as_tibble %>%
  mutate(chain = rep(1:chains, each = iter/2), step = rep(c(1:(iter/2)), chains)) %>%
  gather(var, value, -chain, -step) %>%
  mutate(day = strsplit(var, "V") %>% map_int(~as.integer(.x[2])),
         name = daily_pars[x]) %>%
  select(name, chain, step, day, value)
return(y)
})} %>%
  bind_rows()

# clean variable names in summary
fit_clean <- fit_summary %>%
  rename(lower16 = `16%`, middle = `50%`, upper84 = `84%`)  %>%
  mutate(name = strsplit(var, "\\[|\\]|,") %>% map_chr(~.x[1]),
         index = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2])),
         day = ifelse(name %in% c("beta0","alpha","rho","GPP","ER","NEP","AIR","Flux"), 
                      index, 
                      data$map_days[index])) %>%
  select(name, index, day, middle, lower16, upper84)





#==========
#==========  Likelihoods & LOO
#==========

# log-likelihoods
likelihood <- tibble(analysis = analysis,
                     mean_log_lik = {fit_summary %>%
    filter(str_detect(fit_summary$var, "lik"))}$mean %>% sum(),
    median_log_lik = {fit_summary %>%
        filter(str_detect(fit_summary$var, "lik"))}$`50%` %>% sum())

# extract log Likelihoods and exponentiate
log_lik <- extract_log_lik(fit, merge_chains = FALSE)
r_eff <- relative_eff(exp(log_lik)) 

# LOO
loo <- loo(log_lik, r_eff = r_eff, cores = parallel::detectCores()-2)





#==========
#==========  Export
#==========

# export
# write_csv(fixed_pars, paste0(output_path,"/fixed_pars_full.csv"))
# write_csv(daily, paste0(output_path,"/daily_full.csv"))
# write_csv(fit_clean, paste0(output_path,"/summary_clean.csv"))
# write_csv(likelihood, paste0(output_path,"/log_liks.csv"))
# saveRDS(loo, paste0(output_path,"/loo.rds"))




