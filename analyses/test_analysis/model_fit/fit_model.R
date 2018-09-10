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
data = read_rdump("analyses/test_analysis/model_fit/input/sonde_list.R")





#==========
#========== Priors & initial values
#==========

# priors
priors = list(k0_prior = c(2, 1),
              k1_prior = c(0.2, 0.1),
              gamma_1_prior = c(1.1, 1),
              gamma_2_prior = c(1.1, 1),
              sig_beta0_prior = c(0, 1),
              sig_alpha_prior = c(0, 1),
              sig_rho_prior = c(0, 1),
              sig_proc_prior = c(100, 100),
              log_beta0_prior = c(5, 7),
              log_alpha_prior = c(2, 4),
              log_rho_prior = c(5, 7))

# examine priors (change name and bounds as appropriate)
# priors$k0_prior %>%
#   {rtruncnorm(n = 50000, a = 0, b = Inf, mean = .[1], sd = .[2])} %>%
#   hist

# export priors
as_data_frame(append(list(par = c("mean","sd")), priors)) %>%
  gather(name, value, 2:12) %>%
  spread(par, value) %>%
  write_csv("analyses/test_analysis/model_fit/input/priors.csv")

# add priors to data
data_full = data %>% append(priors)
  
# function for initial values
init_fn = function(){
  list(k0 = runif(n = 1, min = 1, max = 3),
       k1 = runif(n = 1, min = 0.1, max = 0.3),
       gamma_1 = runif(n = 1, min = 1, max = 2),
       gamma_2 = runif(n = 1, min = 1, max = 2),
       sig_beta = runif(n = 1, min = 0, max = 0.2),
       sig_alpha = runif(n = 1, min = 0, max = 0.2),
       sig_rho = runif(n = 1, min = 0, max = 0.2),
       sig_proc = runif(n = 1, min = 50, max = 150),
       log_beta0_init = runif(n = data$Y + 1, min = log(0.1) + 5, 
                           max = log(1.9) + 5),
       log_alpha_init = runif(n = data$Y + 1, min = log(0.1) + 2, 
                           max = log(1.9) + 2),
       log_rho_init = runif(n = data$Y + 1, min = log(0.1) + 5, 
                         max = log(1.9) + 5)
  )
}





#==========
#========== Fit model
#==========

# initial specifications
model = "o2_model"
model_path = paste0("model/",model,".stan")
chains = 4
iter = 2000

# fit model
fit = stan(file = model_path, data = data_full, seed=1, chains = chains,
           init = init_fn, iter = iter)

# summary of fit
fit_summary = summary(fit, probs=c(0.16, 0.5, 0.84))$summary %>% 
{as_data_frame(.) %>%
    mutate(var = rownames(summary(fit)$summary))}

# check Rhat & n_eff
# note that initial values beta0, alpha, and rho sample for an extra year
# these extra values do not contribute to the likelihood
fit_summary %>% filter(Rhat > 1.05) %>% select(Rhat, n_eff, var)
fit_summary %>% filter(n_eff < 0.5*(iter/2)) %>% select(Rhat, n_eff, var)

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
# ggpairs(fixed_pars %>% select(-chain, -step))





#==========
#========== Poterior Predictive Check
#==========

# extract relevant variables
post_pred_v = c("chi_proc_real","chi_proc_sim","chi_obs_real","chi_obs_sim")
post_pred = rstan::extract(fit, pars=post_pred_v) %>%
  lapply(as_data_frame) %>%
  bind_cols() %>%
  mutate(chain = rep(1:chains, each = iter/2), step = rep(c(1:(iter/2)), chains))
names(post_pred) = c(post_pred_v,"chain","step")

# process error
post_pred %>%
  ggplot(aes(chi_proc_real,chi_proc_sim))+
  geom_point()+
  geom_abline(intercept=0, slope=1)+
  theme_bw()

# observation error
post_pred %>%
  ggplot(aes(chi_obs_real,chi_obs_sim))+
  geom_point()+
  geom_abline(intercept=0, slope=1)+
  theme_bw()
  



#==========
#========== Prepare output for export
#==========

# beta0 and rho full
daily_pars = c("beta0","alpha","rho")
daily = rstan::extract(fit, pars=daily_pars) %>% 
{lapply(1:3, function(x){y = .[[x]] %>% as_data_frame %>%
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
         day = ifelse(name %in% c("beta0","alpha","rho","GPP","ER","NEP","AIR","Flux"), index, data$D_M[index])) %>%
  select(name, index, day, middle, lower16, upper84) %>%
  filter(!(name %in% c("log_beta0","log_rho","lp__")))

# export
output_path = "analyses/test_analysis/model_fit/output"
# write_csv(fixed_pars, paste0(output_path,"/fixed_pars_full.csv"))
# write_csv(post_pred, paste0(output_path,"/post_pred_full.csv"))
# write_csv(daily, paste0(output_path,"/daily_full.csv"))
# write_csv(fit_clean, paste0(output_path,"/summary_clean.csv"))



