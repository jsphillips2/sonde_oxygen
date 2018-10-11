#==========
#========== Preliminaries
#==========

# load packages
library(tidyverse)
library(matrixcalc)

# import data and model fit
input_dir = "analyses/full_analysis/model_fit/"
sonde_data = read_csv(paste0(input_dir,"input/sonde_prep.csv"))
params_full = read_csv(paste0(input_dir,"output/sig_obs10/fixed_pars_full.csv"))
daily = read_csv(paste0(input_dir,"output/sig_obs10/daily_full.csv"))





#==========
#========== Prepare Data
#==========

# combine data
comb_data = daily %>%
  spread(name, value) %>%
  full_join(params_full) %>%
  left_join(sonde_data %>%
              rename(day = unique_day) %>%
              select(year,yday,day,hour,par_int,temp) %>%
              na.omit()) %>%
  select(chain,step,year,yday,hour,par_int,temp,beta0,alpha,rho,alpha,gamma_1,gamma_2) %>%
  mutate(temp_ref = mean(temp))

years = unique(sonde_data$year)





#==========
#========== Define functions
#==========

# GPP
gpp_part_fun = function(y, x){
  
  # Select year from data frame of inputs
  x_filt = x %>% filter(year %in% y) 
  
  # Vector of mean values
  gpp_par_mean = x_filt %>%
    select(beta0, alpha, gamma_1, par_int_0, par_int_z, temp_0, temp_z, temp_ref) %>%
    summarise_all(mean)
  
  # Define model function
  gpp_fun = function(beta0,alpha,gamma_1,par_int_0,par_int_z,temp_0,temp_z,temp_ref){
    beta0*gamma_1^((temp_0 + temp_z)-temp_ref)*tanh((alpha/(beta0*gamma_1^((temp_0 + temp_z)-temp_ref)))*(par_int_0 + par_int_z))
  }
  
  # Evaluate at mean value
  gpp_mean = gpp_par_mean %>% 
  {gpp_fun(.$beta0,.$alpha,.$gamma_1,.$par_int_0,.$par_int_z,.$temp_0,.$temp_z,.$temp_ref)}
  
  # Define model gradient
  gpp_grad_fun <-deriv(~beta0*gamma_1^((temp_0 + temp_z)-temp_ref)*tanh((alpha/(beta0*gamma_1^((temp_0 + temp_z)-temp_ref)))*(par_int_0 + par_int_z)), 
                       c("beta0","alpha","gamma_1","par_int_0","par_int_z","temp_0","temp_z","temp_ref"), 
                       function(beta0,alpha,gamma_1,par_int_0,par_int_z,temp_0,temp_z,temp_ref){})
  gpp_grad = t(attributes(gpp_par_mean %>% {gpp_grad_fun(.$beta0,.$alpha,.$gamma_1,.$par_int_0,.$par_int_z,.$temp_0,.$temp_z,.$temp_ref)})$gradient)
  
  # Define variance-covariance matrix
  gpp_vcv = cov(x_filt %>% select(beta0,alpha,gamma_1,par_int_0,par_int_z,temp_0,temp_z,temp_ref)) 
  
  # Calcuate % contributions to variance using delta method
  gpp_mat = hadamard.prod(gpp_vcv, gpp_grad %*% t(gpp_grad))
  gpp_cont = rowSums(gpp_mat)/sum(gpp_mat)
  
  # Extract contribution of variance and covariance terms
  var_cont = diag(gpp_mat)
  gpp_cov = gpp_mat
  diag(gpp_cov) = 0
  cov_cont = colSums(gpp_cov)
  
  # Store in data frame
  yy = data_frame(param = names(gpp_cont), 
                  mean = t(gpp_par_mean)[,1],
                  sd = sqrt(diag(gpp_vcv)),
                  cv = sd/mean,
                  sen = gpp_grad[,1],
                  ela = round((mean/abs(gpp_mean))*sen,2),
                  cont100 = gpp_cont) %>% 
    filter(sd > 0)
  return(yy)
}




# er
er_part_fun = function(y, x){
  
  # Select year from data frame of inputs
  x_filt = x %>% filter(year %in% y) 
  
  # Vector of mean values
  er_par_mean = x_filt %>%
    select(rho, gamma_2, temp_0, temp_z, temp_ref) %>%
    summarise_all(mean)
  
  # Define model function
  er_fun = function(rho, gamma_2, temp_0, temp_z, temp_ref){
    rho*gamma_2^((temp_0 + temp_z)-temp_ref)
  }
  
  # Evaluate at mean value
  er_mean = er_par_mean %>% 
  {er_fun(.$rho, .$gamma_2, .$temp_0, .$temp_z, .$temp_ref)}
  
  # Define model gradient
  er_grad_fun <-deriv(~rho*gamma_2^((temp_0 + temp_z)-temp_ref), c("rho", "gamma_2", "temp_0", "temp_z", "temp_ref"), 
                        function(rho, gamma_2, temp_0, temp_z, temp_ref){})
  er_grad = t(attributes(er_par_mean %>% {er_grad_fun(.$rho, .$gamma_2, .$temp_0, .$temp_z, .$temp_ref)})$gradient)
  
  # Define variance-covariance matrix
  er_vcv = cov(x %>% select(rho, gamma_2, temp_0, temp_z, temp_ref)) 
  
  # Calcuate % contributions to variance using delta method
  er_mat = hadamard.prod(er_vcv, er_grad %*% t(er_grad))
  er_cont = rowSums(er_mat)/sum(er_mat)
  
  # Extract contribution of variance and covariance terms
  var_cont = diag(er_mat)
  er_cov = er_mat
  diag(er_cov) = 0
  cov_cont = colSums(er_cov)
  
  # Store in data frame
  yy = data_frame(param = names(er_cont), 
                  mean = t(er_par_mean)[,1],
                  sd = sqrt(diag(er_vcv)),
                  cv = sd/mean,
                  sen = er_grad[,1],
                  ela = round((mean/abs(er_mean))*sen,2),
                  cont100 = er_cont) %>% 
    filter(sd > 0)
  return(yy)
}



# NEP
nep_part_fun = function(y, x){
  
  # Select year from data frame of inputs
  x_filt = x %>% filter(year %in% y) 
  
  # Vector of mean values
  nep_par_mean = x_filt %>%
    select(beta0, alpha, rho, gamma_1, gamma_2, par_int_0, par_int_z, temp_0, temp_z, temp_ref) %>%
    summarise_all(mean)
  
  # Define model function
  nep_fun = function(beta0,alpha,rho,gamma_1,gamma_2,par_int_0,par_int_z,temp_0,temp_z,temp_ref){
    beta0*gamma_1^((temp_0 + temp_z)-temp_ref)*tanh((alpha/(beta0*gamma_1^((temp_0 + temp_z)-temp_ref)))*(par_int_0 + par_int_z)) 
    - rho*gamma_2^((temp_0 + temp_z)-temp_ref)
  }
  
  # Evaluate at mean value
  nep_mean = nep_par_mean %>% 
  {nep_fun(.$beta0,.$alpha,.$rho,.$gamma_1,.$gamma_2,.$par_int_0,.$par_int_z,.$temp_0,.$temp_z,.$temp_ref)}
  
  # Define model gradient
  nep_grad_fun <-deriv(~beta0*gamma_1^((temp_0 + temp_z)-temp_ref)*tanh((alpha/(beta0*gamma_1^((temp_0 + temp_z)-temp_ref)))*(par_int_0 + par_int_z)) 
                       - rho*gamma_2^((temp_0 + temp_z)-temp_ref), 
                       c("beta0","alpha","rho","gamma_1","gamma_2","par_int_0","par_int_z","temp_0","temp_z","temp_ref"), 
                       function(beta0,alpha,rho,gamma_1,gamma_2,par_int_0,par_int_z,temp_0,temp_z,temp_ref){})
  nep_grad = t(attributes(nep_par_mean %>% {nep_grad_fun(.$beta0,.$alpha,.$rho,.$gamma_1,.$gamma_2,.$par_int_0,.$par_int_z,.$temp_0,.$temp_z,.$temp_ref)})$gradient)
  
  # Define variance-covariance matrix
  nep_vcv = cov(x_filt %>% select(beta0,alpha,rho,gamma_1,gamma_2,par_int_0,par_int_z,temp_0,temp_z,temp_ref)) 
  
  # Calcuate % contributions to variance using delta method
  nep_mat = hadamard.prod(nep_vcv, nep_grad %*% t(nep_grad))
  nep_cont = rowSums(nep_mat)/sum(nep_mat)
  
  # Extract contribution of variance and covariance terms
  var_cont = diag(nep_mat)
  nep_cov = nep_mat
  diag(nep_cov) = 0
  cov_cont = colSums(nep_cov)
  
  # Store in data frame
  yy = data_frame(param = names(nep_cont), 
                  mean = t(nep_par_mean)[,1],
                  sd = sqrt(diag(nep_vcv)),
                  cv = sd/mean,
                  sen = nep_grad[,1],
                  ela = round((mean/abs(nep_mean))*sen,2),
                  cont100 = nep_cont) %>% 
    filter(sd > 0)
  return(yy)
}





#==========
#========== GPP (pool years)
#==========

# the functions are designed for separate with vs. between day variation in par_int/temp
# to collapose these into a single variable, set par_int_0 = par_int and par_int_z = 0, etc.
comb_data_trans = comb_data %>% 
  filter(chain==1) %>%
  group_by(chain, step, year, yday) %>%
  mutate(par_int_0 = par_int,
         par_int_z = par_int - par_int_0,
         temp_0 = temp,
         temp_z = temp - temp_0) %>%
  select(-par_int, -temp) %>%
  ungroup

# partition for each step 
gpp_part = comb_data_trans$step %>% 
  unique %>%
  lapply(function(x){
    d = comb_data_trans %>% filter(step==x)
    gpp_part_fun(years, d)
  }) %>%
  bind_rows()

# partition for each step 
gpp_part_sum = lapply(c(0.16,0.5,0.84), function(x){
  gpp_part %>%
    group_by(param) %>%
    summarize_all(funs(quantile), probs = x) %>%
    mutate(quant = x)
})
names(gpp_part_sum) = c("lower16","middle","upper84")

# export
# gpp_part_sum %>% 
#   bind_rows %>% 
#   select(param, quant, cv, ela, cont100) %>%
#   mutate(cv = round(cv, 2),
#          ela = round(ela, 2),
#          cont100 = round(cont100, 2)) %>%
#   write_csv("analyses/full_analysis/figures/gpp_part.csv")





#==========
#========== ER (pool years)
#==========

# the functions are designed for separate with vs. between day variation in par_int/temp
# to collapose these into a single variable, set par_int_0 = par_int and par_int_z = 0, etc.
comb_data_trans = comb_data %>% 
  filter(chain==1) %>%
  group_by(chain, step, year, yday) %>%
  mutate(par_int_0 = par_int,
         par_int_z = par_int - par_int_0,
         temp_0 = temp,
         temp_z = temp - temp_0) %>%
  select(-par_int, -temp) %>%
  ungroup

# partition for each step 
er_part = comb_data_trans$step %>% 
  unique %>%
  lapply(function(x){
    d = comb_data_trans %>% filter(step==x)
    er_part_fun(years, d)
  }) %>%
  bind_rows()

# partition for each step 
er_part_sum = lapply(c(0.16,0.5,0.84), function(x){
  er_part %>%
    group_by(param) %>%
    summarize_all(funs(quantile), probs = x) %>%
    mutate(quant = x)
})
names(er_part_sum) = c("lower16","middle","upper84")

# export
# er_part_sum %>% 
#   bind_rows %>% 
#   select(param, quant, cv, ela, cont100) %>%
#   mutate(cv = round(cv, 2),
#          ela = round(ela, 2),
#          cont100 = round(cont100, 2)) %>%
#   write_csv("analyses/full_analysis/figures/er_part.csv")






#==========
#========== NEP (pool years)
#==========

# the functions are designed for separate with vs. between day variation in par_int/temp
# to collapose these into a single variable, set par_int_0 = par_int and par_int_z = 0, etc.
comb_data_trans = comb_data %>% 
  filter(chain==1) %>%
  group_by(chain, step, year, yday) %>%
  mutate(par_int_0 = par_int,
         par_int_z = par_int - par_int_0,
         temp_0 = temp,
         temp_z = temp - temp_0) %>%
  select(-par_int, -temp) %>%
  ungroup

# partition for each step 
nep_part = comb_data_trans$step %>% 
  unique %>%
  lapply(function(x){
    d = comb_data_trans %>% filter(step==x)
    nep_part_fun(years, d)
  }) %>%
  bind_rows()

# partition for each step 
nep_part_sum = lapply(c(0.16,0.5,0.84), function(x){
  nep_part %>%
    group_by(param) %>%
    summarize_all(funs(quantile), probs = x) %>%
    mutate(quant = x)
})
names(nep_part_sum) = c("lower16","middle","upper84")

# export
# nep_part_sum %>% 
#   bind_rows %>% 
#   select(param, quant, cv, ela, cont100) %>%
#   mutate(cv = round(cv, 2),
#          ela = round(ela, 2),
#          cont100 = round(cont100, 2)) %>%
#   write_csv("analyses/full_analysis/figures/nep_part.csv")



