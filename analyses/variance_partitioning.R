#==========
#========== Preliminaries
#==========

# load packages
library(tidyverse)
library(matrixcalc)

# import data and model fit
input_dir = "analyses/model_fit/"
sonde_data = read_csv(paste0(input_dir,"input/sonde_prep.csv"))
params_full = read_csv(paste0(input_dir,"output/fixed_pars_full.csv"))
daily = read_csv(paste0(input_dir,"output/daily_full.csv"))





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

# trim data for efficiency
set.seed(1)
comb_data_slim = comb_data %>%
  filter(chain == 1, step %in% sample(1:max(comb_data$step), 500, replace = F)) 

# select years
years = unique(sonde_data$year)





#==========
#========== Define functions
#==========

# GPP gradient and VCOV matrix
gpp_grad = function(x){
  
  # vector of mean values
  par_mean = x %>%
    select(beta0, alpha, gamma_1, par_int, temp, temp_ref) %>%
    summarize_all(mean)
  
  # define model function
  model = function(beta0,alpha,gamma_1,par_int,temp,temp_ref){
    beta0*gamma_1^(temp-temp_ref)*tanh((alpha/(beta0*gamma_1^(temp-temp_ref)))*par_int)
  }
  
  # evaluate at mean value
  model_mean = par_mean %>% 
  {model(.$beta0,.$alpha,.$gamma_1,.$par_int,.$temp,.$temp_ref)}
  
  # define model gradient
  grad_fun <-deriv(~beta0*gamma_1^(temp-temp_ref)*tanh((alpha/(beta0*gamma_1^(temp-temp_ref)))*par_int), 
                   c("beta0","alpha","gamma_1","par_int","temp","temp_ref"), 
                   function(beta0,alpha,gamma_1,par_int,temp,temp_ref){})
  grad = t(attributes(par_mean %>% {grad_fun(.$beta0,.$alpha,.$gamma_1,.$par_int,.$temp,.$temp_ref)})$gradient)
  
  # define variance-covariance matrix
  vcv = cov(x %>% select(beta0,alpha,gamma_1,par_int,temp,temp_ref)) 
  
  # export
  return(list(par_mean = par_mean, model_mean = model_mean, grad = grad, vcv = vcv))
}



# ER gradient and VCOV matrix
er_grad = function(x){
  
  # vector of mean values
  par_mean = x %>%
    select(rho, gamma_2, temp, temp_ref) %>%
    summarise_all(mean)
  
  # define model function
  model = function(rho, gamma_2, temp, temp_ref){
    rho*gamma_2^(temp-temp_ref)
  }
  
  # evaluate at mean value
  model_mean = par_mean %>% 
  {model(.$rho, .$gamma_2, .$temp, .$temp_ref)}
  
  # define model gradient
  grad_fun <-deriv(~rho*gamma_2^(temp-temp_ref), c("rho", "gamma_2", "temp", "temp_ref"), 
                   function(rho, gamma_2, temp, temp_ref){})
  grad = t(attributes(par_mean %>% {grad_fun(.$rho, .$gamma_2, .$temp, .$temp_ref)})$gradient)
  
  # define variance-covariance matrix
  vcv = cov(x %>% select(rho, gamma_2, temp, temp_ref)) 
  
  # export
  return(list(par_mean = par_mean, model_mean = model_mean, grad = grad, vcv = vcv))
} 



# NEP gradient and VCOV matrix
nep_grad = function(x){
  
  # Vector of mean values
  par_mean = x %>%
    select(beta0, alpha, rho, gamma_1, gamma_2, par_int, temp, temp_ref) %>%
    summarise_all(mean)
  
  # Define model function
  model = function(beta0,alpha,rho,gamma_1,gamma_2,par_int,temp,temp_ref){
    beta0*gamma_1^(temp-temp_ref)*tanh((alpha/(beta0*gamma_1^(temp-temp_ref)))*par_int) 
    - rho*gamma_2^(temp-temp_ref)
  }
  
  # Evaluate at mean value
  model_mean = par_mean %>% 
  {model(.$beta0,.$alpha,.$rho,.$gamma_1,.$gamma_2,.$par_int,.$temp,.$temp_ref)}
  
  # Define model gradient
  grad_fun <-deriv(~beta0*gamma_1^(temp-temp_ref)*tanh((alpha/(beta0*gamma_1^(temp-temp_ref)))*par_int) 
                   - rho*gamma_2^(temp-temp_ref), 
                   c("beta0","alpha","rho","gamma_1","gamma_2","par_int","temp","temp_ref"), 
                   function(beta0,alpha,rho,gamma_1,gamma_2,par_int,temp,temp_ref){})
  grad = t(attributes(par_mean %>% 
  {grad_fun(.$beta0,.$alpha,.$rho,.$gamma_1,.$gamma_2,
            .$par_int,.$temp,.$temp_ref)})$gradient)
  
  # Define variance-covariance matrix
  vcv = cov(x %>% select(beta0,alpha,rho,gamma_1,gamma_2,par_int,temp,temp_ref)) 
  
  # export
  return(list(par_mean = par_mean, model_mean = model_mean, grad = grad, vcv = vcv))
}




# partition variance
partition_fun = function(d){
  
  # calculate contribution and relative contribution to variance
  var_cont = hadamard.prod(d$vcv, d$grad %*% t(d$grad))
  rel_cont = rowSums(var_cont)/sum(var_cont)
  
  # define "slim" matrices (keeping variables with non-0 variances)
  keep = which(rowSums(var_cont) != 0)
  cont_slim = var_cont[keep,keep]
  vcv_slim = d$vcv[keep,keep]
  
  # squared coefficient of variation
  cv_2 = diag(vcv_slim)/(d$par_mean[keep]^2)
  
  # squared sensitivity (including covariance terms)
  sen_2 = rowSums(cont_slim)/diag(vcv_slim)*(d$par_mean[keep]^2)
  
  # scale squared sensitivity by total variance
  scale_sen = sen_2/sum(cv_2*sen_2)
  
  return(data.frame(var = names(cv_2), cv_2 = t(cv_2), 
                    scale_sen = t(scale_sen), rel_cont = rel_cont[keep]) %>% tbl_df())
  
}





#==========
#========== GPP
#==========

# partition for each step 
gpp_part = comb_data_slim %>% 
  split(.$step) %>%
  lapply(function(x){partition_fun(gpp_grad(x))}) %>%
  bind_rows()


# partition for each step 
gpp_part_sum = lapply(c(0.16,0.5,0.84), function(x){
  gpp_part %>%
    group_by(var) %>%
    summarize_all(funs(quantile), probs = x) %>%
    mutate(quant = x)
})
names(gpp_part_sum) = c("lower16","middle","upper84")





#==========
#========== ER
#==========

# partition for each step 
er_part = comb_data_slim %>% 
  split(.$step) %>%
  lapply(function(x){partition_fun(er_grad(x))}) %>%
  bind_rows()


# partition for each step 
er_part_sum = lapply(c(0.16,0.5,0.84), function(x){
  er_part %>%
    group_by(var) %>%
    summarize_all(funs(quantile), probs = x) %>%
    mutate(quant = x)
})
names(er_part_sum) = c("lower16","middle","upper84")





#==========
#========== NEP
#==========

# partition for each step 
nep_part = comb_data_slim %>% 
  split(.$step) %>%
  lapply(function(x){partition_fun(nep_grad(x))}) %>%
  bind_rows()


# partition for each step 
nep_part_sum = lapply(c(0.16,0.5,0.84), function(x){
  nep_part %>%
    group_by(var) %>%
    summarize_all(funs(quantile), probs = x) %>%
    mutate(quant = x)
})
names(nep_part_sum) = c("lower16","middle","upper84")



