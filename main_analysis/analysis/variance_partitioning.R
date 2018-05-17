#==========
#========== Preliminaries
#==========

# load packages
library(tidyverse)
library(matrixcalc)

# import data and model fit
sonde_data = read_csv("data/sonde_prep.csv")
params_full = read_csv("main_analysis/model_output/fixed_pars_full.csv")
daily = read_csv("main_analysis/model_output/daily_full.csv")

# base theme
theme_base = theme_bw()+
  theme(panel.grid=element_blank(),
        strip.background=element_blank(),
        text = element_text(size=12),
        strip.text = element_text(size=10),
        axis.text=element_text(size=10, color="black"),
        axis.title.y=element_text(margin=margin(0,15,0,0)),
        axis.title.x=element_text(margin=margin(15,0,0,0)))





#==========
#========== Prepare Data
#==========

# combine data
comb_data = daily %>%
  spread(name, value) %>%
  full_join(params_full) %>%
  left_join(sonde_data %>%
              rename(day = D_M) %>%
              select(year,yday,day,hour,par,temp) %>%
              na.omit()) %>%
  select(chain,step,year,yday,hour,par,temp,beta0,alpha,rho,alpha,gamma_1,gamma_2) %>%
  mutate(temp_ref = mean(temp))





#==========
#========== Define functions
#==========

# GPP
gpp_part_fun = function(y, x){
  
  # Select year from data frame of inputs
  x_filt = x %>% filter(year %in% y) 
  
  # Vector of mean values
  gpp_par_mean = x_filt %>%
    select(beta0, alpha, gamma_1, par_0, par_z, temp_0, temp_z, temp_ref) %>%
    summarise_all(mean)
  
  # Define model function
  gpp_fun = function(beta0,alpha,gamma_1,par_0,par_z,temp_0,temp_z,temp_ref){
    beta0*gamma_1^((temp_0 + temp_z)-temp_ref)*tanh((alpha/(beta0*gamma_1^((temp_0 + temp_z)-temp_ref)))*(par_0 + par_z))
  }
  
  # Evaluate at mean value
  gpp_mean = gpp_par_mean %>% 
  {gpp_fun(.$beta0,.$alpha,.$gamma_1,.$par_0,.$par_z,.$temp_0,.$temp_z,.$temp_ref)}
  
  # Define model gradient
  gpp_grad_fun <-deriv(~beta0*gamma_1^((temp_0 + temp_z)-temp_ref)*tanh((alpha/(beta0*gamma_1^((temp_0 + temp_z)-temp_ref)))*(par_0 + par_z)), 
                       c("beta0","alpha","gamma_1","par_0","par_z","temp_0","temp_z","temp_ref"), 
                       function(beta0,alpha,gamma_1,par_0,par_z,temp_0,temp_z,temp_ref){})
  gpp_grad = t(attributes(gpp_par_mean %>% {gpp_grad_fun(.$beta0,.$alpha,.$gamma_1,.$par_0,.$par_z,.$temp_0,.$temp_z,.$temp_ref)})$gradient)
  
  # Define variance-covariance matrix
  gpp_vcv = cov(x_filt %>% select(beta0,alpha,gamma_1,par_0,par_z,temp_0,temp_z,temp_ref)) 
  
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
                  ela = round((mean/gpp_mean)*sen,2),
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
                  ela = round((mean/er_mean)*sen,2),
                  cont100 = er_cont) %>% 
    filter(sd > 0)
  return(yy)
}




#==========
#========== GPP (pool years)
#==========

# the functions are designed for separate with vs. between day variation in par/temp
# to collapose these into a single variable, set par_0 = par and par_z = 0, etc.
comb_data_trans = comb_data %>% 
  filter(chain==1) %>%
  group_by(chain, step, year, yday) %>%
  mutate(par_0 = par,
         par_z = par - par_0,
         temp_0 = temp,
         temp_z = temp - temp_0) %>%
  select(-par, -temp) %>%
  ungroup

# partition for each step 
gpp_part = comb_data_trans$step %>% 
  unique %>%
  lapply(function(x){
    d = comb_data_trans %>% filter(step==x)
    gpp_part_fun(c(2013:2017), d)
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
# write_csv(gpp_part_sum %>% bind_rows, "analysis/tables/gpp_part.csv")





#==========
#========== ER (pool years)
#==========

# the functions are designed for separate with vs. between day variation in par/temp
# to collapose these into a single variable, set par_0 = par and par_z = 0, etc.
comb_data_trans = comb_data %>% 
  filter(chain==1) %>%
  group_by(chain, step, year, yday) %>%
  mutate(par_0 = par,
         par_z = par - par_0,
         temp_0 = temp,
         temp_z = temp - temp_0) %>%
  select(-par, -temp) %>%
  ungroup

# partition for each step 
er_part = comb_data_trans$step %>% 
  unique %>%
  lapply(function(x){
    d = comb_data_trans %>% filter(step==x)
    er_part_fun(c(2013:2017), d)
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
# write_csv(er_part_sum %>% bind_rows, "analysis/tables/er_part.csv")


