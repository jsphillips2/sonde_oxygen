data {
  // declare variables
  // indices
  int n_obs; // number of observations
  int n_years; // number of years
  int n_days; // number of days
  int n_series; // number of time series
  int map_days[n_obs]; // mapping of observations to days
  int days_per[n_years]; // number of days in each year
  int obs_per[n_series]; // number of steps in each time series
  // actual data
  real<lower=0> o2_obs[n_obs]; // observed oxygen [g m^-3]
  real<lower=0> o2_eq[n_obs]; // equilibrium oxygen [g m^-3] 
  real<lower=0> light[n_obs]; // light [umol-photons m^-2 s^-1]
  real<lower=0> temp[n_obs]; // temperature [C]
  real<lower=0> wspeed[n_obs]; // wind speeed [m s^-1]
  real<lower=0> sch_conv[n_obs]; // Schmidt number conversion
  real<lower=0> z; // mixing depth [m]
  real<lower=0> temp_ref; // reference temperature [C]
  real<lower=0> k2; // gas exhange constant 2
  real<lower=0> sig_obs; // observation error sd [g m^-3]
  // priors
  real k0_prior[2];
  real k1_prior[2];
  real gamma_1_prior[2];
  real gamma_2_prior[2];
  real sig_beta0_prior[2];
  real sig_alpha_prior[2];
  real sig_rho_prior[2];
  real sig_proc_prior[2];
  real log_beta0_prior[2];
  real log_alpha_prior[2];
  real log_rho_prior[2];
}
parameters{
  real<lower=0> k0; // gas exchange constant 0
  real<lower=0> k1; // gas exchange constant 1
  real<lower=1> gamma_1; // scaling of gpp with temperature
  real<lower=1> gamma_2; // scaling of er with temperature
  real<lower=0> sig_beta0; // sd of log_beta0 random walk 
  real<lower=0> sig_alpha; // sd of log_rho random walk 
  real<lower=0> sig_rho; // sd of log_rho random walk 
  real<lower=0> sig_proc; // sd of oxygen state process error
  real o2[n_obs]; // inferred oxygen state [g m^-3]
  real log_beta0_init[n_years]; // initial value for log_beta0
  real log_alpha_init[n_years]; // initial value for log_beta0
  real log_rho_init[n_years]; // initial value for log_rho
  real z_beta0[n_days-n_years]; // z value for non-centered parameterization
  real z_alpha[n_days-n_years]; // z value for non-centered parameterization
  real z_rho[n_days-n_years]; // z value for non-centered parameterization
}
transformed parameters {
  // declare variables
  real log_beta0[n_days]; // max gpp at temp_ref (log scale)
  real log_alpha[n_days]; // max gpp at temp_ref (log scale)
  real log_rho[n_days]; // er at temp_ref (log scale)
  real beta0[n_days]; // max gpp at temp_ref [g m^-2 h^-1]
  real alpha[n_days]; // max gpp at temp_ref [g m^-2 h^-1]
  real rho[n_days]; // er at temp_ref [g m^-2 h^-1]
  real beta[n_obs]; // max gpp at high light [g m^-2 h^-1]
  real gpp[n_obs]; // gpp [g m^-2 h^-1]
  real er[n_obs]; // er [g m^-2 h^-1]
  real nep[n_obs]; // nep [g m^-2 h^-1]
  real air[n_obs]; // oxygen exchange with atmosphere [g m^-2 h^-1]
  real o2_pred[n_obs]; // predicted oxygen [g m^-3]
  // daily parameters
  {
    int pos = 1;
    for (y in 1:n_years) {
      // inital value
      log_beta0[pos] = log_beta0_init[y];
      log_alpha[pos] = log_alpha_init[y];
      log_rho[pos] = log_rho_init[y];
      // random walk
      for (d in (pos+1):(pos+days_per[y]-1)){
        log_beta0[d] = log_beta0[d-1] + sig_beta0*z_beta0[d - y];
        log_alpha[d] = log_alpha[d-1] + sig_alpha*z_alpha[d - y];
        log_rho[d] = log_rho[d-1] + sig_rho*z_rho[d - y]; 
      }
      pos = pos + days_per[y];
    }
  }
  // exp parameters
  beta0 = exp(log_beta0); 
  alpha = exp(log_alpha);
  rho = exp(log_rho);
  // predicted oxygen 
  for (n in 1:n_obs) {
    beta[n] = beta0[map_days[n]]*gamma_1^(temp[n] - temp_ref);
    gpp[n] = beta[n]*tanh((alpha[map_days[n]]/beta[n])*light[n]);
    er[n] = rho[map_days[n]]*gamma_2^(temp[n] - temp_ref);
    nep[n] = gpp[n] - er[n];
    air[n] = ((k0 + k1*wspeed[n]^k2)/100)*sch_conv[n]*(o2_eq[n] - o2[n]);
    o2_pred[n] = o2[n] + (nep[n] + air[n])/z;
  }
}
model {
  // priors
  k0 ~ normal(k0_prior[1], k0_prior[2]) T[0, ]; 
  k1 ~ normal(k1_prior[1], k1_prior[2]) T[0, ]; 
  gamma_1 ~ normal(gamma_1_prior[1], gamma_1_prior[2]) T[1, ]; 
  gamma_2 ~ normal(gamma_2_prior[1], gamma_2_prior[2]) T[1, ]; 
  sig_beta0 ~ normal(sig_beta0_prior[1], sig_beta0_prior[2]) T[0, ]; 
  sig_alpha ~ normal(sig_alpha_prior[1], sig_alpha_prior[2]) T[0, ]; 
  sig_rho ~ normal(sig_rho_prior[1], sig_rho_prior[2]) T[0, ]; 
  sig_proc ~ normal(sig_proc_prior[1], sig_proc_prior[2]) T[0, ]; 
  // initial values
  for(y in 1:n_years){
    log_beta0_init[y] ~ normal(log_beta0_prior[1], log_beta0_prior[2]); 
    log_alpha_init[y] ~ normal(log_alpha_prior[1], log_alpha_prior[2]); 
    log_rho_init[y] ~ normal(log_rho_prior[1], log_rho_prior[2]); 
  }
  // z values for non-centered parameterization
  z_beta0 ~ normal(0, 1);
  z_alpha ~ normal(0, 1);
  z_rho ~ normal(0, 1);
  // state process
  {
    int pos = 1; 
    for (t in 1:n_series) {
      o2[pos] ~ normal(o2_obs[pos], sig_obs);
      o2[(pos+1):(pos+obs_per[t]-1)] 
        ~ normal(o2_pred[pos:(pos+obs_per[t]-2)], sig_proc);
      pos = pos + obs_per[t];
    }
  }
  // observation process
  o2_obs ~ normal(o2, sig_obs);
}
generated quantities{
  real[n_days] GPP; // total daily flux
  real[n_days] ER; // total daily flux
  real[n_days] NEP; // total daily flux
  real[n_days] AIR; // total daily flux
  real[n_days] Flux; // total daily flux
  real GPP_mean; // overall mean flux
  real ER_mean; // overall mean flux
  real NEP_mean; // overall mean flux
  {
    int pos = 1;
    for (d in 1:n_days){
      int pos_max = sum(map_days == d);
      GPP[d] = 24*mean(gpp[pos:pos_max]);
      ER[d] = 24*mean(er[pos:pos_max]);
      NEP[d] = 24*mean(nep[pos:pos_max]);
      AIR[d] = 24*mean(air[pos:pos_max]);
      Flux[d] = (NEP[d] + AIR[d])/z;
      pos = pos + pos_max; // advance position counter
  }
  // mean daily fluxes
    GPP_mean = mean(GPP);
    ER_mean = mean(ER);
    NEP_mean = mean(NEP);  
  }
}
