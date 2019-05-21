data {
  // declare variables
  // indices
  int n_obs; // number of observations
  int n_years; // number of years
  int n_days; // number of days
  int n_series; // number of time series
  int map_days[n_obs]; // mapping of observations to days
  int days_per_year[n_years]; // number of days in each year
  int obs_per_series[n_series]; // number of obervsations in each time series
  int obs_per_day[n_days]; // number of observations in each day
  // actual data
  real<lower=0> o2_obs[n_obs]; // observed oxygen [mg m^-3]
  real<lower=0> o2_eq[n_obs]; // equilibrium oxygen [mg m^-3] 
  real<lower=0> light[n_obs]; // light [umol-photons m^-2 s^-1]
  real<lower=0> temp[n_obs]; // temperature [C]
  real<lower=0> wspeed[n_obs]; // wind speeed [m s^-1]
  real<lower=0> sch_conv[n_obs]; // Schmidt number conversion
  real<lower=0> z; // mixing depth [m]
  real<lower=0> temp_ref; // reference temperature [C]
  real<lower=0> k0; // gas exchange constant 0
  real<lower=0> k1; // gas exchange constant 1
  real<lower=0> k2; // gas exhange constant 2
  real<lower=0> sig_obs; // observation error sd [g m^-3]
}
transformed data {
  // declare variables
  real<lower=0> mu; // mean of DO for entire time series
  real<lower=0> tau; // sd of DO for entire time series
  real x_obs[n_obs]; // scaled observed DO
  real x_eq[n_obs]; // scaled equilibrium DO
  real<lower=0> eta; // mean of light for entire time series
  real<lower=0> lambda[n_obs]; // scaled light 
  real<lower=0> k[n_obs]; // scaled gas exchange constant
  // scale data
  mu = mean(o2_obs);
  tau = sd(o2_obs);
  eta = mean(light);
  for (n in 1:n_obs){
    x_obs[n] = (o2_obs[n] - mu)/tau;
    x_eq[n] = (o2_eq[n] - mu)/tau;
    lambda[n] = light[n]/eta;
    k[n] = (1/z)*((k0 + k1*wspeed[n]^k2)/100)*sch_conv[n];
  }
}
parameters {
  real<lower=1> gamma_1; // scaling of gpp with temperature
  real<lower=1> gamma_2; // scaling of er with temperature
  real<lower=0> sig_proc; // sd of oxygen state process error
  real x[n_obs]; // scaled inferred oxygen state
  real log_b0; // initial value for log_b0
  real log_a; // initial value for log_b0
  real log_r; // initial value for log_r
  real z_b0[n_days-n_years]; // z value for non-centered parameterization
  real z_a[n_days-n_years]; // z value for non-centered parameterization
  real z_r[n_days-n_years]; // z value for non-centered parameterization
}
transformed parameters {
  // declare variables
  real b0; // scaled max gpp at temp_ref
  real a; // scaled initial slope
  real r; // scaled er at temp_ref
  real b[n_obs]; // scaled max gpp at high light
  real chi[n_obs]; // scaled gpp 
  real kappa[n_obs]; // scaled er 
  real phi[n_obs]; // nep [g m^-2 h^-1]
  real nu[n_obs]; // oxygen exchange with atmosphere [g m^-2 h^-1]
  real x_pred[n_obs]; // predicted oxygen [g m^-3]
  // exponentiate parameters
  b0 = exp(log_b0); 
  a = exp(log_a);
  r = exp(log_r);
  // predicted oxygen 
  for (n in 1:n_obs) {
    b[n] = b0*gamma_1^(temp[n] - temp_ref);
    chi[n] = b[n]*tanh((a/b[n])*lambda[n]);
    kappa[n] = r*gamma_2^(temp[n] - temp_ref);
    phi[n] = chi[n] - kappa[n];
    nu[n] = k[n]*(x_eq[n] - x[n]);
    x_pred[n] = x[n] + phi[n] + nu[n];
  }
}
model {
 // priors
  gamma_1 ~ normal(1, 1) T[1, ]; 
  gamma_2 ~ normal(1, 1) T[1, ]; 
  sig_proc ~ normal(0, 1) T[0, ]; 
  // initial values
  log_b0 ~ normal(0, 1); 
  log_a ~ normal(0, 1); 
  log_r ~ normal(0, 1); 
  // z values for non-centered parameterization
  z_b0 ~ normal(0, 1);
  z_a ~ normal(0, 1);
  z_r ~ normal(0, 1);
  // state process
  {
    int pos = 1; 
    for (t in 1:n_series) {
      x[pos] ~ normal(x_obs[pos], sig_obs);
      x[(pos+1):(pos+obs_per_series[t]-1)] 
        ~ normal(x_pred[pos:(pos+obs_per_series[t]-2)], sig_proc);
      pos = pos + obs_per_series[t];
    }
  }
  // observation process
  x_obs ~ normal(x, sig_obs);
}
generated quantities {
  real o2[n_obs]; // inferred oxygen state [mg m^-3]
  real o2_pred[n_obs]; // predicted oxygen [mg m^-3]
  real beta0[n_days]; // max gpp at temp_ref [mg m^-2 h^-1]
  real alpha[n_days]; // initial slope [mg-O2 s umol-photons-1 m-1 h-1]
  real rho[n_days]; // er at temp_ref [mg m^-2 h^-1]
  real beta[n_obs]; // max gpp at high light [mg m^-2 h^-1]
  real gpp[n_obs]; // gpp [mg m^-2 h^-1]
  real er[n_obs]; // er [mg m^-2 h^-1]
  real nep[n_obs]; // nep [mg m^-2 h^-1]
  real air[n_obs]; // oxygen exchange with atmosphere [mg m^-2 h^-1]
  real GPP[n_days]; // total daily flux [g m^-2 d^-1]
  real ER[n_days]; // total daily flux [g m^-2 d^-1]
  real NEP[n_days]; // total daily flux [g m^-2 d^-1]
  real AIR[n_days]; // total daily flux [g m^-2 d^-1]
  real Flux[n_days]; // total daily flux [g m^-2 d^-1]
  real GPP_mean; // overall mean flux [g m^-2 d^-1]
  real ER_mean; // overall mean flux [g m^-2 d^-1]
  real NEP_mean; // overall mean flux [g m^-2 d^-1]
  real log_lik[n_obs]; // log-likelihood
  // back-tranformed scaled variables [g m^-2 d^-1]
  for(n in 1:n_obs){
    o2[n] = tau*x[n] + mu;
    o2_pred[n] = tau*x_pred[n] + mu;
    beta[n] = tau*z*b[n];
    gpp[n] = tau*z*chi[n];
    er[n] = tau*z*kappa[n];
    nep[n] = tau*z*phi[n];
    air[n] = tau*z*nu[n];
    log_lik[n] = normal_lpdf(x_obs[n]|x[n], sig_obs);
  }
  {
    int pos = 1;
    for (d in 1:n_days){
      beta0[d] = tau*z*b0;
      alpha[d] = tau*z*a/eta;
      rho[d] = tau*z*r;
      GPP[d] = 24*mean(gpp[pos:(pos + obs_per_day[d] - 1)]);
      ER[d] = 24*mean(er[pos:(pos + obs_per_day[d] - 1)]);
      NEP[d] = 24*mean(nep[pos:(pos + obs_per_day[d] - 1)]);
      AIR[d] = 24*mean(air[pos:(pos + obs_per_day[d] - 1)]);
      Flux[d] = (NEP[d] + AIR[d])/z;
      pos = pos + obs_per_day[d]; // advance position counter
  }
  // mean daily fluxes
    GPP_mean = mean(GPP);
    ER_mean = mean(ER);
    NEP_mean = mean(NEP);  
  }
}
