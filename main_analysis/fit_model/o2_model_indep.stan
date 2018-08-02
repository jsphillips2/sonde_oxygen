data {
  // declare variables
  // indices
  int N; // number of observations
  int Y; // number of years
  int D; // number of days
  int T_S; // number of time series
  int D_M[N]; // mapping of observations to days
  int K[Y+1]; // number of days in each year
  int S[T_S+1]; // number of steps in each time series
  int o2_st[T_S+1]; // starting positions for each time series
  int dy_st[Y+1]; // starting positions for each day
  // actual data
  vector<lower=0>[N] o2_obs; // observed oxygen [g m^-3]
  vector<lower=0>[N] o2_eq; // equilibrium oxygen [g m^-3] 
  vector<lower=0>[N] light; // light [umol-photons m^-2 s^-1]
  vector<lower=0>[N] temp; // temperature [C]
  vector<lower=0>[N] wspeed; // wind speeed [m s^-1]
  vector<lower=0>[N] sch_conv; // Schmidt number conversion
  real<lower=0> z; // mixing depth [m]
  real<lower=0> temp_ref; // reference temperature [C]
  real<lower=0> k2; // gas exhange constant 2
  real<lower=0> sig_obs; // observation error sd [g m^-3]
}
parameters{
  real<lower=1> gamma_1; // scaling of gpp with temperature
  real<lower=1> gamma_2; // scaling of er with temperature
  real<lower=0> k0; // gas exchange constant 0
  real<lower=0> k1; // gas exchange constant 1
  real<lower=0> sig_proc; // sd of oxygen state process error
  vector[N] o2; // inferred oxygen state [g m^-3]
  vector[D] log_beta0; // mean value for log_beta0
  vector[D] log_alpha; // mean value for log_beta0
  vector[D] log_rho; // mean value for log_rho
}
transformed parameters {
  // declare variables
  vector<lower=0>[D] beta0; // max gpp at temp_ref [g m^-2 h^-1]
  vector<lower=0>[D] alpha; // max gpp at temp_ref [g m^-2 h^-1]
  vector<lower=0>[D] rho; // er at temp_ref [g m^-2 h^-1]
  vector<lower=0>[N] beta; // max gpp at high light [g m^-2 h^-1]
  vector<lower=0>[N] gpp; // gpp [g m^-2 h^-1]
  vector<lower=0>[N] er; // er [g m^-2 h^-1]
  vector[N] nep; // nep [g m^-2 h^-1]
  vector[N] air; // oxygen exchange with atmosphere [g m^-2 h^-1]
  vector[N] o2_pred; // predicted oxygen [g m^-3]
  // daily parameters
    // random walk
    for (d in 1:D){
      beta0[d] = exp(log_beta0[d]);
      alpha[d] = exp(log_alpha[d]);
      rho[d] = exp(log_rho[d]); 
    }
  // predicted oxygen 
  for (n in 1:N) {
    beta[n] = beta0[D_M[n]]*gamma_1^(temp[n] - temp_ref);
    gpp[n] = beta[n]*tanh((alpha[D_M[n]]/beta[n])*light[n]);
    er[n] = rho[D_M[n]]*gamma_2^(temp[n] - temp_ref);
    nep[n] = gpp[n] - er[n];
    air[n] = ((k0 + k1*wspeed[n]^k2)/100)*sch_conv[n]*(o2_eq[n] - o2[n]);
    o2_pred[n] = o2[n] + (nep[n] + air[n])/z;
  }
}
model {
  // priors
  gamma_1 ~ normal(1.1, 0.4) T[1, ]; 
  gamma_2 ~ normal(1.1, 0.4) T[1, ]; 
  k0 ~ normal(2, 2) T[0, ]; 
  k1 ~ normal(0.2, 0.2) T[0, ]; 
  sig_proc ~ normal(100, 100) T[0, ]; 
  // mean values
  log_beta0 ~ normal(5.5, 0.6); 
  log_alpha ~ normal(1, 0.1); 
  log_rho ~ normal(5.5, 0.6); 
  // z values for non-centered parameterization
  // z_beta0 ~ normal(0, 1);
  // z_alpha ~ normal(0, 1);
  // z_rho ~ normal(0, 1);
  // state process
  for (t in 1:T_S) {
    o2[o2_st[t]] ~ normal(o2_obs[o2_st[t]], sig_obs);
    o2[(o2_st[t]+1):(o2_st[t]+S[t]-1)] 
      ~ normal(o2_pred[o2_st[t]:(o2_st[t]+S[t]-2)], sig_proc);
  }
  // observation process
  o2_obs ~ normal(o2, sig_obs);
}
generated quantities{
  // declare variables
  vector[N] error_proc_real;
  vector[N] error_proc_sim;
  vector[N] error_obs_real;
  vector[N] error_obs_sim;
  vector[N] sq_error_proc_real;
  vector[N] sq_error_proc_sim;
  vector[N] sq_error_obs_real;
  vector[N] sq_error_obs_sim;
  real chi_proc_real;
  real chi_proc_sim;
  real chi_obs_real;
  real chi_obs_sim;
  int pos;
  vector[D] GPP; // total daily flux
  vector[D] ER; // total daily flux
  vector[D] NEP; // total daily flux
  vector[D] AIR; // total daily flux
  vector[D] Flux; // total daily flux
  real GPP_mean; // overall mean flux
  real ER_mean; // overall mean flux
  real NEP_mean; // overall mean flux
  for (t in 1:T_S){
    // set initial values to 0 (they don't make sense to calculate)
    error_proc_real[o2_st[t]] = 0;
    error_proc_sim[o2_st[t]] = 0;
    error_obs_real[o2_st[t]] = 0;
    error_obs_sim[o2_st[t]] = 0;
    sq_error_proc_real[o2_st[t]] = 0;
    sq_error_proc_sim[o2_st[t]] = 0;
    sq_error_obs_real[o2_st[t]] = 0;
    sq_error_obs_sim[o2_st[t]] = 0;
    for (n in (o2_st[t]+1):(o2_st[t]+S[t]-1)){
      // errors 
      error_proc_real[n] = o2[n] - o2_pred[n-1];
      error_proc_sim[n] = normal_rng(o2_pred[n-1], sig_proc) - o2_pred[n-1];
      error_obs_real[n] = o2_obs[n] - o2[n];
      error_obs_sim[n] = normal_rng(o2[n], sig_obs) - o2[n];
      // squared errors
      sq_error_proc_real[n] = error_proc_real[n]^2;
      sq_error_proc_sim[n] = error_proc_sim[n]^2;
      sq_error_obs_real[n] = error_obs_real[n]^2;
      sq_error_obs_sim[n] = error_obs_sim[n]^2;
    }
  }
  // chi-squared goodness-of-fit
  chi_proc_real = sum(sq_error_proc_real/(sig_proc^2));
  chi_proc_sim = sum(sq_error_proc_sim/(sig_proc^2));
  chi_obs_real = sum(sq_error_obs_real/(sig_obs^2));
  chi_obs_sim = sum(sq_error_obs_sim/(sig_obs^2));
  // daily fluxes
  pos = 1;
  for (d in 1:D){
    GPP[d] = sum(gpp[pos:(pos+23)]);
    ER[d] = sum(er[pos:(pos+23)]);
    NEP[d] = sum(nep[pos:(pos+23)]);
    AIR[d] = sum(air[pos:(pos+23)]);
    Flux[d] = (NEP[d] + AIR[d])/z;
    pos = pos + 24; // advance position counter
  }
  // mean daily fluxes
  GPP_mean = mean(GPP);
  ER_mean = mean(ER);
  NEP_mean = mean(NEP);
}
