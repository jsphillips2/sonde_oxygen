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
  // values from previous fit
  // vector<lower=0>[D] beta0_f; // max gpp at temp_ref [g m^-2 h^-1]
  vector<lower=0>[N] beta_f; // max gpp at temp_ref [g m^-2 h^-1]
  vector<lower=0>[D] alpha_f; // max gpp at temp_ref [g m^-2 h^-1]
  vector<lower=0>[D] rho_f; // er at temp_ref [g m^-2 h^-1]
  real<lower=1> gamma_1_f; // scaling of gpp with temperature
  real<lower=1> gamma_2_f; // scaling of er with temperature
  real<lower=0> k0_f; // gas exchange constant 0
  real<lower=0> k1_f; // gas exchange constant 1
  real<lower=0> sig_proc_f; // sd of oxygen state process error
}
transformed data{
  // vector<lower=0>[N] beta_s; // max gpp at high light [g m^-2 h^-1]
  vector<lower=0>[N] gpp_s; // simulated gpp
  vector<lower=0>[N] er_s; // simulated er 
  vector[N] nep_s; // simulated nep 
  vector[N] air_s; // simulated exchange with atmosphere
  vector[N] o2_pred_s; // simulated predicted oxygen  
  vector[N+1] o2_s; // simulated "real" oxygen  
  vector[N] o2_obs_s; // simulated observed oxygen  
  // simulate state process
  for (t in 1:T_S) {
    // initial oxygen values
    o2_s[o2_st[t]] = o2_obs[o2_st[t]];
    for(n in (o2_st[t]+1):(o2_st[t]+S[t])){
      // beta_s[n-1] = beta0_f[D_M[n-1]]*gamma_1_f^(temp[n-1] - temp_ref);
      //gpp_s[n-1] = beta_s[n-1]*tanh((alpha_f[D_M[n-1]]/beta_s[n-1])*light[n-1]);
      gpp_s[n-1] = beta_f[n-1]*tanh((alpha_f[D_M[n-1]]/beta_f[n-1])*light[n-1]);
      er_s[n-1] = rho_f[D_M[n-1]]*gamma_2_f^(temp[n-1] - temp_ref);
      nep_s[n-1] = gpp_s[n-1] - er_s[n-1];
      air_s[n-1] = ((k0_f + k1_f*wspeed[n-1]^k2)/100)*sch_conv[n-1]*(o2_eq[n-1] 
        - o2_s[n-1]);
      o2_pred_s[n-1] = o2_s[n-1] + (nep_s[n-1] + air_s[n-1])/z;      
      o2_s[n] = normal_rng(o2_pred_s[n-1], sig_proc_f);
    }
  }
  // simulate observation process
  for (n in 1:N){
    o2_obs_s[n] = normal_rng(o2_s[n], sig_obs);
  }
}
parameters{
  real<lower=1> gamma_1; // scaling of gpp with temperature
  real<lower=1> gamma_2; // scaling of er with temperature
  real<lower=0> k0; // gas exchange constant 0
  real<lower=0> k1; // gas exchange constant 1
  real<lower=0> sig_beta0; // sd of log_beta0 random walk 
  real<lower=0> sig_alpha; // sd of log_rho random walk 
  real<lower=0> sig_rho; // sd of log_rho random walk 
  real<lower=0> sig_proc; // sd of oxygen state process error
  vector[N] o2; // inferred oxygen state [g m^-3]
  vector[Y+1] log_beta0_init; // initial value for log_beta0
  vector[Y+1] log_alpha_init; // initial value for log_beta0
  vector[Y+1] log_rho_init; // initial value for log_rho
  vector[D] z_beta0; // z value for non-centered parameterization
  vector[D] z_alpha; // z value for non-centered parameterization
  vector[D] z_rho; // z value for non-centered parameterization
}
transformed parameters {
  // declare variables
  vector[D] log_beta0; // max gpp at temp_ref (log scale)
  vector[D] log_alpha; // max gpp at temp_ref (log scale)
  vector[D] log_rho; // er at temp_ref (log scale)
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
  for (y in 1:Y) {
    // inital value
    log_beta0[dy_st[y]] = log_beta0_init[y];
    log_alpha[dy_st[y]] = log_alpha_init[y];
    log_rho[dy_st[y]] = log_rho_init[y];
    // random walk
    for (d in (dy_st[y]+1):(dy_st[y]+K[y]-1)){
      log_beta0[d] = log_beta0[d-1] + sig_beta0*z_beta0[d-1];
      log_alpha[d] = log_alpha[d-1] + sig_alpha*z_alpha[d-1];
      log_rho[d] = log_rho[d-1] + sig_rho*z_rho[d-1]; 
    }
  }
  // exp parameters
  beta0 = exp(log_beta0); 
  alpha = exp(log_alpha);
  rho = exp(log_rho);
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
  sig_beta0 ~ normal(0.5, 0.6) T[0, ]; 
  sig_alpha ~ normal(0.5, 0.6) T[0, ]; 
  sig_rho ~ normal(0.5, 0.6) T[0, ]; 
  sig_proc ~ normal(100, 100) T[0, ]; 
  // initial values
  log_beta0_init ~ normal(6, 0.6); 
  log_alpha_init ~ normal(1, 0.1); 
  log_rho_init ~ normal(5.5, 0.6); 
  // z values for non-centered parameterization
  z_beta0 ~ normal(0, 1);
  z_alpha ~ normal(0, 1);
  z_rho ~ normal(0, 1);
  // state process
  for (t in 1:T_S) {
    o2[o2_st[t]] ~ normal(o2_obs_s[o2_st[t]], sig_obs);
    o2[(o2_st[t]+1):(o2_st[t]+S[t]-1)] 
      ~ normal(o2_pred[o2_st[t]:(o2_st[t]+S[t]-2)], sig_proc);
  }
  // observation process
  o2_obs_s ~ normal(o2, sig_obs);
}
generated quantities{
  // declare variables
  int pos;
  vector[D] GPP; // total daily flux
  vector[D] ER; // total daily flux
  vector[D] NEP; // total daily flux
  vector[D] AIR; // total daily flux
  vector[D] Flux; // total daily flux
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
}
