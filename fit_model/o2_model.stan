data {
  // declare variables
  // indeces
  int N; // number of observations
  int Y; // number of years
  int D; // number of days
  int T_S; // number of time series
  int D_M[N]; // mapping of observations to days
  int K[Y]; // number of days in each year
  int S[T_S]; // number of steps in each time series
  int o2_st[T_S]; // starting positions for each time series
  int dy_st[Y]; // starting positions for each day
  // actual data
  vector<lower=0>[N] o2_obs; // observed oxygen [g m^-3]
  vector<lower=0>[N] o2_eq; // equilibrium oxygen [g m^-3] 
  vector<lower=0>[N] light; // light [umol-photons m^-2 s^-1]
  vector<lower=0>[N] temp; // temperature [C]
  vector<lower=0>[N] wspeed; // wind speeed [m s^-1]
  vector<lower=0>[N] sch_conv; // Schmidt number conversion
  real<lower=0> z; // mixing depth [m]
  real<lower=0> temp_ref; // reference temperature [C]
  real<lower=0> k0; // gas exchange constant 1
  real<lower=0> k1; // gas exchange constant 2
  real<lower=0> k3; // gas exhange constant 3
  real<lower=0> sig_obs; // observation error sd [g m^-3]
}
parameters{
  real<lower=0> alpha; // slope of gpp ~ light at low light
  real<lower=1> gamma_1; // scaling of gpp with temperature
  real<lower=1> gamma_2; // scaling of er with temperature
  real<lower=0> sig_beta0; // sd of log_beta0 random walk 
  real<lower=0> sig_rho; // sd of log_rho random walk 
  real<lower=0> sig_proc; // sd of oxygen state process error
  vector[D] log_beta0; // max gpp at temp_ref (log scale)
  vector[D] log_rho; // er at temp_ref (log scale)
  vector[N] o2; // inferred oxygen state [g m^-3]
}
transformed parameters {
  // declare variables
  vector<lower=0>[D] beta0; // max gpp at temp_ref [g m^-2 h^-1]
  vector<lower=0>[D] rho; // er at temp_ref [g m^-2 h^-1]
  vector<lower=0>[N] beta; // max gpp at high light [g m^-2 h^-1]
  vector<lower=0>[N] gpp; // gpp [g m^-2 h^-1]
  vector<lower=0>[N] er; // er [g m^-2 h^-1]
  vector[N] nep; // nep [g m^-2 h^-1]
  vector[N] air; // oxygen exchange with atmosphere [g m^-2 h^-1]
  vector[N] o2_pred; // predicted oxygen [g m^-3]
  // exp parameters
  beta0 = exp(log_beta0); 
  rho = exp(log_rho);
  // predicted oxygen 
  for (n in 1:N) {
    beta[n] = beta0[D_M[n]]*gamma_1^(temp[n] - temp_ref);
    gpp[n] = beta[n]*tanh((alpha/beta[n])*light[n]);
    er[n] = rho[D_M[n]]*gamma_2^(temp[n] - temp_ref);
    nep[n] = gpp[n] - er[n];
    air[n] = ((k0 + k1*wspeed[n]^k3)/100)*sch_conv[n]*(o2_eq[n] - o2[n]);
    o2_pred[n] = o2[n] + (nep[n] + air[n])/z;
  }
}
model {
  // priors
  alpha ~ normal(3, 1.5) T[0, ]; 
  gamma_1 ~ normal(1.1, 0.4) T[1, ]; 
  gamma_2 ~ normal(1.1, 0.4) T[1, ]; 
  sig_beta0 ~ normal(0.5, 0.6) T[0, ]; 
  sig_rho ~ normal(0.5, 0.6) T[0, ]; 
  sig_proc ~ normal(100, 100) T[0, ]; 
  // random walk for daily parameters
  for (y in 1:Y) {
    log_beta0[dy_st[y]] ~ normal(6, 0.6); 
    log_rho[dy_st[y]] ~ normal(5.5, 0.6); 
    log_beta0[(dy_st[y]+1):(dy_st[y]+K[y]-1)] 
      ~ normal(log_beta0[dy_st[y]:(dy_st[y]+K[y]-2)], sig_beta0);
    log_rho[(dy_st[y]+1):(dy_st[y]+K[y]-1)] 
      ~ normal(log_rho[dy_st[y]:(dy_st[y]+K[y]-2)], sig_rho);
  }
  // state process
  for (t in 1:T_S) {
    o2[o2_st[t]] ~ normal(o2_obs[o2_st[t]], sig_obs);
    o2[(o2_st[t]+1):(o2_st[t]+S[t]-1)] 
      ~ normal(o2_pred[o2_st[t]:(o2_st[t]+S[t]-2)], sig_proc);
  }
  // observation process
  o2_obs ~ normal(o2, sig_obs);
}
