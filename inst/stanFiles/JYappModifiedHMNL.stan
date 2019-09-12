// original code for HMNL from J Yapp
//https://discourse.mc-stan.org/t/speeding-up-a-hierarchical-multinomial-logit-model/1538/5
// modified  tmn YW references commented out and YB changed to Y
// ALSO Now know that setting coef 1 to 0 is WRONG so sampling Theta[K] directly.
// https://discourse.mc-stan.org/t/speeding-up-a-hierarchical-multinomial-logit-model/1538/9?u=tnearey
data {
  int<lower=2> C; // Number of alternatives (choices) in each scenario
  int<lower=1> K; // Number of alternatives
  int<lower=1> R; // Number of respondents
  int<lower=1> S; // Number of scenarios per respondent
  int<lower=1,upper=C> Y[R, S]; // best choices
  // int<lower=1,upper=C> YW[R, S]; // worst choices
  matrix[C, K] X[R, S]; // matrix of attributes for each obs
}

parameters {
  vector[K] Beta[R];
  vector[K] Theta;
  cholesky_factor_corr[K] L_Omega;
  // vector<lower=0, upper=pi()/2>[K] L_sigma_unif;
  vector <lower=0> [K] L_sigma_unif; //  we're making this simpler
}

transformed parameters {
  vector<lower=0>[K] L_sigma;
  matrix[K, K] L_Sigma;
  vector[C] XB[R, S];
  // vector[K] Theta;

  for (k in 1:K) {
    // L_sigma[k] = 2.5 * tan(L_sigma_unif[k]);
    L_sigma[k] = L_sigma_unif[k];
  }

  L_Sigma = diag_pre_multiply(L_sigma, L_Omega);

  // Theta[1] = 0;
  // for (k in 1:(K-1)) {
  //   Theta[k + 1] = Theta_raw[k];
  // }

  for (r in 1:R) {
    for (s in 1:S) {
      XB[r,s] = X[r,s] * Beta[r];
    }
  }
}

model {
  //priors
  // Theta_raw ~ normal(0, 10);
  // L_Omega ~ lkj_corr_cholesky(4);

 // Theta_raw ~ normal(0, .2);
 Theta ~ normal(0,1);
 L_Omega ~ lkj_corr_cholesky(50);
 L_sigma_unif ~ lognormal(0,1); # TMN

  //likelihood
  Beta ~ multi_normal_cholesky(Theta, L_Sigma);
  for (r in 1:R) {
    for (s in 1:S) {
       Y[r,s] ~ categorical_logit(XB[r,s]);
      // YW[r,s] ~ categorical_logit(-XB[r,s]);
    }
  }
}
