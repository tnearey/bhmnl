// BayesMNL_s2zNC.stan
// This one goes too far. We move in this direction
// ?/ from  BayesMNL_s2z thru BayesMNL_s2zM<x>.stan
// Non centered everywhere SUG_MVHPriors
// We had problems with looser priors.
// Tryingn ***
// tauHN ~ norm(0,1) instead of tau ~ 2.5*tan(uniform(0,pi))
//  L_Omega ~ lkj_corr_cholesky(5); (5 instead of 2 as lkj par)
// theta prior sd at 3 instead of 5. theta ~ normal(0,3)
// ***
// https://mc-stan.org/docs/2_19/stan-users-guide/multivariate-hierarchical-priors-section.html
// Sum to zero coefficients ('effects coding of response categories;)
//Adapted from From: https://discourse.mc-stan.org/t/speeding-up-a-hierarchical-multinomial-logit-model/1538/4
//  Justin Yap's final 2x faster model
// Yap partly follows  Stan Users Guide  2_19 [SUGMVHP]
// We try to take it here as far as we can.
//https://mc-stan.org/docs/2_19/stan-users-guide/multivariate-hierarchical-priors-section.html
// and binary hierarchial logistic models.
data {
  int<lower=0> N;              // num individuals <<total trials nT>>
  int<lower=1> D;              // num ind predictors
  int<lower=1> J;              // num groups ( Subjects)
  int<lower=2> C;  // number of response categories
  // int<lower=1> L;              // num group predictors we don't have any
  int<lower=1,upper=J> jj[N];  // group <subject> for individual <trial>
  int<lower=0,upper=1> dbgLev;

  // matrix[N, D] x;               // individual predictors <upper cased TMN>
  matrix[C, D] X[N];
  // row_vector[L] u[J];          // group predictors
  // vector[N] y;                 // outcomes
  //Changed:

  int<lower=1,upper=C> Y[N];

}
parameters {
  vector[D] z[J];
  cholesky_factor_corr[D] L_Omega;
  // ***vector<lower=0,upper=pi()/2>[D] tau_unif;
  // ***TMN Half normal tau --> tauHN
  vector<lower=0>[D] tauHN;
  // matrix[L, D] gamma;                         // group coeffs
  vector[D] theta;
  real<lower=0> sigma;                       // prediction error scale
}
transformed parameters {
  // matrix[J, D] beta;
  vector[D] beta [J];
  //*** vector<lower=0>[D] tau;     // prior scale
  vector[C] XB[N]; // was vector[nRespCat] XB[nT] XB[nSj, S];

  // for (k in 1:D) tau[k] = 2.5 * tan(tau_unif[k]);
   // This is pretty tricky don't follow args in SUG_MVHPriors
  // beta = u * gamma + (diag_pre_multiply(tau,L_Omega) * z)';

  for (j in 1:J){
    //*** beta[j] = (diag_pre_multiply(tau,L_Omega) * z[j]);
  beta[j] = (diag_pre_multiply(tauHN,L_Omega) * z[j]);
 if (dbgLev>0){
     for( k in 1:D)
       print("beta ", j," ",k," ", beta[j][k]);
 }
  }
  // we're adding theta here

  for (t in 1:N) {
    // XB[t] = X[t] * Beta[SjID[t]] ; // was :XB[r,s] = X[r,s] * Beta[r];
    // Note xB has done all its j indexings
    XB[t] =X[t]*(theta+beta[jj[t]]) ;
  }
}

model {
  for (j in 1:J){
    z[j] ~ normal(0,1);
  }
  // to_vector(z) ~ std_normal();
  //*** L_Omega ~ lkj_corr_cholesky(2);
  // *** Constraining closer to identity Omega
  L_Omega ~ lkj_corr_cholesky(5);

  // to_vector(theta) ~ normal(0, 5);
  tauHN ~ normal(0,1);
  theta ~ normal(0,3);
  for (t in 1:N) {   // was for (r in 1:nSj) {
    Y[t]  ~ categorical_logit(XB[t]);
  }
  // y ~ normal(rows_dot_product(beta[jj] , x), sigma);
}
