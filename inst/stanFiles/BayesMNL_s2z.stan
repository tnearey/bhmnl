// BayesMNL_s2z.stan
// Simple multinomial (non hierarchal)
// We're sticking closer to Stan users guide 2-18  section 1.
// We're using C instead of D as number of categories.
// and our former D predictors are now D dimensions of design matrix
// but we're allowing for custom (alternative-specific like)
// design matrices to accomodate factorial responses
data {
  int<lower=0> N;              // num individuals <<total trials nT>>
  int<lower=1> D;              // num ind predictors
  int<lower=1> J;              // num groups ( Subjects) // NOT USED
  int<lower=2> C;  // number of response categories
  int<lower=1,upper=J> jj[N];  // NOT YET USED group <subject> for individual <trial>
  int<lower=0,upper=1> dbgLev;

  // matrix[N, D] x;               // individual predictors <upper cased TMN>
  matrix[C,D] X[N];
  // row_vector[L] u[J];          // group predictors
  // vector[N] y;                 // outcomes
  //Changed:

  int<lower=1,upper=C> Y[N];

}
parameters {
  // vector[D] z[J];
  // cholesky_factor_corr[D] L_Omega;
  // ***vector<lower=0,upper=pi()/2>[D] tau_unif;
  // ***TMN Half normal tau --> tauHN
  // vector<lower=0>[D] tauHN;
  // matrix[L, D] gamma;                         // group coeffs
  vector[D] theta;
  // real<lower=0> sigma;                       // prediction error scale
}
transformed parameters {
  // matrix[J, D] beta;
  vector[D] beta [J];
  //*** vector<lower=0>[D] tau;     // prior scale
  vector[C] XB[N]; // was vector[nRespCat] XB[nT] XB[nSj, S];

  // for (k in 1:D) tau[k] = 2.5 * tan(tau_unif[k]);
   // This is pretty tricky don't follow args in SUG_MVHPriors
  // beta = u * gamma + (diag_pre_multiply(tau,L_Omega) * z)';

  // for (j in 1:J){
    //*** beta[j] = (diag_pre_multiply(tau,L_Omega) * z[j]);
  // beta[j] = (diag_pre_multiply(tauHN,L_Omega) * z[j]);
 // if (dbgLev>0){
 //     for( k in 1:D)
 //       print("beta ", j," ",k," ", beta[j][k]);
 // }
  // }
  for (t in 1:N) {
    // XB[t] = X[t] * Beta[SjID[t]] ; // was :XB[r,s] = X[r,s] * Beta[r];
    // Note xB has done all its j indexings
    // XB[t] =X[t]*(theta+beta[jj[t]]) ; NOT YET
    XB[t]= X[t]*theta;
  }
}

model {
  // for (j in 1:J){
  //   z[j] ~ normal(0,1);
  // }
  // to_vector(z) ~ std_normal();
  //*** L_Omega ~ lkj_corr_cholesky(2);
  // *** Constraining closer to identity Omega
  // L_Omega ~ lkj_corr_cholesky(5);

  // to_vector(theta) ~ normal(0, 5);
  // tauHN ~ normal(0,1);
  theta ~ normal(0,3);
  for (t in 1:N) {   // was for (r in 1:nSj) {
    Y[t]  ~ categorical_logit(XB[t]);
  }
}
