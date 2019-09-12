// BayesMNLs2zM1.stan
// Mark 1-- made beta deviation from theta
// See also ...M1A wth transformed vaariabls block moved to Model block
//
// N << nT  total number of trials
// -> J; nSj
// jj[] << SjID[] ->
// -> C nRespCat
//
//
// Sum to zero coefficients ('effects coding of response categories;)
//Adapted from From: https://discourse.mc-stan.org/t/speeding-up-a-hierarchical-multinomial-logit-model/1538/4
//  Justin Yap's final 2x faster model
// This seems to follow Stan Users Guide  2_19 [SUGMVHP]
//https://mc-stan.org/docs/2_19/stan-users-guide/multivariate-hierarchical-priors-section.html
// and binary hierarchial logistic models.
// See yegBaselineHMNL.stan
// Code from YAP is changed as noted. Unmodified code is in YapSpeededHMNLReference.stan
// Feeding data. Should work if we provide data as per // https://mc-stan.org/rstan/reference/stan.html
//       When an element is of type list, it is supposed to make it easier to pass data
//      for those declared in Stan code such as "vector[J] y1[I]" and "matrix[J,D] y2[I]".
//    Using the latter as an example, we can use a list for y2 if the list has "I" elements,
//    each of which is an array (matrix) of dimension "J*D". However, it is not possible
//     to pass a list for data declared such as "vector[D] y3[I,J]"; the only way for it is to
//      use an array with dimension "I*J*D".
// J Yapp's oiginal Data section
// from :  https://discourse.mc-stan.org/t/speeding-up-a-hierarchical-multinomial-logit-model/1538/5
// data {
  //   int<lower=2> C; // Number of alternatives (choices) in each scenario
  //   int<lower=1> D; // Number of // variables in design matrix
  //   int<lower=1> R; // Number of respondents
  //   int<lower=1> S; // Number of scenarios per respondent
  //   int<lower=1,upper=C> YB[R, S]; // best choices
  //   int<lower=1,upper=C> YW[R, S]; // worst choices
  //   matrix[C, D] X[R, S]; // matrix of attributes for each obs
  // }

  // variable changes Yapp -> this file
  // C -> #C -> C  -- # of alternatives -> responseCategories ;
  // R -> nSj -> J respondent -- # of respondents- subjects/participants/persons
  // D is number of columns in design matrix (same as Yapp )
  // S  is not used.
  // Added
  // nT -> N  total number of observation trials. Should be equal to  S x R in Yapp's model
  // for fully balianced model (every respondent (participant) gets exactly same set of Scenarios)
  //  where each observation  trial is a single trial for a single participant.
  // we have to identify the Subject from an additoinal data variable
  data {
    int<lower=2> C; // Number of response category alternatives (choices) in each scenario
    int<lower=1> D; // Number of columns in each design matrix list X[t]
    int<lower=1> J; // Number of subjects (respondants/ agents/ persons/perceivers)
    int <lower=1> N; // total number of observation trials  (grand trials) //// int<lower=1> S; // Number of scenarios per respondent
    int<lower=1,upper=C> Y[N]; // YB[J, S]; // best choices
    // int<lower=1,upper=C> YW[J, S]; // worst choices
    matrix[C, D] X[N]; // was   matrix[C, D] X[J, S]; // matrix of attributes for each obs
    int<lower=1, upper=J> jj[N]; // Added tmn. Serial identifier for participant on each observation trial
    int<lower=0,upper=1> dbgLev;
  }

  parameters {
    // vector[D] beta[J];
    vector[D]  beta[J];
    // vector[D - 1] Theta_raw;
    vector[D] theta;
    cholesky_factor_corr[D] L_Omega;
    vector<lower=0, upper=pi()/2>[D] L_sigma_unif;
  }

  transformed parameters {
    vector[D] DZeros;
    vector<lower=0>[D] L_sigma;
    matrix[D, D] L_Sigma;
    vector[C] XB[N]; // was vector[C] XB[N] XB[J, S];
    // vector[D] theta;
    // vector[D] beta[J];
    // DZeros = rep_vector(0.0,D); NOT EEDED
    // for (k in 1:D) {
    //   L_sigma[k] = 2.5 * tan(L_sigma_unif[k]);
    // } Does  vectorized version work? M1
    L_sigma=2.5*tan(L_sigma_unif);

    L_Sigma = diag_pre_multiply(L_sigma, L_Omega);
    // theta[1] = 0;
    // for (k in 1:(D-1)) {
      //   theta[k + 1] = Theta_raw[k];
      // }
      // Assume effects coding(sum to zero)
      // theta = Theta_raw;
      // for (iSj in 1:J){
        //   beta[iSj] = RM * Beta_raw[iSj];
        // }
        for (t in 1:N) {
            // XB[t] = X[t] * beta[jj[t]] ; // was :XB[r,s] = X[r,s] * beta[r];
              XB[t] = X[t] * (beta[jj[t]]+theta) ; // Mark1
              // if( t==1)
              //   for (c in 1:C)
              //       print("t ", t, "XB[t][] ",c," ", XB[t][c])
        }
  }

  model {
    //priors (narrowed dowen tmnb)
    theta ~ normal(0, 5); //  M1 used 2 changed from norm(0,10)
    L_Omega ~ lkj_corr_cholesky(4);  //  M1 used 10 changed from lkj_corr_cholesky(4)
    // Not necessary see SUGMVHP the limits work
    // L_sigma_unif ~ uniform(0,pi()/2);
    //likelihood
    // beta ~ multi_normal_cholesky(theta, L_Sigma);
    beta ~ multi_normal_cholesky(rep_vector(0.0,D), L_Sigma);
    for (t in 1:N) { // was for (r in 1:J) {
      Y[t]  ~ categorical_logit(XB[t]);
      // if (t==1)
      //   print("t ", t, "Y[t] ", Y[t])

    }
  }
  // was YB[r,s] ~ categorical_logit(XB[r,s]);
      // was:    for (s in 1:S) {
      // was:  YW[r,s] ~ categorical_logit(-XB[r,s]);
      // was: }
  // DONE LIST
  //  1. Remove YW (worst choice) references lines 11 and 56
  //  2. Change YB references to just Y (best/only choice)
  // Replace [J,S] dimensioning of X and Y  and XB witn N
  //     where N is total number of scenarios for all subjects.
  //     For balanced complete design N would equal J * X
  //
  // Provide a respondant index  PID[N] that selects current beta[r] as beta[PID[n]]


  // provide an PID[N] vector of indices for respondend.
  // XB
