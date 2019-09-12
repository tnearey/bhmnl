// BayesMNLs2z.stan
// Sum to zero coefficients ('effects coding of response categories;)
//Adapted from From: https://discourse.mc-stan.org/t/speeding-up-a-hierarchical-multinomial-logit-model/1538/4
//  Justin Yap's final 2x faster model
// See yegBaselineHMNL.stan
// Code from YAP is changed as noted. Unmodified code is in YapSpeededHMNLReference.stan
// Feeding data. Should work if we provide data as per // https://mc-stan.org/rstan/reference/stan.html
//       When an element is of type list, it is supposed to make it easier to pass data
//      for those declared in Stan code such as "vector[J] y1[I]" and "matrix[J,K] y2[I]".
//    Using the latter as an example, we can use a list for y2 if the list has "I" elements,
//    each of which is an array (matrix) of dimension "J*K". However, it is not possible
//     to pass a list for data declared such as "vector[K] y3[I,J]"; the only way for it is to
//      use an array with dimension "I*J*K".
// J Yapp's oiginal Data section
// from :  https://discourse.mc-stan.org/t/speeding-up-a-hierarchical-multinomial-logit-model/1538/5
// data {
  //   int<lower=2> C; // Number of alternatives (choices) in each scenario
  //   int<lower=1> K; // Number of // variables in design matrix
  //   int<lower=1> R; // Number of respondents
  //   int<lower=1> S; // Number of scenarios per respondent
  //   int<lower=1,upper=C> YB[R, S]; // best choices
  //   int<lower=1,upper=C> YW[R, S]; // worst choices
  //   matrix[C, K] X[R, S]; // matrix of attributes for each obs
  // }

  // variable changes Yapp -> this file
  // C -> #nRespCat  -- # of alternatives -> responseCategories ;
  // R -> nSj respondent -- # of respondents- subjects/participants/persons
  // K is number of columns in design matrix (same as Yapp )
  // S  is not used.
  // Added
  // nT -> total number of observation trials. Should be equal to  S x R in Yapp's model
  // for fully balianced model (every respondent (participant) gets exactly same set of Scenarios)
  //  where each observation  trial is a single trial for a single participant.
  // we have to identify the Subject from an additoinal data variable
  data {
    int<lower=2> nRespCat; // Number of response category alternatives (choices) in each scenario
    int<lower=1> K; // Number of columns in each design matrix list X[t]
    int<lower=1> nSj; // Number of subjects (respondants/ agents/ persons/perceivers)
    int <lower=1> nT; // total number of observation trials  (grand trials) //// int<lower=1> S; // Number of scenarios per respondent
    int<lower=1,upper=nRespCat> Y[nT]; // YB[nSj, S]; // best choices
    // int<lower=1,upper=nRespCat> YW[nSj, S]; // worst choices
    matrix[nRespCat, K] X[nT]; // was   matrix[nRespCat, K] X[nSj, S]; // matrix of attributes for each obs
    int<lower=1, upper=nSj> SjID[nT]; // Added tmn. Serial identifier for participant on each observation trial
  }

  parameters {
    // vector[K] Beta[nSj];
    vector[K]  Beta[nSj];
    // vector[K - 1] Theta_raw;
    vector[K] Theta;
    cholesky_factor_corr[K] L_Omega;
    vector<lower=0, upper=pi()/2>[K] L_sigma_unif;
  }

  transformed parameters {
    vector<lower=0>[K] L_sigma;
    matrix[K, K] L_Sigma;
    vector[nRespCat] XB[nT]; // was vector[nRespCat] XB[nT] XB[nSj, S];
    // vector[K] Theta;
    // vector[K] Beta[nSj];

    for (k in 1:K) {
      L_sigma[k] = 2.5 * tan(L_sigma_unif[k]);
    }

    L_Sigma = diag_pre_multiply(L_sigma, L_Omega);
    // Theta[1] = 0;
    // for (k in 1:(K-1)) {
      //   Theta[k + 1] = Theta_raw[k];
      // }
      // Assume effects coding(sum to zero)
      // Theta = Theta_raw;
      // for (iSj in 1:nSj){
        //   Beta[iSj] = RM * Beta_raw[iSj];
        // }
        for (t in 1:nT) {// was for (r in 1:nSj) {
          // was:    for (s in 1:S) {
            XB[t] = X[t] * Beta[SjID[t]] ; // was :XB[r,s] = X[r,s] * Beta[r];
            // was:    }
        }
  }

  model {
    //priors (narrowed dowen tmnb)
    Theta ~ normal(0, 2); // changed from norm(0,10)
    L_Omega ~ lkj_corr_cholesky(10);  // changed from lkj_corr_cholesky(4)
    // Was next line missing? -- or is it just defaulted
    L_sigma_unif ~ uniform(0,pi()/2);
    //likelihood
    Beta ~ multi_normal_cholesky(Theta, L_Sigma);
    for (t in 1:nT) { // was for (r in 1:nSj) {
      Y[t]  ~ categorical_logit(XB[t]);

    }
  }
  // was YB[r,s] ~ categorical_logit(XB[r,s]);
      // was:    for (s in 1:S) {
      // was:  YW[r,s] ~ categorical_logit(-XB[r,s]);
      // was: }
  // DONE LIST
  //  1. Remove YW (worst choice) references lines 11 and 56
  //  2. Change YB references to just Y (best/only choice)
  // Replace [nSj,S] dimensioning of X and Y  and XB witn nT
  //     where nT is total number of scenarios for all subjects.
  //     For balanced complete design nT would equal nSj * X
  //
  // Provide a respondant index  PID[nT] that selects current Beta[r] as Beta[PID[n]]


  // provide an PID[nT] vector of indices for respondend.
  // XB
