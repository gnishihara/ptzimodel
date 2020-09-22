// generated with brms 2.13.9
functions {
  /* zero-inflated beta log-PDF of a single response 
   * Args: 
   *   y: the response value 
   *   mu: mean parameter of the beta distribution
   *   phi: precision parameter of the beta distribution
   *   zi: zero-inflation probability
   * Returns:  
   *   a scalar to be added to the log posterior 
   */ 
   real zero_inflated_beta_lpdf(real y, real mu, real phi, real zi) {
     row_vector[2] shape = [mu * phi, (1 - mu) * phi];
     if (y == 0) { 
       return bernoulli_lpmf(1 | zi); 
     } else { 
       return bernoulli_lpmf(0 | zi) +  
              beta_lpdf(y | shape[1], shape[2]); 
     } 
   }
  /* zero-inflated beta log-PDF of a single response 
   * logit parameterization of the zero-inflation part
   * Args: 
   *   y: the response value 
   *   mu: mean parameter of the beta distribution
   *   phi: precision parameter of the beta distribution
   *   zi: linear predictor for zero-inflation part
   * Returns:  
   *   a scalar to be added to the log posterior 
   */ 
   real zero_inflated_beta_logit_lpdf(real y, real mu, real phi, real zi) {
     row_vector[2] shape = [mu * phi, (1 - mu) * phi];
     if (y == 0) { 
       return bernoulli_logit_lpmf(1 | zi); 
     } else { 
       return bernoulli_logit_lpmf(0 | zi) +  
              beta_lpdf(y | shape[1], shape[2]); 
     } 
   }
  // zero-inflated beta log-CCDF and log-CDF functions
  real zero_inflated_beta_lccdf(real y, real mu, real phi, real zi) { 
    row_vector[2] shape = [mu * phi, (1 - mu) * phi];
    return bernoulli_lpmf(0 | zi) + beta_lccdf(y | shape[1], shape[2]); 
  }
  real zero_inflated_beta_lcdf(real y, real mu, real phi, real zi) { 
    return log1m_exp(zero_inflated_beta_lccdf(y | mu, phi, zi));
  }

  real fvfmmodel (real ps, real ha, real eta, real ktopt, real temperature) {
    real inverse_kelvin = 1.0 / (temperature + 273.15);
    real gas_constant = 8.314/1000.0;
    real hd = ha * eta;
    real x = (1.0 / ktopt - inverse_kelvin);
    real numerator = ha / gas_constant * x;
    real denominator = hd / gas_constant * x;
    return ps * hd * exp(numerator) / (hd - ha * (1.0 - exp(denominator)));
  }

}
data {
  int<lower=1> N;  // number of observations
  vector[N] Y;  // response variable
  int<lower=1> K_PS;  // number of population-level effects
  matrix[N, K_PS] X_PS;  // population-level design matrix
  int<lower=1> K_HA;  // number of population-level effects
  matrix[N, K_HA] X_HA;  // population-level design matrix
  int<lower=1> K_ET;  // number of population-level effects
  matrix[N, K_ET] X_ET;  // population-level design matrix
  int<lower=1> K_KT;  // number of population-level effects
  matrix[N, K_KT] X_KT;  // population-level design matrix
  // covariate vectors for non-linear functions
  vector[N] C_1;
  int<lower=1> K_zi;  // number of population-level effects
  matrix[N, K_zi] X_zi;  // population-level design matrix
  // data for group-level effects of ID 1
  int<lower=1> N_1;  // number of grouping levels
  int<lower=1> M_1;  // number of coefficients per level
  int<lower=1> J_1[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_1_PS_1;
  // data for group-level effects of ID 2
  int<lower=1> N_2;  // number of grouping levels
  int<lower=1> M_2;  // number of coefficients per level
  int<lower=1> J_2[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_2_HA_1;
  // data for group-level effects of ID 3
  int<lower=1> N_3;  // number of grouping levels
  int<lower=1> M_3;  // number of coefficients per level
  int<lower=1> J_3[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_3_ET_1;
  // data for group-level effects of ID 4
  int<lower=1> N_4;  // number of grouping levels
  int<lower=1> M_4;  // number of coefficients per level
  int<lower=1> J_4[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_4_KT_1;
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
  int Kc_zi = K_zi - 1;
  matrix[N, Kc_zi] Xc_zi;  // centered version of X_zi without an intercept
  vector[Kc_zi] means_X_zi;  // column means of X_zi before centering
  for (i in 2:K_zi) {
    means_X_zi[i - 1] = mean(X_zi[, i]);
    Xc_zi[, i - 1] = X_zi[, i] - means_X_zi[i - 1];
  }
}
parameters {
  vector<lower=0>[K_PS] b_PS;  // population-level effects
  vector<lower=0>[K_HA] b_HA;  // population-level effects
  vector<lower=1>[K_ET] b_ET;  // population-level effects
  vector<lower=273.15>[K_KT] b_KT;  // population-level effects
  real<lower=0> phi;  // precision parameter
  vector[Kc_zi] b_zi;  // population-level effects
  real Intercept_zi;  // temporary intercept for centered predictors
  vector<lower=0>[M_1] sd_1;  // group-level standard deviations
  vector[N_1] z_1[M_1];  // standardized group-level effects
  vector<lower=0>[M_2] sd_2;  // group-level standard deviations
  vector[N_2] z_2[M_2];  // standardized group-level effects
  vector<lower=0>[M_3] sd_3;  // group-level standard deviations
  vector[N_3] z_3[M_3];  // standardized group-level effects
  vector<lower=0>[M_4] sd_4;  // group-level standard deviations
  vector[N_4] z_4[M_4];  // standardized group-level effects
}
transformed parameters {
  vector[N_1] r_1_PS_1;  // actual group-level effects
  vector[N_2] r_2_HA_1;  // actual group-level effects
  vector[N_3] r_3_ET_1;  // actual group-level effects
  vector[N_4] r_4_KT_1;  // actual group-level effects
  r_1_PS_1 = (sd_1[1] * (z_1[1]));
  r_2_HA_1 = (sd_2[1] * (z_2[1]));
  r_3_ET_1 = (sd_3[1] * (z_3[1]));
  r_4_KT_1 = (sd_4[1] * (z_4[1]));
}
model {
  // initialize linear predictor term
  vector[N] nlp_PS = X_PS * b_PS;
  // initialize linear predictor term
  vector[N] nlp_HA = X_HA * b_HA;
  // initialize linear predictor term
  vector[N] nlp_ET = X_ET * b_ET;
  // initialize linear predictor term
  vector[N] nlp_KT = X_KT * b_KT;
  // initialize non-linear predictor term
  vector[N] mu;
  // initialize linear predictor term
  vector[N] zi = Intercept_zi + Xc_zi * b_zi;
  for (n in 1:N) {
    // add more terms to the linear predictor
    nlp_PS[n] += r_1_PS_1[J_1[n]] * Z_1_PS_1[n];
  }
  for (n in 1:N) {
    // add more terms to the linear predictor
    nlp_HA[n] += r_2_HA_1[J_2[n]] * Z_2_HA_1[n];
  }
  for (n in 1:N) {
    // add more terms to the linear predictor
    nlp_ET[n] += r_3_ET_1[J_3[n]] * Z_3_ET_1[n];
  }
  for (n in 1:N) {
    // add more terms to the linear predictor
    nlp_KT[n] += r_4_KT_1[J_4[n]] * Z_4_KT_1[n];
  }
  for (n in 1:N) {
    // compute non-linear predictor values
    mu[n] = fvfmmodel(nlp_PS[n] , nlp_HA[n] , nlp_ET[n] , nlp_KT[n] , C_1[n]);
  }
  // priors including all constants
  target += normal_lpdf(b_PS | 0.608571428571429, 0.122982676324587)
    - 1 * normal_lccdf(0 | 0.608571428571429, 0.122982676324587);
  target += normal_lpdf(b_HA | 28.522619047619, 33.3237790334597)
    - 1 * normal_lccdf(0 | 28.522619047619, 33.3237790334597);
  target += normal_lpdf(b_ET | 24.1204761904762, 70.4801524372116)
    - 1 * normal_lccdf(1 | 24.1204761904762, 70.4801524372116);
  target += normal_lpdf(b_KT | 290.461904761905, 4.8592482797617)
    - 1 * normal_lccdf(273.15 | 290.461904761905, 4.8592482797617);
  target += student_t_lpdf(phi | 3, 0, 2)
    - 1 * student_t_lccdf(0 | 3, 0, 2);
  target += student_t_lpdf(b_zi | 3, 0, 2);
  target += student_t_lpdf(Intercept_zi | 3, 0, 2);
  target += student_t_lpdf(sd_1 | 3, 0, 2)
    - 1 * student_t_lccdf(0 | 3, 0, 2);
  target += std_normal_lpdf(z_1[1]);
  target += student_t_lpdf(sd_2 | 3, 0, 2)
    - 1 * student_t_lccdf(0 | 3, 0, 2);
  target += std_normal_lpdf(z_2[1]);
  target += student_t_lpdf(sd_3 | 3, 0, 2)
    - 1 * student_t_lccdf(0 | 3, 0, 2);
  target += std_normal_lpdf(z_3[1]);
  target += student_t_lpdf(sd_4 | 3, 0, 2)
    - 1 * student_t_lccdf(0 | 3, 0, 2);
  target += std_normal_lpdf(z_4[1]);
  // likelihood including all constants
  if (!prior_only) {
    for (n in 1:N) {
      target += zero_inflated_beta_logit_lpdf(Y[n] | mu[n], phi, zi[n]);
    }
  }
}
generated quantities {
  // actual population-level intercept
  real b_zi_Intercept = Intercept_zi - dot_product(means_X_zi, b_zi);
}
