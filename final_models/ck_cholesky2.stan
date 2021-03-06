# same as ck_cholesky but only GMRF prior on bias (not responsiveness)

functions {
  /** 
  * Transform standard normal to cauchy(location, scale).
  *
  * @param location Location parameter for Cauchy distribution
  * @param scale Scale parameter for Cauchy distribution
  * @param noise Standard normal variable (as declared in parameters block)
  */
    real cauchy_trans_lp(real location, real scale, real noise) {
      noise ~ normal(0,1) ;
      return location + scale * tan(pi() * (Phi_approx(noise) - 0.5)) ;
    }
  
  /** 
  * Transform standard normal to normal(location, scale).
  *
  * @param location Location parameter for normal distribution
  * @param scale Scale parameter for normal distribution
  * @param noise Standard normal variable (as declared in parameters block)
  */
    real normal_trans_lp(real location, real scale, real noise) {
      noise ~ normal(0,1) ;
      return location + scale * noise ;
    }
    
  /** 
  * Transform vector of standard normals to vector of normal(location, scale).
  *
  * @param location Location parameter for normal distribution
  * @param scale Scale parameter for normal distribution
  * @param noise Standard normal variable (as declared in parameters block)
  */
    vector normal_trans_vec_lp(real location, real scale, vector noise) {
      noise ~ normal(0,1) ;
      return location + scale * noise ;
    }
  
}

data {
  // dimensions 
  int<lower=1>              N ; # number of observations 
  int<lower=1>              C ; # number of congresses/periods
  
  // variables
  int<lower=1,upper=C>      congress[N] ; 
  int<lower=1,upper=56>     nvotes[N] ;
  int<lower=0,upper=55>     majps[N] ;
  real                      lnmajvavg[N] ;
  
  // stuff for spatial-smoothing priors
  matrix[C,C]               SIGMA ; # covariance matrix
}
transformed data {
  matrix[C,C] OMEGA ;
  OMEGA <- cholesky_decompose(SIGMA) ;
}
parameters {
  vector[C]             bias_noise;
  real                  resp_noise;
  real<lower=0>         tau_b_noise ;
  real<lower=0>         phi ;
  # real                  Const_noise ;
}
transformed parameters {
  vector[C]             b_bias ;
  real                  b_resp ;
  real<lower=0>         tau_b ;
  #  real                  Const ;
  
  # Const <- normal_trans_lp(0, 5, Const_noise) ;
  tau_b <- cauchy_trans_lp(0, 2.5, tau_b_noise) ;
  b_resp <- normal_trans_lp(0, 2.5, resp_noise) ;
  b_bias  <- (tau_b * OMEGA) * bias_noise ;
}
model {
  // local variables
  vector<lower=0>[N]  alphas ;   
  vector<lower=0>[N]  betas ;    
  
  bias_noise ~ normal(0,1) ;
  phi ~ gamma(0.0001, 0.0001) ;
  
  // likelihood
  for (n in 1:N) {
    real theta_n ;
    theta_n <- inv_logit(b_bias[congress[n]] + b_resp*lnmajvavg[n]) ;    
    alphas[n] <- theta_n * phi ;
    betas[n]  <- (1 - theta_n) * phi ;
  }
  
  majps ~ beta_binomial(nvotes, alphas, betas) ;
}
generated quantities {
  real        bias[C] ;
  vector[N]   log_lik ; # for computing information criteria, etc. 
  int         y_rep[N] ; # for posterior predictive checking
  
  for (c in 1:C) 
    bias[c] <- inv_logit(b_bias[c]) - 0.5 ;
  
  for (n in 1:N) {
    real alpha_n ;
    real beta_n ;
    alpha_n <- phi * inv_logit(b_bias[congress[n]] + b_resp*lnmajvavg[n]) ;
    beta_n  <- phi * (1 - inv_logit(b_bias[congress[n]] + b_resp*lnmajvavg[n])) ;
    log_lik[n] <- beta_binomial_log(majps[n], nvotes[n], alpha_n, beta_n) ;
    y_rep[n] <- beta_binomial_rng(nvotes[n], alpha_n, beta_n) ;
  }
  
}




