# same as ck_cholesky but with global intercept

data {
  // dimensions 
  int<lower=1>              N ; # number of observations 
  int<lower=1>              C ; # number of congresses (time periods)
  
  // variables
  int<lower=1,upper=C>      cong[N] ;   # maps between observations & congresses
  int<lower=1,upper=56>     nVotes[N] ; # numerber of votes
  int<lower=0,upper=55>     nWins[N] ;  # number of majority party victories
  real                      vRatio[N] ; # vRATIO
  
  // inverse of penalty matrix for GMRF prior
  matrix[C,C]               Pinverse ; 
}
transformed data {
  real<lower=0> tau_scale ;   # scale for Cauchy priors on taus
  real<lower=0> tau_loc ;     # location for Cauchy priors on taus
  matrix[C,C] cholPinverse ;  # Cholesky decomposition of Pinverse
  
  tau_loc <- 0.0 ;
  tau_scale <- 2.5 ;
  cholPinverse <- cholesky_decompose(Pinverse) ;
}
parameters {
  real                  Const ; # common intercept
  real<lower=0>         phi ;
  vector[C]             bias_noise ;
  vector[C]             resp_noise ;
  real<lower=0>         tau_bias_noise ;
  real<lower=0>         tau_resp_noise ;
}
transformed parameters {
  vector[C]             b_bias ;
  vector[C]             b_resp ;
  real<lower=0>         tau_bias ;
  real<lower=0>         tau_resp ;
  
  tau_bias <- tau_loc + tau_scale * tan(pi() * (Phi_approx(tau_bias_noise) - 0.5)) ;
  tau_resp <- tau_loc + tau_scale * tan(pi() * (Phi_approx(tau_resp_noise) - 0.5)) ;
  
  b_bias  <- (tau_b * cholPinverse) * bias_noise ;
  b_resp  <- (tau_r * cholPinverse) * resp_noise ;
}
model {
  // local variables
  real                logLik ;    # log likelihood
  real                logPrior ;  # log prior
  vector<lower=0>[N]  alphas ;   
  vector<lower=0>[N]  betas ;    
  
  logPrior <- (
    gamma_log(phi, 0.0001, 0.0001) +
    normal_log(bias_noise, 0, 1) +  // vectorized
    normal_log(resp_noise, 0, 1) +  // vectorized   
    normal_log(tau_bias_noise, 0, 1) +
    normal_log(tau_resp_noise, 0, 1)
  ) ;
  
  // likelihood
  for (n in 1:N) {
    real theta_n ;
    theta_n   <- inv_logit(Const + b_bias[cong[n]] + b_resp[cong[n]] * vRatio[n]) ;    
    alphas[n] <- theta_n * phi ;
    betas[n]  <- (1 - theta_n) * phi ;
  }
  logLik <- beta_binomial_log(nWins, nVotes, alphas, betas) ; // vectorized
  
  increment_log_prob(logPrior + logLik) ; 
}
generated quantities {
  real        bias[C] ;
  int         y_rep[N] ; # for posterior predictive checking
  
  for (c in 1:C) 
    bias[c] <- inv_logit(b_bias[c]) - 0.5 ;
  
  for (n in 1:N) {
    real theta_n ;
    real alpha_n ;
    real beta_n ;
    theta_n   <- inv_logit(Const + b_bias[cong[n]] + b_resp[cong[n]] * vRatio[n]) ;
    alpha_n   <- phi * theta_n ;
    beta_n    <- phi * (1 - theta_n) ;
    y_rep[n]  <- beta_binomial_rng(nVotes[n], alpha_n, beta_n) ;
  }
}

