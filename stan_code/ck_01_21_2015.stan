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
}

data {
  int<lower=1>            N ;
  int<lower=1>            C ;
  int<lower=0>            congress[N] ;
  int<lower=0,upper=56>   nvotes[N] ;
  int<lower=0,upper=55>   majps[N] ;
  real                    lnmajvavg[N] ;
  matrix[C,C]             PREC ;
}

parameters {
  real lambda_bar_noise ;
  real<lower=0> rho_bar ;
  vector[C - 2] free_lambdas ;
  vector[C - 2] free_rhos ;
  real<lower=0> tau_lambda_noise ;
  real<lower=0> tau_rho_noise ; 
  real<lower=0> phi ;
}

transformed parameters {
  real lambda_bar ;
  real<lower=0> tau_lambda ;
  real<lower=0> tau_rho ;
  vector[2] pinned_lambdas ;
  vector[2] pinned_rhos ;

  lambda_bar <- cauchy_trans_lp(0.0, 2.5, lambda_bar_noise) ;
  tau_lambda <- cauchy_trans_lp(0.0, 2.5, tau_lambda_noise) ;
  tau_rho <- cauchy_trans_lp(0.0, 2.5, tau_rho_noise) ;

  pinned_lambdas[1] <- 2 * free_lambdas[C - 2] - free_lambdas[C - 3] ;
  pinned_lambdas[2]  <- 3 * free_lambdas[C - 2] - 2 * free_lambdas[C - 3] ;
  pinned_rhos[1] <- 2 * free_rhos[C - 2] - free_rhos[C - 3] ;
  pinned_rhos[2]  <- 3 * free_rhos[C - 2] - 2 * free_rhos[C - 3] ;
}

model {
  vector[C] all_lambdas ;
  vector[C] all_rhos ;
  vector[N] prior_samp_sizes1 ;  // for beta_binomial dist   
  vector[N] prior_samp_sizes2 ;  // for beta_binomial dist
  
  // priors on phi and rho_bar
  phi ~ gamma(0.0001, 0.0001) ;
  rho_bar ~ gamma(1.0,1.0) ;

  // fill in all_lambdas and all_rhos
#  all_lambdas <- append_row(free_lambdas, pinned_lambdas) ;
#  all_rhos <- append_row(free_rhos, pinned_rhos) ;

  for (c in 1:(C - 2)) {
    all_rhos[c] <- free_rhos[c] ;
    all_lambdas[c] <- free_lambdas[c] ;
  }
  all_rhos[C - 1] <- pinned_rhos[1] ;
  all_rhos[C] <- pinned_rhos[2] ;
  all_lambdas[C - 1] <- pinned_lambdas[1] ;
  all_lambdas[C] <- pinned_lambdas[2] ;


  // multivariate normal priors for lambda and rho
  increment_log_prob(- 0.5 * tau_lambda * quad_form(PREC, all_lambdas)) ;
  increment_log_prob(- 0.5 * tau_rho * quad_form(PREC, all_rhos)) ;
  

  // likelihood
  for (n in 1:N) { 
    real theta ;
    theta <- inv_logit(lambda_bar + all_lambdas[congress[n]] + (rho_bar + all_rhos[congress[n]]) * lnmajvavg[n]) ;
    prior_samp_sizes1[n] <- theta * phi ;
    prior_samp_sizes2[n] <- (1 - theta) * phi ;
  }

  majps ~ beta_binomial(nvotes, prior_samp_sizes1, prior_samp_sizes2) ;

}

generated quantities {
  real        bias[C];
  vector[N]   log_lik ;  # for computing information criteria, etc. 
  int         y_rep[N] ; # for posterior predictive checking
  
  { # local
  vector[C] lambdas_temp ;
  vector[C] rhos_temp ;

  for (c in 1:(C - 2)) {
    rhos_temp[c] <- free_rhos[c] ;
    lambdas_temp[c] <- free_lambdas[c] ;
  }
  rhos_temp[C - 1] <- pinned_rhos[1] ;
  rhos_temp[C] <- pinned_rhos[2] ;
  lambdas_temp[C - 1] <- pinned_lambdas[1] ;
  lambdas_temp[C] <- pinned_lambdas[2] ;



  for (c in 1:C) 
    bias[c] <- inv_logit(lambda_bar + lambdas_temp[c]) - 0.5;

  for (n in 1:N) {
    real theta_n ;
    real alpha_n ;
    real beta_n ;
    theta_n <- inv_logit(lambda_bar + lambdas_temp[congress[n]] + (rho_bar + rhos_temp[congress[n]]) * lnmajvavg[n]) ;
    alpha_n <- phi * theta_n  ;
    beta_n  <- phi * (1 - theta_n) ;
    log_lik[n] <- beta_binomial_log(majps[n], nvotes[n], alpha_n, beta_n) ;
    y_rep[n] <- beta_binomial_rng(nvotes[n], alpha_n, beta_n) ;
  }

  } # end local

}
