# estimating omega in (D - omega*A)

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
  * Transform vector of standard normals to vector of Cauchys 
  * with same location and scale.
  *
  * @param location Location parameter for Cauchy distribution
  * @param scale Scale parameter for Cauchy distribution
  * @param noise Vector of standard normals (as declared in parameters block)
  */
    vector cauchy_trans_vec1_lp(real location, real scale, vector noise) {
      vector[num_elements(noise)] out ;
      noise ~ normal(0,1) ;
      for (j in 1:num_elements(out)) {
        out[j] <- location + scale * tan(pi() * (Phi_approx(noise[j]) - 0.5)) ;
      }
      return out ;
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
  matrix[C,C]               A ; # adjacency matrix
  matrix[C,C]               D ; # degree matrix
}
parameters {
  real                  Const ; # global intercept
  vector[C]             b_bias ;
  vector[C]             b_resp ;
  real<lower=0>         tau_b_noise ;
  real<lower=0>         tau_r_noise ;
  real<lower=0>         phi ;
  real<lower=0,upper=1> omega ;
}
transformed parameters {
  real<lower=0>         tau_b ;
  real<lower=0>         tau_r ;

  tau_b <- cauchy_trans_lp(0, 2.5, tau_b_noise) ;
  tau_r <- cauchy_trans_lp(0, 2.5, tau_r_noise) ;
}
model {
  // local variables
  vector<lower=0>[N]  alphas ;   
  vector<lower=0>[N]  betas ;  
  matrix[C,C] PREC ;
  
  PREC <- D - omega * A ;
  b_bias ~ multi_normal_prec(rep_vector(0,C), PREC/tau_b) ;
  b_resp ~ multi_normal_prec(rep_vector(0,C), PREC/tau_r) ;
  
  phi ~ gamma(0.0001, 0.0001) ;
  omega ~ beta(100, 1) ;
  
  // likelihood
  for (n in 1:N) {
    real theta_n ;
    theta_n <- inv_logit(Const + b_bias[congress[n]] + b_resp[congress[n]]*lnmajvavg[n]) ;    
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
    real theta_n ;
    real alpha_n ;
    real beta_n ;
    theta_n <- inv_logit(Const + b_bias[congress[n]] + b_resp[congress[n]]*lnmajvavg[n]) ;
    alpha_n <- phi * theta_n ;
    beta_n  <- phi * (1 - theta_n) ;
    log_lik[n] <- beta_binomial_log(majps[n], nvotes[n], alpha_n, beta_n) ;
    y_rep[n] <- beta_binomial_rng(nvotes[n], alpha_n, beta_n) ;
  }
  
}

