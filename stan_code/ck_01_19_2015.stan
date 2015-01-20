# notes: 2 rhos, 1 phi
# added constant
# RW1

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
  matrix[C, C]              A ;         # adjacency matrix
  vector[C]                 d ;         # diagonal vector of degree matrix
  int                       A_N ;       # number of adjacent region pairs
  int                       A1[A_N] ;   # 1st half of adjacency pairs (rows)
  int                       A2[A_N] ;   # 2nd half of adjacency pairs (cols)
}

transformed data { # deterministic functions of data to be called only once
  real                  b_bias_mean ;
  real                  b_resp_mean ;
  vector[C] sqrt_d ; # sqrt(diag elements of degree mtrx)

  b_bias_mean <- 0.0 ;
  b_resp_mean <- 0.0 ;
  for (c in 1:C) {
    sqrt_d[c] <- sqrt(d[c]) ;
  }
}
parameters { 
  vector[C]             b_bias ;
  vector[C]             b_resp ;
  real<lower=0>         tau_sq_b ;
  real<lower=0>         tau_sq_r ;
  real<lower=0>         phi ;
  real<lower=0,upper=1> rho_bias ;
  real<lower=0,upper=1> rho_resp ;
  real<lower=-100,upper=100>  Const ;
}
model {
// local variables
  vector<lower=0>[N]  alphas ;   
  vector<lower=0>[N]  betas ;    
  matrix[C,C]         PREC_bias ;
  matrix[C,C]         PREC_resp ;
  vector[C]           b_m_mean ;  # b_bias - b_bias_mean
  vector[C]           r_m_mean ;  # b_resp - b_resp_mean
  row_vector[C]       b_m_mean_t_A ;
  row_vector[C]       r_m_mean_t_A ;

// priors
  #rhos implied Unif(0,1)
  #const implied Unif(-100,100)
  phi      ~ gamma(0.0001, 0.0001) ;
  tau_sq_b  ~ inv_gamma(0.001, 0.001) ;
  tau_sq_r  ~ inv_gamma(0.001, 0.001) ;

// spatial smoothing priors
  b_m_mean      <- b_bias - b_bias_mean  ;
  b_m_mean_t_A  <- rep_vector(0.0, C)' ;
  r_m_mean      <- b_resp - b_resp_mean ;
  r_m_mean_t_A  <- rep_vector(0.0, C)' ;
    
  for (i in 1:A_N) {
    b_m_mean_t_A[A1[i]] <- (b_m_mean_t_A[A1[i]] + b_m_mean[A2[i]]) ;
    b_m_mean_t_A[A2[i]] <- (b_m_mean_t_A[A2[i]] + b_m_mean[A1[i]]) ;
    r_m_mean_t_A[A1[i]] <- (r_m_mean_t_A[A1[i]] + r_m_mean[A2[i]]) ;
    r_m_mean_t_A[A2[i]] <- (r_m_mean_t_A[A2[i]] + r_m_mean[A1[i]]) ;
  }
    
  // increment log probability manually for multi_normal_prec  
  increment_log_prob(-0.5*tau_sq_b*dot_self(sqrt_d .* b_m_mean) ) ;
  increment_log_prob(-0.5*tau_sq_r*dot_self(sqrt_d .* r_m_mean) ) ;
  increment_log_prob(0.5*tau_sq_b*rho_bias*dot_product(b_m_mean_t_A, b_m_mean)) ;
  increment_log_prob(0.5*tau_sq_r*rho_resp*dot_product(r_m_mean_t_A, r_m_mean)) ;
  PREC_bias <- -rho_bias * A ;
  PREC_resp <- -rho_resp * A ;
  for(c in 1:C) {
    PREC_bias[c,c] <- PREC_bias[c,c] + d[c] ;
    PREC_resp[c,c] <- PREC_resp[c,c] + d[c] ;
  }
  increment_log_prob(0.5 * log_determinant(PREC_bias)) ;
  increment_log_prob(0.5 * log_determinant(PREC_resp)) ;
  increment_log_prob(0.5 * C * log( tau_sq_b )) ;
  increment_log_prob(0.5 * C * log( tau_sq_r)) ;


// likelihood
  for (n in 1:N) {
    real theta ;
	  theta <- inv_logit(Const + b_bias[congress[n]] + b_resp[congress[n]]*lnmajvavg[n]) ;
	  alphas[n] <- theta * phi ;
	  betas[n]  <- (1 - theta) * phi ;
  }

  majps ~ beta_binomial(nvotes, alphas, betas) ;

}
generated quantities {
  vector[N]   log_lik ;
  real        bias[C] ;
  int         y_rep[N] ;

  for (c in 1:C) {
    bias[c] <- inv_logit(b_bias[c]) - 0.5 ;
  }

  for (n in 1:N) {
    real alpha_n ;
    real beta_n ;
    alpha_n <- phi * inv_logit(Const + b_bias[congress[n]] + b_resp[congress[n]]*lnmajvavg[n]) ;
    beta_n  <- phi * (1 - inv_logit(Const + b_bias[congress[n]] + b_resp[congress[n]]*lnmajvavg[n])) ;
    log_lik[n] <- beta_binomial_log(majps[n], nvotes[n], alpha_n, beta_n) ;
    y_rep[n] <- beta_binomial_rng(nvotes[n], alpha_n, beta_n) ;
  }

}
