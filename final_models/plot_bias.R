library(ggplot2)
library(gridExtra)
library(plyr)
library(rstan)

load("ck_data.RData")
source("final_models/plot_code.R")
prep_data <- function(fit, ck_data, cred_lev = 0.95) {
  if (!inherits(fit, "stanfit")) {
    quants <- fit
  } else {
    bias_posterior <- rstan::extract(fit, pars = "bias")[[1]]
    alpha <- 1 - (1 - cred_lev)/2
    quants <- apply(bias_posterior, 2, quantile, probs = c(1-alpha, 0.5, alpha))
  }
  
  congress_data <- plyr::ddply(ck_data, "congress", summarise, majority = mean(demmaj))
  Majority <- factor(congress_data$majority, labels = c("Republican", "Democrat"))
  data.frame(Congress = congress_data$congress, Majority, LB = quants[1,], Bias = quants[2, ], UB = quants[3,])
}


load("final_models/ck_cholesky_reformat.RData")
df <- prep_data(fit = ck_cholesky_reformat, ck_data = ck_data, cred_lev = 0.95)
gmrfplot <- make_plot()




# make same graph for cox and katz results
source("ben_apsa/APSA2012_dump.R")
library(foreign)
library(rstan)
library(mvtnorm)

RC <- read.dta("ben_apsa/subgroup-roll-call4606.dta")

Congresses <- unique(RC$congress)
out <- sapply(Congresses, simplify = FALSE, FUN = function(C) {
  model <- glm(cbind(majps, nvotes - majps) ~ lnmajvavg, data = RC, family = binomial(),
               subset = congress >= C - 3 & congress <= C + 3)
  betas <- rmvnorm(10000, mean = coef(model), sigma = vcov(model))
  betas[,1] <- plogis(betas[,1]) - 0.5
  return(betas)
})

qs <- mat.or.vec(nr = 3, nc = length(Congresses))
for(i in seq_along(out)) {
  qs[,i] <- quantile(out[[i]][,"(Intercept)"], probs = c(.025, 0.5, .975))
}

df <- prep_data(fit = qs, ck_data = ck_data)
ckplot <- make_plot()

pdf("final_models/ck_replication.pdf", w = 8, h = 10)
grid_arrange_shared_legend(ckplot, gmrfplot)
graphics.off()





