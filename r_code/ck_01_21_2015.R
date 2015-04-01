
# setwd("/Users/jgabry/Desktop/COLUMBIA/Stuff_for_Wawro/Rsync/cox_katz")
library(rstan)
library(foreign)
# data <- read.dta("cox_katz_replication/data/subgroup-roll-call4606.dta")
data <- read.dta("subgroup-roll-call4606.dta")

period <- with(data, congress - 45)

N <- nrow(data)
C <- length(unique(data$congress))

PREC <- matrix(NA, C, C)
for (i in 1:C) {
  for(j in 1:C) {
    PREC[i,j] <- 0.0;
  }
  PREC[i,i] <- 6.0;
}
PREC[1,1] <- 1.0;
PREC[2,2] <- 5.0;
PREC[C-1,C-1] <- 5.0;
PREC[C,C] <- 1.0;
PREC[1,2] <- -2.0;
PREC[2,1] <- -2.0;
PREC[C,C-1] <- -2.0;
PREC[C-1,C] <- -2.0;
PREC[1,3] <- 1.0;
PREC[3,1] <- 1.0;
PREC[C,C-2] <- 1.0;
PREC[C-2,C] <- 1.0;
for (i in 2:(C - 2)) {
  PREC[i,i+1] <- -4.0;
  PREC[i+1,i] <- -4.0;
  PREC[i,i+2] <- 1.0;
  PREC[i+2,i] <- 1.0;
}

stan_data <- list(N = N,
                  C = C,
                  PREC = PREC,
                  congress = period,
                  nvotes = data$nvotes, 
                  majps = data$majps, 
                  lnmajvavg = data$lnmajvavg)



fit_compile <- stan(file = "stan_tests/01_21_2015/ck_01_21_2015.stan", 
                    data = stan_data, chains = 1, iter = 10)
fit_compile <- stan(file = "01_21_2015/ck_01_21_2015.stan", 
                    data = stan_data, chains = 1, iter = 10)

ck_01_21_2015 <- stan(fit = fit_compile, data = stan_data, refresh = 10)

save(ck_01_21_2015, file = "01_21_2015/ck_01_21_2015.RData")


source("WAIC.R")
WAIC(ck_01_21_2015)

bias <- extract(ck_01_21_2015, pars = "bias")[[1]]
bias_q <- apply(bias, 2, quantile, probs = c(0.025, 0.5, 0.975))

plot(0, xlim = c(1,C), ylim = range(bias_q), pch = 20, type = "n", 
     axes = FALSE, xlab = "Congress", ylab = "Bias")
abline(h = 0, col = "red")
segments(1:C, bias_q[1,], 1:C, bias_q[3,], lwd = 6, col = "lightgray")
lines(1:C, bias_q[2,], ylim = range(bias_q), type = "o", pch = 20)
axis(1, lwd = 3, lwd.ticks = .5)
axis(2, lwd = 3, lwd.ticks = .5)

