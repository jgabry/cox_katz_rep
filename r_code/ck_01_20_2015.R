
# setwd("/Users/jgabry/Desktop/COLUMBIA/Stuff_for_Wawro/Rsync/cox_katz")
library(rstan)
library(foreign)
# data <- read.dta("cox_katz_replication/data/subgroup-roll-call4606.dta")
load("data/ck_data.RData")
data <- ck_data

period <- with(data, congress - 45)

N <- nrow(data)
C <- length(unique(data$congress))

make_A_matrix <- function(dim, type = c("breslow", "rw1", "rw2", "rw3")){
  require(Matrix)
  require(spdep)
  if (type == "breslow") {
    diag <- c(1,5, rep(6, dim-4), 5, 1)
    band1 <- c(-2, rep(-4, dim-3), -2)
    band2 <- rep(1, dim-2)
    Diags <- list(diag, band1, band2)
    k <- c(0,1,2)
    R <- bandSparse(n = dim, k = k, diag = Diags, symm=TRUE)
    return(as.matrix(R))
  }
  if (type == "rw1") {
    nb <- cell2nb(1, dim)
  }
  if (type == "rw2") {
    nb <- cell2nb(1, dim)
    nb[[1]] <- c(nb[[1]], 3)
    nb[[2]] <- c(nb[[2]], 4)
    for(i in 3:(dim-2)) {
      nb[[i]] <- c(i - 2, nb[[i]], i + 2)
    }
    nb[[dim-1]] <- c(nb[[dim-1]], dim-3)
    nb[[dim]] <- c(nb[[dim]], dim-2) 
  }
  if (type == "rw3") {
    nb <- cell2nb(1, dim)
    nb[[1]] <- c(nb[[1]], 3, 4)
    nb[[2]] <- c(nb[[2]], 4, 5)
    nb[[3]] <- c(1, nb[[3]], 5,6)
    for(i in 4:(dim-3)) {
      nb[[i]] <- c(i - 3, i - 2, nb[[i]], i + 2, i + 3)
    }
    nb[[dim-2]] <- c(dim - 5, dim - 4, nb[[dim-2]], dim)
    nb[[dim-1]] <- c(dim - 4, dim - 3, nb[[dim-1]])
    nb[[dim]] <- c(dim - 3, dim-2, nb[[dim]]) 
  }
  A <- matrix(NA, dim,dim)
  for(i in 1:dim) {
    for(j in 1:dim) {
      A[i,j] <- ifelse(i %in% nb[[j]], 1, 0)
    }
  }
  return(A)
}


library(Matrix)

A <- make_A_matrix(C, type = "rw1")
D <- diag(rowSums(A))

PREC <- D - 0.99*A
ld_PREC <- log(det(PREC))

SIGMA <- solve(PREC)

stan_data <- list(N = N,
                  C = C,
                  SIGMA = SIGMA,
                  congress = period,
                  nvotes = data$nvotes, 
                  majps = data$majps, 
                  lnmajvavg = data$lnmajvavg)


# inits <-   list(
#   tau_b_noise = 1.0, 
#   tau_r_noise = 1.0, 
#   phi = 10,
#   noise = array(rep(runif(C, -.01, .01), 2), dim = c(2,C)),
# )
fit_compile <- stan(file = "01_20_2015/ck_01_20_2015.stan", 
                    data = stan_data, chains = 1, iter = 10)

jk_01_20_2015 <- stan(fit = fit_compile, data = stan_data, refresh = 10)
save(jk_01_20_2015, file = "01_20_2015/jk_01_20_2015.RData")


source("WAIC.R")
WAIC(ck_01_20_2015)

bias <- extract(ck_01_20_2015, pars = "bias")[[1]]
bias_q <- apply(bias, 2, quantile, probs = c(0.025, 0.5, 0.975))

plot(0, xlim = c(1,C), ylim = range(bias_q), pch = 20, type = "n", 
     axes = FALSE, xlab = "Congress", ylab = "Bias")
abline(h = 0, col = "red")
segments(1:C, bias_q[1,], 1:C, bias_q[3,], lwd = 6, col = "lightgray")
lines(1:C, bias_q[2,], ylim = range(bias_q), type = "o", pch = 20)
axis(1, lwd = 3, lwd.ticks = .5)
axis(2, lwd = 3, lwd.ticks = .5)

