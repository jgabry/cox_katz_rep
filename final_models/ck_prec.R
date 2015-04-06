library(rstan)
library(foreign)

# Load data
# data <- read.dta("cox_katz_replication/data/subgroup-roll-call4606.dta")
data <- read.dta("subgroup-roll-call4606.dta")

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

# Adjacency matrix
A <- make_A_matrix(dim = C, type = "rw1")

# Degree matrix
D <- diag(rowSums(A))

stan_data <- list(N = N,
                  C = C,
                  A = A,
                  D = D,
                  congress = period,
                  nvotes = data$nvotes, 
                  majps = data$majps, 
                  lnmajvavg = data$lnmajvavg)


fit_compile <- stan(file = "final_models/ck_prec.stan", 
                    data = stan_data, chains = 1, iter = 10)


# save_pars <- fit_compile@model_pars[-grep("noise", fit_compile@model_pars)]

ck_prec <- stan(fit = fit_compile, 
                data = stan_data, 
                iter = 500, 
                chains = 6,
                # pars = save_pars,
                refresh = 10)

save(ck_prec, file = "final_models/ck_prec.RData", compress = "xz")


majps <- data$majps
library(shinyStan)
launch_shinystan(ck_cholesky)

