source("ck_pp_check.R")
source("WAIC.R")
library(rstan)
library(foreign)
load("data/ck_data.RData")
data <- ck_data

load("stan_tests/01_20_2015/ck_01_20_2015.RData")

ck_pp_check(stanfit = ck_01_15_2015, data = data)
ck_pp_check(stanfit = ck_01_19_2015, data = data)
ck_pp_check(stanfit = ck_01_20_2015, data = data)

WAIC(ck_01_15_2015)
WAIC(ck_01_19_2015)
WAIC(ck_01_20_2015)




