library(pomp)
library(haitipkg)
library(dplyr)

#### Build Models ##############################################################
#### Goal: test an assortment of possible model builds/configurations

## default model: no vaccinations
novac <- haiti1_joint()

## 3 department vaccination campaigns
vac3 <- haiti1_joint(vacscen = "id4")

## national vaccination campaigns
vacnat <- haiti1_joint(vacscen = "id3")


#### Tests #####################################################################

#### Test 1: check that models have the appropriate subset of data according to
####         time period
test_that("check time period", {
  ## no vacs
  expect_equal(head(novac@times, n = 1), 1) ## first time in data is 1 week
  expect_equal(tail(novac@times, n = 1), 430) ## last time in data is 430 weeks
  ## 3 dept vacs, epi
  expect_equal(head(vac3@times, n = 1), 1) ## first time in data is 1 week
  expect_equal(tail(vac3@times, n = 1), 430) ## last time in data is 430 weeks
  ## nat vacs, end
  expect_equal(head(vacnat@times, n = 1), 1) ## first time in data is 1 week
  expect_equal(tail(vacnat@times, n = 1), 430) ## last time in data is 430
})

#### Test 2: check that models have the appropriate initial time
test_that("check initial time", {
  ## no vacs
  expect_equal(novac@t0, 0) ## initial time is 0 weeks
  ## 3 dept vacs
  expect_equal(vac3@t0, 0) ## initial time is 0 weeks
  ## nat vacs
  expect_equal(vacnat@t0, 0) ## initial time is 0 weeks
})

#### Test 3: check that models have the appropriate parameters
test_that("check parameters", {
  novac_params <- c("rho_epi", "rho_end", "tau_epi", "tau_end", "sig_sq_epi",
                    "sig_sq_end", "beta1", "beta2", "beta3", "beta4", "beta5",
                    "beta6", "nu", "gamma", "sigma", "theta0", "alpha", "mu",
                    "delta", "kappa", "S_0", "E_0", "I_0", "A_0", "R_0", "pop_0")
  vac3_params <- c("rho_epi", "rho_end", "tau_epi", "tau_end", "sig_sq_epi",
                   "sig_sq_end", "beta1", "beta2", "beta3", "beta4", "beta5",
                   "beta6", "nu", "gamma", "sigma", "theta0", "alpha", "mu",
                   "delta", "kappa", "S_0", "E_0", "I_0", "A_0", "R_0", "pop_0",
                   paste0("S", 1:3, "_0"),
                   paste0("E", 1:3, "_0"),
                   paste0("I", 1:3, "_0"),
                   paste0("A", 1:3, "_0"),
                   paste0("R", 1:3, "_0"))
  vacnat_params <- c("rho_epi", "rho_end", "tau_epi", "tau_end", "sig_sq_epi",
                     "sig_sq_end", "beta1", "beta2", "beta3", "beta4", "beta5",
                     "beta6", "nu", "gamma", "sigma", "theta0", "alpha", "mu",
                     "delta", "kappa", "S_0", "E_0", "I_0", "A_0", "R_0", "pop_0",
                     paste0("S", 1:10, "_0"),
                     paste0("E", 1:10, "_0"),
                     paste0("I", 1:10, "_0"),
                     paste0("A", 1:10, "_0"),
                     paste0("R", 1:10, "_0"))
  ## no vacs
  expect_equal(names(novac@params), novac_params)
  ## 3 dept vacs
  expect_equal(names(vac3@params), vac3_params)
  ## nat vacs
  expect_equal(names(vacnat@params), vacnat_params)
})
