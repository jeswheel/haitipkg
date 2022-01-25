library(pomp)
library(haitipkg)
library(dplyr)

#### Build Models ##############################################################
#### Goal: test an assortment of possible model builds/configurations

## default model: no vaccinations, epidemic period
novac_epi <- haiti1()

## no vaccinations, endemic period
novac_end <- haiti1(period = "endemic")

## 3 department vaccination campaigns, epidemic period
vac3_epi <- haiti1(vacscen = "id4", period = "epidemic")

## national vaccination campaigns, endemic period
vacnat_end <- haiti1(vacscen = "id3", period = "endemic")


#### Tests #####################################################################

#### Test 1: check that models have the appropriate subset of data according to
####         time period
test_that("check time period", {
  ## no vacs, epi
  expect_equal(head(novac_epi@times, n = 1), 1) ## first time in data is 1 week
  expect_equal(tail(novac_epi@times, n = 1), 232) ## last time in data is 232 weeks
  ## no vacs, end
  expect_equal(head(novac_end@times, n = 1), 233) ## first time in data is 233 weeks
  expect_equal(tail(novac_end@times, n = 1), 430) ## last time in data is 430
  ## 3 dept vacs, epi
  expect_equal(head(vac3_epi@times, n = 1), 1) ## first time in data is 1 week
  expect_equal(tail(vac3_epi@times, n = 1), 232) ## last time in data is 232 weeks
  ## nat vacs, end
  expect_equal(head(vacnat_end@times, n = 1), 233) ## first time in data is 233 weeks
  expect_equal(tail(vacnat_end@times, n = 1), 430) ## last time in data is 430
})

#### Test 2: check that models have the appropriate initial time
test_that("check initial time", {
  ## no vacs, epi
  expect_equal(novac_epi@t0, 0) ## initial time is 0 weeks
  ## no vacs, end
  expect_equal(novac_end@t0, 232) ## initial time is 232 weeks
  ## 3 dept vacs, epi
  expect_equal(vac3_epi@t0, 0) ## initial time is 0 weeks
  ## nat vacs, end
  expect_equal(vacnat_end@t0, 232) ## initial time is 232 weeks
})

#### Test 3: check that models have the appropriate parameters
test_that("check parameters", {
  novac_params <- c("rho", "tau", "beta1", "beta2", "beta3", "beta4", "beta5",
                    "beta6", "nu", "gamma", "sigma", "theta0", "alpha", "mu",
                    "delta", "sig_sq", "S_0", "E_0", "I_0", "A_0", "R_0", "pop_0")
  vac3_params <- c("rho", "tau", "beta1", "beta2", "beta3", "beta4", "beta5",
                   "beta6", "nu", "gamma", "sigma", "theta0", "alpha", "mu", "delta",
                   "sig_sq", "S_0","E_0","I_0","A_0","R_0", "pop_0", "kappa",
                   paste0("S", 1:3, "_0"),
                   paste0("E", 1:3, "_0"),
                   paste0("I", 1:3, "_0"),
                   paste0("A", 1:3, "_0"),
                   paste0("R", 1:3, "_0"))
  vacnat_params <- c("rho", "tau", "beta1", "beta2", "beta3", "beta4", "beta5",
                     "beta6", "nu", "gamma", "sigma", "theta0", "alpha", "mu", "delta",
                     "sig_sq", "S_0","E_0","I_0","A_0","R_0", "pop_0", "kappa",
                     paste0("S", 1:10, "_0"),
                     paste0("E", 1:10, "_0"),
                     paste0("I", 1:10, "_0"),
                     paste0("A", 1:10, "_0"),
                     paste0("R", 1:10, "_0"))
  ## no vacs, epi
  expect_equal(names(novac_epi@params), novac_params)
  ## no vacs, end
  expect_equal(names(novac_end@params), novac_params)
  ## 3 dept vacs, epi
  expect_equal(names(vac3_epi@params), vac3_params)
  ## nat vacs, end
  expect_equal(names(vacnat_end@params), vacnat_params)
})
