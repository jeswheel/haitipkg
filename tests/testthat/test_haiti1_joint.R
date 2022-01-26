library(pomp)
library(haitipkg)
library(dplyr)
library(testthat)

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
  novac_params <- unlist(MODEL1_INPUT_PARAMETERS$joint_pars)
  par_names <- names(novac_params)
  depts <- 3
  vac3_params <- c(novac_params, rep(0.0, 5 * depts))
  names(vac3_params) <- c(par_names,
                          paste0("S", 1:depts, "_0"),
                          paste0("E", 1:depts, "_0"),
                          paste0("I", 1:depts, "_0"),
                          paste0("A", 1:depts, "_0"),
                          paste0("R", 1:depts, "_0"))
  depts <- 10
  vacnat_params <- c(novac_params, rep(0.0, 5 * depts))
  names(vacnat_params) <- c(par_names,
                            paste0("S", 1:depts, "_0"),
                            paste0("E", 1:depts, "_0"),
                            paste0("I", 1:depts, "_0"),
                            paste0("A", 1:depts, "_0"),
                            paste0("R", 1:depts, "_0"))
  ## no vacs
  expect_identical(novac@params, novac_params)
  ## 3 dept vacs
  expect_identical(vac3@params, vac3_params)
  ## nat vacs
  expect_identical(vacnat@params, vacnat_params)
})
