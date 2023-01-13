#### Build Models ##############################################################
#### Goal: test an assortment of possible model builds/configurations

## default model: no vaccinations
novac <- haiti1_joint(
  rho_flag = FALSE,
  tau_flag = TRUE,
  sig_sq_flag = TRUE,
  beta_flag = FALSE,
  nu_flag = FALSE
)

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
