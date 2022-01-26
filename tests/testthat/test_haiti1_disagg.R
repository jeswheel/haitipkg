library(pomp)
library(haitipkg)
library(dplyr)
library(testthat)

#### Build Model ###############################################################

## there is only one type: no vaccinations, 10 departments
mod_disagg <- haiti1_disagg()

#### Tests #####################################################################

#### Test 1: check that model has the appropriate times
test_that("check time period", {
  expect_equal(head(mod_disagg@times, n = 1), 1) ## first time in data is 1 week
  expect_equal(tail(mod_disagg@times, n = 1), 430) ## last time in data is 430 weeks
})

#### Test 2: check that model has the appropriate initial time
test_that("check initial time", {
  expect_equal(mod_disagg@t0, 0) ## initial time is 0 weeks
})

#### Test 3: check that models have the appropriate parameters
test_that("check parameters", {
  ## get parameter values
  pars <- unlist(MODEL1_INPUT_PARAMETERS$dep_params)
  depts_names <- c("Artibonite",
                   "Centre",
                   "Grand_Anse",
                   "Nippes",
                   "Nord",
                   "Nord_Est",
                   "Nord_Ouest",
                   "Ouest",
                   "Sud",
                   "Sud_Est")
  rho_epis <- paste0("rho_epi_", depts_names)
  rho_ends <- paste0("rho_end_", depts_names)
  tau_epis <- paste0("tau_epi_", depts_names)
  tau_ends <- paste0("tau_end_", depts_names)
  sig_sq_epis <- paste0("sig_sq_epi_", depts_names)
  sig_sq_ends <- paste0("sig_sq_end_", depts_names)
  beta1s <- paste0("beta1_", depts_names)
  beta2s <- paste0("beta2_", depts_names)
  beta3s <- paste0("beta3_", depts_names)
  beta4s <- paste0("beta4_", depts_names)
  beta5s <- paste0("beta5_", depts_names)
  beta6s <- paste0("beta6_", depts_names)
  nus <- paste0("nu_", depts_names)
  S0s <- paste0("S_", depts_names, "_0")
  E0s <- paste0("E_", depts_names, "_0")
  I0s <- paste0("I_", depts_names, "_0")
  A0s <- paste0("A_", depts_names, "_0")
  R0s <- paste0("R_", depts_names, "_0")
  pop0s <- paste0("pop0_", depts_names)
  mob_cs <- paste0("mob_c_", depts_names)
  par_names <- c(rho_epis, rho_ends, tau_epis, tau_ends, sig_sq_epis,
                 sig_sq_ends, beta1s, beta2s, beta3s, beta4s, beta5s, beta6s,
                 nus, S0s, E0s, I0s, A0s, R0s, pop0s, mob_cs,
                 "gamma", "sigma", "theta0", "alpha", "mu", "delta", "kappa")
  rho_epis <- pars[grepl('rho_epi', names(pars))]
  rho_ends <- pars[grepl('rho_end', names(pars))]
  tau_epis <- pars[grepl('tau_epi', names(pars))]
  tau_ends <- pars[grepl('tau_end', names(pars))]
  sig_sq_epis <- pars[grepl('sig_sq_epi', names(pars))]
  sig_sq_ends <- pars[grepl('sig_sq_end', names(pars))]
  beta1s <- pars[grepl('beta1', names(pars))]
  beta2s <- pars[grepl('beta2', names(pars))]
  beta3s <- pars[grepl('beta3', names(pars))]
  beta4s <- pars[grepl('beta4', names(pars))]
  beta5s <- pars[grepl('beta5', names(pars))]
  beta6s <- pars[grepl('beta6', names(pars))]
  nus <- pars[grepl('nu', names(pars))]
  S0s <- pars[grepl('S_0', names(pars))]
  E0s <- pars[grepl('E_0', names(pars))]
  I0s <- pars[grepl('I_0', names(pars))]
  A0s <- pars[grepl('A_0', names(pars))]
  R0s <- pars[grepl('R_0', names(pars))]
  pop0s <- pars[grepl('pop_0', names(pars))]
  mob_cs <- pars[grepl('mob_c', names(pars))]
  par_vals <- c(rho_epis, rho_ends, tau_epis, tau_ends, sig_sq_epis,
                sig_sq_ends, beta1s, beta2s, beta3s, beta4s, beta5s, beta6s,
                nus, S0s, E0s, I0s, A0s, R0s, pop0s, mob_cs,
                3.5, 5.0, 0.0, 0.00239726, 0.0004287149, 0.000143317, 0.0) %>%
    unlist()
  names(par_vals) <- par_names

  expect_identical(mod_disagg@params, par_vals)
})
