#' Fit Model 1
#'
#' This function fits Model 1 by performing calculating a Monte-Carlo adjusted
#' profile for the linear trend in transmission parameter.
#'
#' @param NP Number of particles to use for the Iterated Filter
#' @param NMIF Number of MIF iterations
#' @param NUM_TREND Number of points on profile to try for the trend parameter.
#' @param NPROF Number of global search for each value of the trend parameter.
#' @param NREPS_EVAL Number of particle filter replications to evaluate the
#'        likelihood.
#' @param NP_EVAL Number of particles used for likelihood evaluation.
#' @param ncores Number of cores to use for parallelization.
#'
#' @examples
#' fit_m1 <- fit_haiti1()
#'
#' @import pomp
#' @import foreach
#' @import doRNG
#'
#' @importFrom magrittr %>%
#'
#' @export
fit_haiti1 <- function(
    NP = 50, NMIF = 3, NUM_TREND = 3, NPROF = 3,
    NREPS_EVAL = 3, NP_EVAL = 50, ncores = 3
    # rho_flag = TRUE, tau_flag = TRUE, sig_sq_flag = TRUE,
    # beta_flag = FALSE, nu_flag = FALSE
) {

  doParallel::registerDoParallel(ncores)
  set.seed(636813)

  # Load the model, allow breakpoint in rho, tau, and sig_sq, but not beta or nu
  mod_1 <- haiti1_joint(
    rho_flag = FALSE,
    tau_flag = TRUE,
    sig_sq_flag = TRUE,
    beta_flag = FALSE,
    nu_flag = FALSE
  )

  # Minimum value for positive parameters
  min_param_val <- 5e-8

  # Create table for bounds used for starting values in global search
  bounds <- tibble::tribble(
    ~param, ~lower, ~upper,
    "beta1",         .75, 1.75,
    "beta2",         .75, 1.75,
    "beta3",         .75, 1.75,
    "beta4",         .75, 1.75,
    "beta5",         .75, 1.75,
    "beta6",         .75, 1.75,
    "tau_epi",       100,  850,
    "tau_end",       100,  850,
    "rho",          0.15,    1,
    "nu",           0.95,    1,
    "sig_sq_epi",   0.05, 0.15,
    "sig_sq_end",   0.05, 0.15,
    "E_0", min_param_val, 2e-3,
    "I_0", min_param_val, 3e-3
  )

  lower <- bounds$lower
  names(lower) <- bounds$param

  upper <- bounds$upper
  names(upper) <- bounds$param

  # Create starting value for the profile search
  guesses <- profile_design(
    betat = seq(-0.15, 0.05, length = NUM_TREND),
    lower = lower,
    upper = upper,
    nprof = NPROF,
    type = "runif"
  )

  # Default the endemic starting values to the same as epidemic values
  guesses$tau_end <- guesses$tau_epi
  guesses$rho_end <- guesses$rho_epi
  guesses$sig_sq_end <- guesses$sig_sq_epi

  fixed_params <- coef(mod_1)[!names(coef(mod_1)) %in% c(colnames(guesses), "S_0")]

  # Set rw.sd
  rw_sd <- rw.sd(
    beta1 = .02,
    beta2 = .02,
    beta3 = .02,
    beta4 = .02,
    beta5 = .02,
    beta6 = .02,
    rho   = 0.02,
    nu    = 0.01,
    E_0   = ivp(0.2),
    I_0   = ivp(0.2),
    tau_epi    = ifelse(time < 232, 0.02, 0),
    tau_end    = ifelse(time >= 232, 0.02, 0),
    sig_sq_epi = ifelse(time < 232, 0.02, 0),
    sig_sq_end = ifelse(time >= 232, 0.02, 0)
  )

  # Make reproducible results
  registerDoRNG(17495987)

  # Conduct profile search
  foreach(
    i = 1:nrow(guesses),
    .packages = c('pomp'),
    .combine = c
  ) %dopar% {

    r_params <- unlist(guesses[i, ])
    S_0 <- unname(1 - r_params['E_0'] - r_params["I_0"])
    names(S_0) <- "S_0"
    coef(mod_1) <- c(r_params, fixed_params, S_0)[names(coef(mod_1))]
    mif2(
      mod_1,
      Np = NP,
      Nmif = NMIF,
      cooling.fraction.50 = 0.5,  # Cooling set so that we will reach smaller rw.sd at the end than the end of the unit3 search.
      rw.sd = rw_sd
    ) -> m2
  } -> MIF2_search

  # Just removing everything so that there is no issue with renaming stuff,
  # and so that we are being a bit more memory efficient.
  rm(NMIF, NP, lower, fixed_params, min_param_val,
     upper, rw_sd, bounds)
  gc()

  # pfilter the results -----------------------------------------------------

  ll_matrix <- matrix(
    nrow = NPROF * NUM_TREND, ncol = NREPS_EVAL
  )

  for (j in 1:length(MIF2_search)) {

    mif_params <- coef(MIF2_search[[j]])
    registerDoRNG((j * 687383921) %% 7919)

    ll_evals <- foreach(i=1:NREPS_EVAL, .combine = c, .packages = 'pomp') %dopar% {
      logLik(pfilter(mod_1, params = mif_params, Np = NP_EVAL))
    }

    ll_matrix[j, ] <- ll_evals
  }

  lls <- as.data.frame(t(apply(ll_matrix, 1, logmeanexp, se = TRUE)))
  results <- as.data.frame(t(sapply(MIF2_search, coef)))
  results$ll <- lls$V1
  results$ll.se <- lls$se

  # See model1/haiti1_VaccinationScenarios.Rmd for further processing
  results
}
