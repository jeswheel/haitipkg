#' Fit Model 3
#'
#' This function fits model 3 to the available haiti cholera data
#'
#' As of Feb 16, 2022, Model 3 can be written in 3 different forms:
#' \itemize{
#'    \item Coupled `pomp` object (`haiti3_correct()`).
#'    \item Independent `panelPomp` object (`haiti3_panel()`).
#'    \item Coupled `spatPomp` object (`haiti3_spatPomp()`).
#' }
#'
#' Despite being represented in so many different forms, the model is only fit
#' using `panelPomp` methods. This is because it best replicates what was done
#' in Lee et al. by not using `spatPomp` (which was not available at the time).
#' The modeling fitting is improved from the fitting done by Lee et al by using
#' `panelPomp`, which came out just before the Lee et al. paper did.
#'
#' Because of the complex nature of the `panelPomp` model, we found that
#' performing the model fitting in three different stages (1. Performing a
#' global search of all model parameters, 2. Fix the unit-specific parameters
#' found in the global search, and perform another global search just on the
#' shared parameters, 3. Fix the shared parameters found in the global search
#' and perform a local search on the unit-specific parameters). Note that adding
#' more additional stages or and changing any of the three stages we used may
#' have resulted in a slightly different estimate.
#'
#'  The following is a more in depth description of the `RUN_LEVEL` argument:
#' \describe{
#'    \item{`RUN_LEVEL = 1`}{Stage 1: `Np = 50` and `Nmif = 3` for global search,
#'    with 3 starting locations for the parameters used. Stage 2: `Np = 50` and
#'    `Nmif = 3` and 3 starting locations for global search. Stage 3: `Np = 50`
#'    and `Nmif = 3`, repeating the local search 3 times. At the end of each
#'    stage, 50 particles are used to estimate the likelihood, and 3 replications
#'    of the estimation are conducted.}
#'    \item{`RUN_LEVEL = 2`}{Stage 1: `Np = 400` and `Nmif = 15` for global search,
#'    with 10 starting locations for the parameters used. Stage 2: `Np = 400` and
#'    `Nmif = 15` and 10 starting locations for global search. Stage 3: `Np = 400`
#'    and `Nmif = 15`, repeating the local search 3 times. At the end of each
#'    stage, 400 particles are used to estimate the likelihood, and 5 replications
#'    of the estimation are conducted.}
#'    \item{`RUN_LEVEL = 3`}{Stage 1: `Np = 1000` and `Nmif = 50` for global search,
#'    with 72 starting locations for the parameters used. Stage 2: `Np = 1000` and
#'    `Nmif = 50` and 72 starting locations for global search. Stage 3: `Np = 1500`
#'    and `Nmif = 50`, repeating the local search 36 times. At the end of each
#'    stage, 36 replications are conducted with `Np = 2000, 2500, 3000` at stages
#'    1, 2, and 3, respectively}
#' }
#'
#' @param RUN_LEVEL parameter in `c(1, 2, 3)`. The different RUN_LEVELS
#'    correspond to the computational effort that is used to fit the model, with
#'    `RUN_LEVEL = 1` being the least ammount of computation used to fit and
#'    `RUN_LEVEL = 3` being the largest ammount of computation used to fit the
#'    model. Note that `RUN_LEVEL = 3` is used to obtain our final results, and
#'    `RUN_LEVEL = 1` is primarily used for debugging purposes. More precise
#'    information about each of the run levels is provided in the Details
#'    Section below.
#' @param ncores Number of cores used to fit the model. The code is written
#'    so that the optimal number of cores with `RUN_LEVEL = 3` is 36.
#'
#' @import pomp
#' @import panelPomp
#' @import foreach
#' @import doRNG
#'
#' @importFrom magrittr %>%
#' @importFrom foreach %do%
#' @importFrom foreach foreach
#'
#' @export
fit_haiti3 <- function(RUN_LEVEL = 1, ncores = 1) {


  doParallel::registerDoParallel(ncores)

  NP_GLOBAL    <- switch(RUN_LEVEL, 50, 4e3,  1000)
  NP_LOCAL     <- switch(RUN_LEVEL, 50, 4e3,  1500)
  NMIF         <- switch(RUN_LEVEL,  3,  15,    50)
  NREPS_GLOBAL <- switch(RUN_LEVEL,  3,  10,    72)
  NREPS_LOCAL  <- switch(RUN_LEVEL,  3,  10,    36)
  NREPS_EVAL   <- switch(RUN_LEVEL,  3,   5, cores)
  NP_EVAL1     <- switch(RUN_LEVEL, 50, 4e3,  2000)
  NP_EVAL2     <- switch(RUN_LEVEL, 50, 4e3,  2500)
  NP_EVAL3     <- switch(RUN_LEVEL, 50, 4e3,  3000)

  # STEP 1: Global all ------------------------------------------------------

  SIRB_panel <- haiti3_panel(start_time = "2010-10-23",
                             B0 = TRUE)

  # Smallest value positive parameters will start at
  min_param_val <- 5e-8

  unit_bounds <- tibble::tribble(
    ~param, ~lower, ~upper,
    "betaB", min_param_val, 20,
    "foi_add", min_param_val, 1e-5,
    "B0", 0.15, 0.4
  )

  original_unit <- SIRB_panel@specific
  fixed_unit <- SIRB_panel@specific[c('H', 'D'), ]
  shared_params <- SIRB_panel@shared

  deps <- colnames(original_unit)

  # From the table above, create lower and upper bounds for each parameter
  lb_unit <- unlist(unit_bounds[, 'lower'])
  names(lb_unit) <- unlist(unit_bounds[, 'param'], use.names = FALSE)
  ub_unit <- unlist(unit_bounds[, 'upper'])
  names(ub_unit) <- unlist(unit_bounds[, 'param'], use.names = FALSE)

  set.seed(3178689)
  guesses_unit <- runif_design(
    lower = lb_unit,
    upper = ub_unit,
    nseq = (NREPS_GLOBAL - 1) * 10
  )

  guess_list_unit <- list()
  for (i in 1:(NREPS_GLOBAL - 1)) {
    Betas <- guesses_unit$betaB[(10 * i - 9):(10 * i)]
    Fois <- guesses_unit$foi_add[(10 * i - 9):(10 * i)]
    B0s <- guesses_unit$B0[(10 * i - 9):(10 * i)]
    params <- rbind(Betas, Fois, B0s)
    colnames(params) <- deps
    rownames(params) <- c('betaB', 'foi_add', 'B0')
    params <- rbind(params, fixed_unit)
    guess_list_unit[[i]] <- params
  }

  guess_list_unit[[NREPS_GLOBAL]] <- original_unit

  # Table copied from run_mif_haitiOCV.R, but bounds where changed a bit
  parameter_bounds <- tibble::tribble(
    ~param, ~lower, ~upper,
    "mu_B", 25, 250,
    "XthetaA", min_param_val, 0.6,
    "thetaI", min_param_val, 5e-3,
    "lambdaR", min_param_val, 5,
    "r", min_param_val, 2,
    "std_W", min_param_val, 0.15,
    "epsilon", .25, 1,
    "k", 10, 1000,# hard to get negbin like this, sobol in log scale -5 et 4 TODO IF ENABLE: UNCOMMENT ID2314
    "sigma", 0.01, 0.5
  )

  original_fixed_shared <- SIRB_panel@shared[names(SIRB_panel@shared) %in% parameter_bounds$param]
  fixed_shared <- SIRB_panel@shared[!names(SIRB_panel@shared) %in% parameter_bounds$param]
  fixed_shared['cas_def'] <- 1

  # All parameters are fixed that aren't in the table above
  # fixed_params <- sirb_cholera@params[!names(sirb_cholera@params) %in% parameter_bounds$param]

  # From the table above, create lower and upper bounds for each parameter
  lb <- unlist(parameter_bounds[, 'lower'], use.names = FALSE)
  names(lb) <- unlist(parameter_bounds[, 'param'], use.names = FALSE)
  ub <- unlist(parameter_bounds[, 'upper'], use.names = FALSE)
  names(ub) <- unlist(parameter_bounds[, 'param'], use.names = FALSE)

  # Using the bounds defined above, create a grid of parameters
  # to search globally.
  set.seed(7869381)
  guesses <- runif_design(
    lower = lb,
    upper = ub,
    nseq = NREPS_GLOBAL - 1
  )

  guesses <- rbind(guesses, original_fixed_shared[colnames(guesses)])

  chol_rw <- rw.sd(
    betaB = 0.02,
    mu_B = 0.02,
    thetaI = 0.02,
    XthetaA = 0.02,
    lambdaR = 0.02,
    r = 0.02,
    std_W = 0.02,
    epsilon = 0.02,
    k = 0.02,
    sigma = 0.02,
    foi_add = 0.02,
    B0 = ivp(0.1)
  )

  # Get all of the cores available.
  registerDoRNG(1851563)
  cat('\nStarting Global Search...\n')

  # Global MIF at "MLE"

  # Run global MIF chol_Nreps_global times
  no_trend_global <- foreach(
    i=1:NREPS_GLOBAL,
    .packages = c('panelPomp'),
    .combine = c
  ) %dopar% {
    r_shared_params <- unlist(guesses[i, ])
    r_unit_params <- guess_list_unit[[i]]
    mif2(
      SIRB_panel,
      Np = NP_GLOBAL,
      Nmif = NMIF,
      cooling.fraction.50 = 0.5,
      rw.sd = chol_rw,
      cooling.type = 'geometric',
      shared.start = c(fixed_shared, r_shared_params),
      specific.start = r_unit_params,
      block = TRUE
    )
  }

  cat('Finished!\n\n')

  # Just removing everything so that there is no issue with renaming stuff,
  # and so that we are being a bit more memory efficient.
  rm(
    chol_rw, fixed_unit, guess_list_unit, guesses, guesses_unit, original_unit,
    parameter_bounds, params, unit_bounds, Betas, fixed_shared, Fois, i, lb,
    lb_unit, original_fixed_shared, shared_params, ub, ub_unit
  )

  gc()


  # Step 1.5: PFilter results -----------------------------------------------

  cat('\nPFiltering Global All...\n')
  mif_logLik <- data.frame(
    'logLik' = rep(0, length(no_trend_global)),
    'se' = rep(0, length(no_trend_global)),
    'which' = 1:length(no_trend_global)
  )

  mif_results <- foreach(i=1:length(no_trend_global)) %do% {
    mf <- no_trend_global[[i]]
    list(logLik = logLik(mf), params = coef(mf))
  }

  for (j in 1:length(no_trend_global)) {
    mf <- no_trend_global[[j]]
    mif_params <- mif_results[[j]]$params

    # registerDoParallel(36)
    registerDoRNG((j * 38763911) %% 7919)

    library(tictoc)

    tic()
    pf3_loglik_matrix <- foreach(i=1:NREPS_EVAL, .combine = rbind) %dopar% {
      library(panelPomp)
      unitlogLik(pfilter(mf, params = mif_params, Np = NP_EVAL1))
    }
    toc()

    mif_logLik[mif_logLik$which == j, 1:2] <- panel_logmeanexp(pf3_loglik_matrix, MARGIN = 2, se = TRUE)
  }

  cat('Finished!\n\n')

  # save(mif_logLik, file = 'output/global_all_correctPanel_PF.rda')

  rm(
    i, j, mf, mif_params, mif_results, pf3_loglik_matrix
  )

  gc()

  # Step 2: Global Search of Shared Parameters ------------------------------

  cat('Starting Global Search of Shared...\n')

  best_m <- mif_logLik %>%
    arrange(-logLik) %>%
    slice_head(n = 1) %>%
    pull(which)

  params_shared <- no_trend_global[[best_m]]@shared
  params_unit <- no_trend_global[[best_m]]@specific

  rm(no_trend_global, mif_logLik)
  gc()

  # Table copied from run_mif_haitiOCV.R, but bounds where changed a bit
  parameter_bounds <- tribble(
    ~param, ~lower, ~upper,
    "mu_B", 20, 200,
    "XthetaA", min_param_val, 1,
    "thetaI", min_param_val, 5e-04,
    "lambdaR", min_param_val, 10,
    "r", min_param_val, 5e3,
    "std_W", min_param_val, 0.9,
    "epsilon", .25, 1,
    "k", 10, 1000
  )

  original_fixed_shared <- params_shared[names(params_shared) %in% parameter_bounds$param]
  fixed_shared <- params_shared[!names(params_shared) %in% parameter_bounds$param]

  # From the table above, create lower and upper bounds for each parameter
  lb <- unlist(parameter_bounds[, 'lower'])
  names(lb) <- unlist(parameter_bounds[, 'param'], use.names = FALSE)
  ub <- unlist(parameter_bounds[, 'upper'])
  names(ub) <- unlist(parameter_bounds[, 'param'], use.names = FALSE)

  # Using the bounds defined above, create a grid of parameters
  # to search globally.
  set.seed(17839)
  guesses <- runif_design(
    lower = lb,
    upper = ub,
    nseq = NREPS_LOCAL - 1
  )

  guesses <- rbind(guesses, original_fixed_shared)

  chol_rw <- rw.sd(
    mu_B = 0.02,
    thetaI = 0.02,
    XthetaA = 0.02,
    lambdaR = 0.02,
    r = 0.02,
    std_W = 0.02,
    epsilon = 0.02,
    sigma = 0.02,
    k = 0.02
  )

  registerDoRNG(3481569)

  # Local MIF at "MLE"
  # Run local MIF chol_Nreps_local times
  no_trend_shared <- foreach(
    guess=iter(guesses, 'row'),
    .packages = c('panelPomp'),
    .combine = c
  ) %dopar% {
    r_shared_params <- unlist(guess)
    mif2(
      SIRB_panel,
      Np = NP_GLOBAL,
      Nmif = NMIF,
      cooling.fraction.50 = 0.5,
      rw.sd = chol_rw,
      cooling.type = 'geometric',
      shared.start = c(fixed_shared, r_shared_params),
      specific.start = params_unit,
      block = TRUE
    )
  }

  cat("Finished! \n\n")

  # Just removing everything so that there is no issue with renaming stuff,
  # and so that we are being a bit more memory efficient.
  rm(
    chol_rw, guesses, parameter_bounds, params_unit, best_m, fixed_shared, lb,
    ub, params_shared, original_fixed_shared
  )

  gc()


  # Step 2.5: PFilter results -----------------------------------------------

  cat('\nPFiltering Global Shared...\n')
  mif_logLik <- data.frame(
    'logLik' = rep(0, length(no_trend_shared)),
    'se' = rep(0, length(no_trend_shared)),
    'which' = 1:length(no_trend_shared)
  )

  mif_results <- foreach(i=1:length(no_trend_shared)) %do% {
    mf <- no_trend_shared[[i]]
    list(logLik = logLik(mf), params = coef(mf))
  }

  for (j in 1:length(no_trend_shared)) {
    mf <- no_trend_shared[[j]]
    mif_params <- mif_results[[j]]$params

    # registerDoParallel(36)
    registerDoRNG((j * 38763911) %% 7919)

    library(tictoc)

    tic()
    pf3_loglik_matrix <- foreach(i=1:NREPS_EVAL, .combine = rbind) %dopar% {
      library(panelPomp)
      unitlogLik(pfilter(mf, params = mif_params, Np = NP_EVAL2))
    }
    toc()

    mif_logLik[mif_logLik$which == j, 1:2] <- panel_logmeanexp(pf3_loglik_matrix, MARGIN = 2, se = TRUE)
  }

  cat('Finished!\n\n')
  # save(mif_logLik, file = 'output/global_shared_correctPanel_PF.rda')

  rm(
    i, j, mf, mif_params, mif_results, pf3_loglik_matrix
  )

  gc()

  # Step 3: Local Search of Unit Params -------------------------------------

  cat('Starting Local Search of Unit params...\n')

  params_shared <- no_trend_shared[[which.max(mif_logLik$logLik)]]@shared
  params_unit <- no_trend_shared[[which.max(mif_logLik$logLik)]]@specific

  rm(no_trend_shared, mif_logLik)
  gc()

  chol_rw <- rw.sd(
    betaB = 0.02,
    foi_add = 0.02,
    B0 = ivp(0.1)  # Change this to 0.1.
  )

  # Get all of the cores available.
  registerDoRNG(40765101)

  # Local MIF at "MLE"
  # Run local MIF chol_Nreps_local times
  no_trend_local <- foreach(
    i=1:NREPS_LOCAL,
    .packages = c('panelPomp'),
    .combine = c
  ) %dopar% {
    mif2(
      SIRB_panel,
      Np = NP_LOCAL,
      Nmif = NMIF,
      cooling.fraction.50 = 0.5,
      rw.sd = chol_rw,
      cooling.type = 'geometric',
      shared.start = params_shared,
      specific.start = params_unit,
      block = TRUE
    )
  }

  cat("Finished! \n\n")

  # Just removing everything so that there is no issue with renaming stuff,
  # and so that we are being a bit more memory efficient.
  rm(
    chol_rw, params_unit, params_shared
  )

  gc()


  # Step 3.5: PFilter of Local results --------------------------------------

  cat('\nPFiltering Local Unit...\n')
  mif_logLik <- data.frame(
    'logLik' = rep(0, length(no_trend_local)),
    'se' = rep(0, length(no_trend_local)),
    'which' = 1:length(no_trend_local)
  )

  mif_results <- foreach(i=1:length(no_trend_local)) %do% {
    mf <- no_trend_local[[i]]
    list(logLik = logLik(mf), params = coef(mf))
  }

  for (j in 1:length(no_trend_local)) {
    mf <- no_trend_local[[j]]
    mif_params <- mif_results[[j]]$params

    # registerDoParallel(36)
    registerDoRNG((j * 38763911) %% 7919)

    library(tictoc)

    tic()
    pf3_loglik_matrix <- foreach(i=1:NREPS_EVAL, .combine = rbind) %dopar% {
      library(panelPomp)
      unitlogLik(pfilter(mf, params = mif_params, Np = NP_EVAL2))
    }
    toc()

    mif_logLik[mif_logLik$which == j, 1:2] <- panel_logmeanexp(pf3_loglik_matrix, MARGIN = 2, se = TRUE)
  }

  cat('Finished!\n\n')
  # save(mif_logLik, file = 'output/local_unit_correctPanel_PF.rda')

  rm(
    i, j, mf, mif_params, mif_results, pf3_loglik_matrix
  )

  gc()

  best_m <- mif_logLik %>%
    arrange(-logLik) %>%
    slice_head(n = 1) %>%
    pull(which)

  coef(no_trend_local[[best_m]])

}
