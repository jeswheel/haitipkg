#' Fit Model 3
#'
#' This function fits the panelPOMP version of Model 3 to cholera incidence
#' data from Oct 2010 - Jan 2019 in Haiti. The data are measured at a
#' weekly timescale.
#'
#' Because of the complex structure of the model, and the high-dimensional nature
#' of panelPomp models, Model 3 requires several rounds of Panel Iterated
#' Filtering in order to obtain the best possible fits. This function breaks the
#' model fitting process into 6 different phases:
#' \itemize{
#'    \item{Phase 1: }{Global search of all model parameters.}
#'    \item{Phase 2: }{Local search of all unit specific parameters, with starting
#'    values derived from the top \emph{TOP_N} results of the global search.}
#'    \item{Phases 3-4: }{Local search of all unit-specific parameters, starting
#'    from the top \emph{TOP_N} results of the previous local search.}
#'    \item{Phase 5: }{Local search for all parameters.}
#' }
#'
#' The user has some control over each of these phases. Importantly, by setting
#' \emph{n_searches}, a user may make the fitting process shorter by avoiding
#' some of the later phases. Additional MIF2 search parameters that can be set
#' for each phase are the following:
#' \itemize{
#'    \item{TOP_N: }{Number of results from the previous search to use as a
#'    starting point for the next search. Note that this argument is ignored for
#'    the first global search, and for searches after local unit searches.}
#'    \item{NP: }{Number of particles to use in the Panel Iterated Filtering.}
#'    \item{NMIF: }{Number of MIF iterations to use in the Panel Iterated Filtering.}
#'    \item{NREPS: }{For each starting point, how many replicated searches should
#'    be conducted. Here, we recommend that NREPS * TOP_N be a multiple of
#'    ncores.}
#'    \item{NREPS_EVAL: }{Number of times to replicate the evaluation using a
#'    particle filter.}
#'    \item{NP_EVAL: }{Number of particles to use in model evaluation.}
#'    \item{RW_SD: }{Specification of the rw.sd used in the MIF2 algorithm.}
#' }
#'
#' @param n_searches integer number of searches to conduct. See details below.
#' @param search[i] list containing parameters used to fit the model. See details.
#' @param ncores Number of cores used to fit the model. The code is written
#'    so that the optimal number of cores with `RUN_LEVEL = 3` is 36.
#'
#' @import pomp
#' @import panelPomp
#' @import foreach
#' @import doRNG
#'
#' @importFrom magrittr %>%
#' @importFrom foreach %dopar%
#'
#' @export
fit_haiti3 <- function(
    n_searches = 3L,
    search1 = list(
      NP = 20,
      NMIF = 3,
      NREPS = 3,
      NREPS_EVAL = 3,
      NP_EVAL = 25
    ),
    search2 = list(
      TOP_N = 1,
      NP = 20,
      NMIF = 3,
      NREPS = 3,
      NREPS_EVAL = 3,
      NP_EVAL = 25
    ),
    search3 = list(
      TOP_N = 1,
      NP = 20,
      NMIF = 3,
      NREPS = 3,
      NREPS_EVAL = 3,
      NP_EVAL = 25
    ),
    ncores = 3
) {

  # Set of names required for each search. Throw error if missing
  req_names <- c('NP', 'NMIF', 'NREPS', 'NREPS_EVAL', 'NP_EVAL')

  # Throw error if missing names in search1
  if (any(!req_names %in% names(search1))) {
    missing <- req_names[which(!req_names %in% names(search1))]
    stop(
      paste0(
        "search1 missing required argument(s): ",
        paste0(missing, collapse = ', '),
        "."
      )
    )
  }

  # Throw error if missing names in search2
  if (n_searches >= 2 && any(!req_names %in% names(search2))) {
    missing <- req_names[which(!req_names %in% names(search2))]
    stop(
      paste0(
        "search2 missing required argument(s): ",
        paste0(missing, collapse = ', '),
        "."
      )
    )
  }

  # Throw error if missing names in search3
  if (n_searches >= 3 && any(!req_names %in% names(search3))) {
    missing <- req_names[which(!req_names %in% names(search3))]
    stop(
      paste0(
        "search3 missing required argument(s): ",
        paste0(missing, collapse = ', '),
        "."
      )
    )
  }

  # TOP_N is not required, but if missing, default to 1.
  if (!"TOP_N" %in% names(search2)) {
    search2$TOP_N <- 1
  }

  if (!"TOP_N" %in% names(search3)) {
    search4$TOP_N <- 1
  }

  # rw.sd is not required, but if not supplied, use these defualts:
  if (!"RW_SD" %in% names(search1)) {
    # Set random walk standard deviation
    chol_rw1 <- pomp::rw.sd(
      betaB = 0.02,  # Unit Specific
      mu_B = 0.02,   # Unit Specific
      thetaI = 0.02,
      XthetaA = 0.02,
      lambdaR = 0.02,
      r = 0.01,
      std_W = 0.02,
      epsilon = 0.02,
      k = 0.02,
      sigma = 0.02,
      foi_add = 0.02
    )
  } else {
    chol_rw1 <- search1$RW_SD
  }

  if (!"RW_SD" %in% names(search2)) {
    chol_rw2 <- pomp::rw.sd(
      betaB = 0.0075,  # Unit Specific
      mu_B = 0.0075,   # Unit Specific
      thetaI = 0.0075,
      XthetaA = 0.0075,
      lambdaR = 0.0075,
      r = 0.005,
      std_W = 0.0075,
      epsilon = 0.0075,
      k = 0.0075,
      sigma = 0.0075,
      foi_add = 0.0075
    )
  } else {
    chol_rw2 <- search2$RW_SD
  }

  if (n_searches >= 3 && !"RW_SD" %in% names(search3)) {
    # Set rw sd
    chol_rw3 <- pomp::rw.sd(
      betaB = 0.005,
      foi_add = 0.005
    )
  } else {
    chol_rw3 <- search3[['RW_SD']]
  }

  #### Load model
  SIRB_panel <- haiti3_panel(
    start_time = "2010-10-23",
    B0 = TRUE
  )

  #### Register cores
  doParallel::registerDoParallel(ncores)

  #### Create a list to store final results
  results = list()

  ####
  #### Start search 1: Global search of all parameters
  ####

  min_param_val <- 5e-8  # Min param value for positive parameters

  # Create data.frame of unit-specific parameter bounds
  unit_bounds <- data.frame(
    "param" = c('betaB', 'foi_add'),
    "lower" = c(min_param_val, min_param_val),
    "upper" = c(50, 1e-5)
  )

  # Get the default values saved in Model 3, to use as a starting point
  original_unit <- SIRB_panel@specific
  fixed_unit <- SIRB_panel@specific[c('H', 'D', 'B0'), ]
  shared_params <- SIRB_panel@shared

  # Get names of all of the departments
  deps <- colnames(original_unit)

  # From the table above, create lower and upper bounds for each parameter
  lb_unit <- unlist(unit_bounds[, 'lower'])
  names(lb_unit) <- unlist(unit_bounds[, 'param'], use.names = FALSE)
  ub_unit <- unlist(unit_bounds[, 'upper'])
  names(ub_unit) <- unlist(unit_bounds[, 'param'], use.names = FALSE)

  # Create random starting points
  set.seed(3178689)
  guesses_unit <- pomp::runif_design(
    lower = lb_unit,
    upper = ub_unit,
    nseq = (search1$NREPS - 1) * 10
  )

  # Create list that will contain unit-specific and shared parameters
  guess_list_unit <- list()
  for (i in 1:(search1$NREPS - 1)) {
    Betas <- guesses_unit$betaB[(10 * i - 9):(10 * i)]
    Fois <- guesses_unit$foi_add[(10 * i - 9):(10 * i)]
    params <- rbind(Betas, Fois)
    colnames(params) <- deps
    rownames(params) <- c('betaB', 'foi_add')
    params <- rbind(params, fixed_unit)
    guess_list_unit[[i]] <- params
  }

  # Save final set of parameters with the default unit specific params
  guess_list_unit[[search1$NREPS]] <- original_unit

  # Create data.frame object for the shared parameter bounds
  parameter_bounds <- data.frame(
    "param" = c("mu_B", "XthetaA", "thetaI", "lambdaR", "r",
                "std_W", "epsilon", "k", "sigma"),
    "lower" = c(5, min_param_val, min_param_val, min_param_val,
                min_param_val, min_param_val, 0.25, 5, 0.01),
    "upper" = c(300, 1, 2e-3, 5, 1.2, 0.15, 1, 1000, 0.5)
  )

  # Save all default values for starting parameters later
  original_fixed_shared <- SIRB_panel@shared[names(SIRB_panel@shared) %in% parameter_bounds$param]
  fixed_shared <- SIRB_panel@shared[!names(SIRB_panel@shared) %in% parameter_bounds$param]
  fixed_shared['cas_def'] <- 1  # No change in reporting rate

  # From the table above, create lower and upper bounds for each parameter
  lb <- unlist(parameter_bounds[, 'lower'], use.names = FALSE)
  names(lb) <- unlist(parameter_bounds[, 'param'], use.names = FALSE)
  ub <- unlist(parameter_bounds[, 'upper'], use.names = FALSE)
  names(ub) <- unlist(parameter_bounds[, 'param'], use.names = FALSE)

  # Using the bounds defined above, create a grid of parameters
  # to search globally.
  set.seed(7869381)
  guesses <- pomp::runif_design(
    lower = lb,
    upper = ub,
    nseq = search1$NREPS - 1
  )

  # Set last set of shared parameters to model defaults.
  guesses <- rbind(guesses, original_fixed_shared[colnames(guesses)])

  # set RNG for reproducible parallelization
  doRNG::registerDoRNG(1851563)

  # Loop through each starting point and perform MIF2 search.
  foreach(
    i=1:search1$NREPS,
    .packages = c('panelPomp'),
    .combine = c
  ) %dopar% {
    r_shared_params <- unlist(guesses[i, ])
    r_unit_params <- guess_list_unit[[i]]
    mif2(
      SIRB_panel,
      Np = search1$NP,
      Nmif = search1$NMIF,
      cooling.fraction.50 = 0.5,
      rw.sd = chol_rw1,
      cooling.type = 'geometric',
      shared.start = c(fixed_shared, r_shared_params),
      specific.start = r_unit_params,
      block = TRUE
    )
  } -> no_trend_global

  # Some clutter clean-up
  rm(
    chol_rw1, fixed_unit, guess_list_unit, guesses, guesses_unit, original_unit,
    parameter_bounds, params, unit_bounds, Betas, fixed_shared, Fois, i, lb,
    lb_unit, original_fixed_shared, shared_params, ub, ub_unit
  )

  # Garbage collector to free unused memory after clean-up.
  gc()

  #### pfilter global search results

  mif_logLik <- matrix(nrow = length(no_trend_global), ncol = 10)
  colnames(mif_logLik) <- names(unitobjects(SIRB_panel))
  mif_logLik <- as.data.frame(mif_logLik)
  # mif_logLik$which <- 1:length(no_trend_global)

  # Loop through PIF results and perform particle filters to evaluate likelihood
  for (j in 1:length(no_trend_global)) {
    mf <- no_trend_global[[j]]
    mif_params <- coef(mf)

    doRNG::registerDoRNG((j * 38763911) %% 7919)

    pf3_loglik_matrix <- foreach(i=1:search1$NREPS_EVAL, .combine = rbind) %dopar% {
      library(panelPomp)
      unitlogLik(pfilter(mf, params = mif_params, Np = search1$NP_EVAL))
    }

    mif_logLik[j, 1:10] <- apply(pf3_loglik_matrix, 2, logmeanexp)
  }

  # Save results from the global search
  search1_results <- list()
  search1_results$logLiks <- mif_logLik
  search1_results$params  <- t(sapply(no_trend_global, coef))

  # Save search results in final output list
  results$search1 <- search1_results

  # clean-up and garbage collection
  rm(mf, mif_logLik, no_trend_global, pf3_loglik_matrix, j, mif_params, search1)
  gc()

  if (n_searches == 1L) {
    return(results)
  }

  ####
  #### End of global search
  ####

  #### Start local search of top half of the searches

  top_n_global <- order(-apply(search1_results$logLiks, 1, sum))[1:search2$TOP_N]
  params <- search1_results$params[rep(top_n_global, each = search2$NREPS), ]

  registerDoRNG(987153547)

  foreach(
    i = 1:(search2$NREPS * search2$TOP_N),
    .packages = c('panelPomp'),
    .combine = c
  ) %dopar% {

    start_params <- params[i, ]

    mif2(
      SIRB_panel,
      Np = search2$NP,
      Nmif = search2$NMIF,
      cooling.fraction.50 = 0.5,
      rw.sd = chol_rw2,
      cooling.type = 'geometric',
      start = start_params
    )
  } -> local_MIF2_search

  rm(
    chol_rw2, params
  )

  gc()

  mif_logLik <- matrix(nrow = search2$TOP_N * search2$NREPS, ncol = 10)
  colnames(mif_logLik) <- names(unitobjects(SIRB_panel))
  mif_logLik <- as.data.frame(mif_logLik)
  mif_logLik$starting_set <- rep(top_n_global, each = search2$NREPS)

  for (j in 1:length(local_MIF2_search)) {
    mf <- local_MIF2_search[[j]]
    mif_params <- coef(mf)

    # registerDoParallel(36)
    registerDoRNG((j * 38763911) %% 7919)

    pf3_loglik_matrix <- foreach(i=1:search2$NREPS_EVAL, .combine = rbind,
                                 .packages = 'panelPomp') %dopar% {
                                   unitlogLik(pfilter(mf, params = mif_params, Np = search2$NP_EVAL))
                                 }

    mif_logLik[j, 1:10] <- apply(pf3_loglik_matrix, 2, logmeanexp)
  }

  # Save results from the search
  search2_results <- list()
  search2_results$logLiks <- mif_logLik
  search2_results$params  <- t(sapply(local_MIF2_search, coef))

  results$search2 <- search2_results

  rm(j, mf, mif_params, local_MIF2_search, mif_logLik, search2, pf3_loglik_matrix,
     top_n_global, search1_results)
  gc()

  if (n_searches == 2L) {
    return(results)
  }

  #### Start local search of unit specific parameters

  top_n_local <- order(-apply(search2_results$logLiks, 1, sum))[1:search3$TOP_N]
  params <- search2_results$params[rep(top_n_local, each = search3$NREPS), ]

  registerDoRNG(987153547)

  foreach(
    i = 1:(search3$NREPS * search3$TOP_N),
    .packages = c('panelPomp'),
    .combine = c
  ) %dopar% {

    start_params <- params[i, ]

    mif2(
      SIRB_panel,
      Np = search3$NP,
      Nmif = search3$NMIF,
      cooling.fraction.50 = 0.5,
      rw.sd = chol_rw3,
      cooling.type = 'geometric',
      start = start_params
    )
  } -> local_MIF2_search

  rm(
    chol_rw3, params
  )

  gc()

  mif_logLik <- matrix(nrow = length(local_MIF2_search), ncol = 10)
  colnames(mif_logLik) <- names(unitobjects(SIRB_panel))
  mif_logLik <- as.data.frame(mif_logLik)
  mif_logLik$starting_set <- rep(top_n_local, each = search3$NREPS)

  for (j in 1:length(local_MIF2_search)) {
    mf <- local_MIF2_search[[j]]
    mif_params <- coef(mf)

    # registerDoParallel(36)
    registerDoRNG((j * 38763911) %% 7919)

    pf3_loglik_matrix <- foreach(i=1:search3$NREPS_EVAL, .combine = rbind,
                                 .packages = 'panelPomp') %dopar% {
      unitlogLik(pfilter(mf, params = mif_params, Np = search3$NP_EVAL))
    }

    mif_logLik[j, 1:10] <- apply(pf3_loglik_matrix, 2, logmeanexp)
  }

  # Save results from the search
  search3_results <- list()
  search3_results$logLiks <- mif_logLik
  search3_results$params  <- t(sapply(local_MIF2_search, coef))

  # Save the maximum log-likes for each unit within each grouped shared parameter
  search3_results$group_logLiks <- mif_logLik %>%
    dplyr::group_by(starting_set) %>%
    dplyr::summarize(dplyr::across(Artibonite:`Sud-Est`, max, na.rm = TRUE))

  # Record which set of parameters within a group resulted in best unit parameter
  search3_results$which_logLiks <- mif_logLik %>%
    dplyr::group_by(starting_set) %>%
    dplyr::summarize(dplyr::across(Artibonite:`Sud-Est`, which.max))

  # For convenience, save the shared starting values (did not change after this search)
  shared_params_finish <- t(sapply(
    unique(search3_results$which_logLiks$starting_set),
    function(x) {
      panelPomp::pParams(results$search2$params[x, ])$shared
    }
  ))

  search3_results$shared_params_finish <- shared_params_finish

  # Get matrix of all unit parameters after doing MIF
  unit_params_finish <- lapply(
    1:nrow(search3_results$params),
    function(x) {
      panelPomp::pParams(search3_results$params[x, ])$specific
    }
  )

  best_unit_combo <- list()
  best_ll <- numeric(nrow(search3_results$which_logLiks))

  # Loop through each group (starting shared parameters)
  for (i in 1:nrow(search3_results$which_logLiks)) {

    # Find which likelihoods are the largest for group g
    which_unit_max <- unlist(search3_results$which_logLiks[i, -1])

    # Create a list of only the unit specific parameters for this particular group
    group_unit_params <- unit_params_finish[((i - 1) * search3$NREPS + 1):(i * search3$NREPS)]

    # Create a matrix of unit specific parameters where the parameters for unit u
    # are those that maximized the likelihood within the group.
    best_pp <- sapply(1:10, function(col) group_unit_params[[which_unit_max[col]]][, col])
    colnames(best_pp) <- colnames(group_unit_params[[1]])

    best_unit_combo[[i]] <- best_pp
    best_ll[i] <- sum(search3_results$group_logLiks[i, -1])
  }

  search3_results$best_unit_combo <- best_unit_combo
  search3_results$best_ll_combo <- best_ll

  results$search3 <- search3_results

  rm(j, mf, mif_params, local_MIF2_search, mif_logLik, search3, pf3_loglik_matrix,
     top_n_local, search2_results, best_pp, best_unit_combo, group_unit_params,
     unit_params_finish, best_ll, i, which_unit_max, shared_params_finish)
  gc()

  if (n_searches == 3L) {
    return(results)
  }

  results
}
