#' Fit Model 3 SpatPomp
#'
#' This function fits the spatPomp version of Model 3 to cholera incidence
#' data from Oct 2010 - Jan 2019 in Haiti. The data are measured at a
#' weekly timescale.
#'
#' @param nsearches integer number of searches to conduct. See details below.
#' @param search[i] list containing parameters used to fit the model. See details.
#' @param ncores Number of cores used to fit the model. The code is written
#'    so that the optimal number of cores with `RUN_LEVEL = 3` is 36.
#'
#' @import foreach
#' @import doParallel
#' @import doRNG
#' @import spatPomp
#'
#' @importFrom magrittr %>%
#'
#' @export
fit_coupled_haiti3 <- function(
    search1 = list(
      NBPF = 5,
      NP = 50,
      SPAT_REGRESSION = 0.5,
      NREPS = 6,
      NP_EVAL = 100,
      NREPS_EVAL = 6
    ),
    search2 = list(
      TOP_N = 2,
      NBPF = 5,
      NP = 50,
      SPAT_REGRESSION = 0.5,
      NREPS = 3,
      NP_EVAL = 100,
      NREPS_EVAL = 6
    ),
    search3 = list(
      TOP_N = 2,
      NBPF = 5,
      NP = 100,
      SPAT_REGRESSION = 0.5,
      NREPS = 3,
      NP_EVAL = 100,
      NREPS_EVAL = 6
    ),
    ncores = 3,
    nsearches = 2
    ) {

  #
  ##
  ### Search1: Global search
  ##
  #

  # Create the model that will be fit to cholera incidence data
  h3_spat <- haiti3_spatPomp()

  # Create a list to save all of the results.
  results <- list()

  # Create vectors for the unit and shared parameters
  unit_specific_names <- c("betaB", "foi_add")
  shared_param_names <- c(
    "mu_B", "XthetaA", "thetaI", "lambdaR", "r", "std_W",
    "epsilon", "k", "sigma"
  )
  est_param_names <- c(
    unit_specific_names, shared_param_names
  )

  # Add unit numbers to each parameter
  est_param_names_expanded <- paste0(rep(est_param_names, each = 10), 1:10)

  # Create rw.sd for each parameter, for search 1
  reg_rw_1.sd <- lapply(est_param_names_expanded, function(x) 0.01)
  names(reg_rw_1.sd) <- est_param_names_expanded
  chol_rw_1.sd <- do.call(rw.sd, reg_rw_1.sd)

  # Create rw.sd for each parameter, for search 2
  reg_rw_2.sd <- lapply(est_param_names_expanded, function(x) 0.00175)
  names(reg_rw_2.sd) <- est_param_names_expanded
  chol_rw_2.sd <- do.call(rw.sd, reg_rw_2.sd)

  # Create rw.sd for each parameter, for search 3
  reg_rw_3.sd <- lapply(est_param_names_expanded, function(x) 0.001)
  names(reg_rw_3.sd) <- est_param_names_expanded
  chol_rw_3.sd <- do.call(rw.sd, reg_rw_3.sd)

  # Get lower bound for unit parameters (global search)
  min_val <- 1e-8
  unit_lb <- rep(c(min_val, min_val), each = 10)
  names(unit_lb) <- paste0(rep(unit_specific_names, each = 10), 1:10)

  # Get upper bound for unit parameters (global search)
  unit_ub <- rep(c(50, 1e-5), each = 10)
  names(unit_ub) <- paste0(rep(unit_specific_names, each = 10), 1:10)

  # Get lower bound for shared parameters (global search)
  shared_lb <- rep(c(5, min_val, min_val, min_val,
                     min_val, min_val, 0.25, 5, 0.01))
  names(shared_lb) <- shared_param_names

  # Get upper bound for shared parameters (global search)
  shared_ub <- c(300, 1, 2e-3, 5, 1.2, 0.15, 1, 1000, 0.5)
  names(shared_ub) <- shared_param_names

  # Create data.frame with random unit parameters
  guesses_unit <- pomp::runif_design(
    lower = unit_lb,
    upper = unit_ub,
    nseq = search1$NREPS
  )

  # Create data.frame with random shared parameters
  guesses_shared <- pomp::runif_design(
    lower = shared_lb,
    upper = shared_ub,
    nseq = search1$NREPS
  )

  # Need to duplicate each of the shared parameter columns
  guesses_shared <- guesses_shared[, rep(1:length(shared_param_names), each = 10)]
  colnames(guesses_shared) <- paste0(rep(shared_param_names, each = 10), 1:10)

  # Combine the unit and shared parameters
  guesses <- cbind(guesses_shared, guesses_unit)

  # We need to add fixed parameters
  all_params <- coef(h3_spat)
  fixed_params <- all_params[!names(all_params) %in% colnames(guesses)]
  fixed_mat <- matrix(
    rep(fixed_params, search1$NREPS_EVAL),
    byrow = TRUE, nrow = search1$NREPS_EVAL
  )
  colnames(fixed_mat) <- names(all_params[!names(all_params) %in% colnames(guesses)])
  guesses_all <- cbind(guesses, fixed_mat)[names(coef(h3_spat))]

  # Memory clean-up
  rm(guesses_unit, guesses_shared, fixed_mat, min_val, shared_lb,
     shared_ub, unit_lb, unit_ub, reg_rw_1.sd, reg_rw_2.sd, reg_rw_3.sd,
     all_params, fixed_params, est_param_names, est_param_names_expanded)
  gc()

  doParallel::registerDoParallel(ncores)
  registerDoRNG(2198635)

  t_global <- system.time(
    foreach(
      i = 1:search1$NREPS,
      .packages = c("spatPomp"),
      .combine = c
    ) %dopar% {
      r_params <- unlist(guesses_all[i, ])
      coef(h3_spat) <- r_params

      ibpf_out <- ibpf(
        h3_spat,
        Nbpf = search1$NBPF,
        Np = search1$NP,
        sharedParNames = shared_param_names,
        unitParNames = unit_specific_names,
        spat_regression = search1$SPAT_REGRESSION,
        rw.sd = chol_rw_1.sd,
        cooling.fraction.50 = 0.35,
        block_size = 1
        # params = r_params  Not yet implemented in spatPomp
      )
    } -> Global_ibpf
  )

  rm(chol_rw_1.sd, guesses_all)

  ibpf_logLik <- data.frame(
    'logLik' = rep(0, length(Global_ibpf)),
    'se' = rep(0, length(Global_ibpf)),
    'which' = 1:length(Global_ibpf)
  )

  t_global_bpf <- system.time(
    for (j in 1:length(Global_ibpf)) {
      ibpf_parms <- coef(Global_ibpf[[j]])

      doRNG::registerDoRNG((j * 38763911) %% 7919)

      coef(h3_spat) <- ibpf_parms

      h3_bpf <- foreach(
        i = 1:search1$NREPS_EVAL, .combine = c, .packages = 'spatPomp'
      ) %dopar% {
        bpfilter(h3_spat, Np = search1$NP_EVAL, block_size = 1,
                 params = ibpf_parms)
      }

      lls <- logLik(h3_bpf)
      ibpf_logLik[j, 1:2] <- logmeanexp(lls[!is.na(lls)], se = TRUE)
    }
  )

  search1_results <- list()
  search1_results$logLiks <- ibpf_logLik
  search1_results$params <- t(coef(Global_ibpf))
  search1_results$ibpf_time <- t_global
  search1_results$bpf_time <- t_global_bpf

  results$search1 <- search1_results

  if (nsearches == 1) {
    return(results)
  }

  rm(ibpf_parms, j, lls, ibpf_logLik, h3_bpf, guesses,
     Global_ibpf, t_global, t_global_bpf)
  gc()

  #
  ##
  ### Search 2: Local Search
  ##
  #

  top_n_global <- order(-search1_results$logLiks$logLik)[1:search2$TOP_N]
  params <- search1_results$params[rep(top_n_global, each = search2$NREPS), ]

  registerDoRNG(327498615)

  t_ibpf_local <- system.time(
    foreach(
      i = 1:(search2$TOP_N * search2$NREPS),
      .packages = c("spatPomp"),
      .combine = c
    ) %dopar% {
      r_params <- params[i, ]
      coef(h3_spat) <- r_params

      ibpf_out <- ibpf(
        h3_spat,
        Nbpf = search2$NBPF,
        Np = search2$NP,
        sharedParNames = shared_param_names,
        unitParNames = unit_specific_names,
        spat_regression = search2$SPAT_REGRESSION,
        rw.sd = chol_rw_2.sd,
        cooling.fraction.50 = 0.5,
        block_size = 1
        # params = r_params  Not yet implemented in spatPomp
      )
    } -> local_ibpf
  )

  rm(chol_rw_2.sd, params)

  ibpf_logLik <- data.frame(
    'logLik' = rep(0, length(local_ibpf)),
    'se' = rep(0, length(local_ibpf)),
    'starting_set' = rep(top_n_global, each = search2$NREPS)
  )

  t_local_bpf <- system.time(
    for (j in 1:length(local_ibpf)) {
      ibpf_parms <- coef(local_ibpf[[j]])

      doRNG::registerDoRNG((j * 38763911) %% 7919)

      coef(h3_spat) <- ibpf_parms

      h3_bpf <- foreach(
        i = 1:search2$NREPS_EVAL, .combine = c, .packages = 'spatPomp'
      ) %dopar% {
        bpfilter(h3_spat, Np = search2$NP_EVAL, block_size = 1,
                 params = ibpf_parms)
      }

      lls <- logLik(h3_bpf)
      ibpf_logLik[j, 1:2] <- logmeanexp(lls[!is.na(lls)], se = TRUE)
    }
  )

  search2_results <- list()
  search2_results$logLiks <- ibpf_logLik
  search2_results$params <- t(coef(local_ibpf))
  search2_results$ibpf_time <- t_ibpf_local
  search2_results$bpf_time <- t_local_bpf

  results$search2 <- search2_results

  if (nsearches == 1) {
    return(results)
  }

  rm(search1, search1_results, ibpf_parms, j, lls, top_n_global, h3_bpf,
     ibpf_logLik, local_ibpf)

  gc()

  #
  ##
  ### Search 3: Local Search
  ##
  #

  top_n_local <- order(-search2_results$logLiks$logLik)[1:search3$TOP_N]
  params <- search2_results$params[rep(top_n_local, each = search3$NREPS), ]

  registerDoRNG(327498615)

  t_ibpf_local <- system.time(
    foreach(
      i = 1:(search3$TOP_N * search3$NREPS),
      .packages = c("spatPomp"),
      .combine = c
    ) %dopar% {
      r_params <- params[i, ]
      coef(h3_spat) <- r_params

      ibpf_out <- ibpf(
        h3_spat,
        Nbpf = search3$NBPF,
        Np = search3$NP,
        sharedParNames = shared_param_names,
        unitParNames = unit_specific_names,
        spat_regression = search3$SPAT_REGRESSION,
        rw.sd = chol_rw_3.sd,
        cooling.fraction.50 = 0.5,
        block_size = 1
        # params = r_params  Not yet implemented in spatPomp
      )
    } -> local_ibpf
  )

  rm(chol_rw_3.sd, params)

  ibpf_logLik <- data.frame(
    'logLik' = rep(0, length(local_ibpf)),
    'se' = rep(0, length(local_ibpf)),
    'starting_set' = rep(top_n_local, each = search3$NREPS)
  )

  t_local_bpf <- system.time(
    for (j in 1:length(local_ibpf)) {
      ibpf_parms <- coef(local_ibpf[[j]])

      doRNG::registerDoRNG((j * 38763911) %% 7919)

      coef(h3_spat) <- ibpf_parms

      h3_bpf <- foreach(
        i = 1:search2$NREPS_EVAL, .combine = c, .packages = 'spatPomp'
      ) %dopar% {
        bpfilter(h3_spat, Np = search2$NP_EVAL, block_size = 1,
                 params = ibpf_parms)
      }

      lls <- logLik(h3_bpf)
      ibpf_logLik[j, 1:2] <- logmeanexp(lls[!is.na(lls)], se = TRUE)
    }
  )

  search3_results <- list()
  search3_results$logLiks <- ibpf_logLik
  search3_results$params <- t(coef(local_ibpf))
  search3_results$ibpf_time <- t_ibpf_local
  search3_results$bpf_time <- t_local_bpf

  results$search3 <- search3_results


  results
}
