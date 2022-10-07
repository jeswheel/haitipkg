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
      NREPS_EVAL = 6,
      RW_SD = NULL,
      COOLING = 0.5,
      KEEP_TRACES = FALSE
    ),
    search2 = list(
      TOP_N = 2,
      NBPF = 5,
      NP = 50,
      SPAT_REGRESSION = 0.5,
      NREPS = 3,
      NP_EVAL = 100,
      NREPS_EVAL = 6,
      RW_SD = NULL,
      COOLING = 0.5,
      KEEP_TRACES = FALSE
    ),
    search3 = list(
      TOP_N = 2,
      NBPF = 5,
      NP = 100,
      SPAT_REGRESSION = 0.5,
      NREPS = 3,
      NP_EVAL = 100,
      NREPS_EVAL = 6,
      RW_SD = NULL,
      COOLING = 0.5,
      KEEP_TRACES = FALSE
    ),
    ncores = 3,
    nsearches = 2,
    gamma = 182.625
    ) {

  #
  ##
  ### Search1: Global search
  ##
  #

  if (is.null(search1$COOLING)) {
    search1$COOLING <- 0.5
  }

  if (is.null(search2$COOLING)) {
    search2$COOLING <- 0.5
  }

  if (is.null(search3$COOLING)) {
    search3$COOLING <- 0.5
  }

  # Create the model that will be fit to cholera incidence data
  h3_spat <- haiti3_spatPomp()

  h3_spat@params[paste0('gamma', 1:10)] <- gamma

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

  if (is.null(search1$RW_SD)) {
    # Create rw.sd for each parameter, for search 1
    reg_rw_1.sd <- lapply(est_param_names_expanded, function(x) 0.01)
    names(reg_rw_1.sd) <- est_param_names_expanded
    chol_rw_1.sd <- do.call(rw.sd, reg_rw_1.sd)
    rm(reg_rw_1.sd)
  } else {
    chol_rw_1.sd <- search1$RW_SD
  }


  if (is.null(search2$RW_SD)) {
    # Create rw.sd for each parameter, for search 2
    reg_rw_2.sd <- lapply(est_param_names_expanded, function(x) 0.00175)
    names(reg_rw_2.sd) <- est_param_names_expanded
    chol_rw_2.sd <- do.call(rw.sd, reg_rw_2.sd)
    rm(reg_rw_2.sd)
  } else {
    chol_rw_2.sd <- search2$RW_SD
  }


  if (is.null(search3$RW_SD)) {
    # Create rw.sd for each parameter, for search 3
    reg_rw_3.sd <- lapply(est_param_names_expanded, function(x) 0.001)
    names(reg_rw_3.sd) <- est_param_names_expanded
    chol_rw_3.sd <- do.call(rw.sd, reg_rw_3.sd)
    rm(reg_rw_3.sd)
  } else {
    chol_rw_3.sd <- search3$RW_SD
  }


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
     shared_ub, unit_lb, unit_ub,  all_params, fixed_params,
     est_param_names)
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
        cooling.fraction.50 = search1$COOLING,
        block_size = 1
        # params = r_params  Not yet implemented in spatPomp
      )
    } -> Global_ibpf
  )

  rm(guesses_all)

  pfilterLikes <- data.frame(
    "ll" = rep(0, search1$NREPS_EVAL*search1$NREPS),
    "which" = rep(1:search1$NREPS, each = search1$NREPS_EVAL)
  )

  t_global_bpf <- system.time(
    likMat <- foreach(  # Get log-likelihood for each unit and set of parameters, NREPS_EVAL times each
      i = 1:(search1$NREPS_EVAL*search1$NREPS), .combine = rbind, .packages = 'spatPomp'
    ) %dopar% {
      p3 <- coef(Global_ibpf[[(i-1) %/% search1$NREPS_EVAL + 1]])
      coef(h3_spat) <- p3
      apply(bpfilter(
        h3_spat, Np = search1$NP_EVAL,
        block_size = 1, parms = p3
      )@block.cond.loglik, 1, sum)
    }
  )

  # Condense unit likelihoods into model likelihood
  pfilterLikes$ll <- apply(likMat, 1, sum)

  # Group by model parameter set.
  ibpf_logLik <- pfilterLikes %>%
    dplyr::group_by(which) %>%
    dplyr::summarize(logLik = logmeanexp(ll),
              se = logmeanexp(ll, se = TRUE)[2])

  search1_results <- list()
  search1_results$logLiks <- ibpf_logLik
  search1_results$params <- t(coef(Global_ibpf))
  search1_results$ibpf_time <- t_global
  search1_results$bpf_time <- t_global_bpf
  search1_results$likMat <- likMat

  if (isTRUE(search1$KEEP_TRACES)) {
    search1_results$traces <- lapply(Global_ibpf, function(x) x@traces[, c('loglik', names(chol_rw_1.sd@call)[-1])])
  }

  results$search1 <- search1_results

  if (nsearches == 1) {
    return(results)
  }

  rm(ibpf_logLik, guesses, likMat, pfilterLikes,
     Global_ibpf, chol_rw_1.sd, t_global, t_global_bpf)
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
        cooling.fraction.50 = search2$COOLING,
        block_size = 1
        # params = r_params  Not yet implemented in spatPomp
      )
    } -> local_ibpf
  )

  rm(params)

  pfilterLikes <- data.frame(
    "ll" = 0,
    "starting_set" = rep(top_n_global, each = search2$NREPS * search2$NREPS_EVAL),
    "which" = rep(1:(search2$NREPS * search2$TOP_N), each = search1$NREPS_EVAL)
  )

  t_local_bpf <- system.time(
    likMat <- foreach(  # Get log-likelihood for each unit and set of parameters, NREPS_EVAL times each
      i = 1:(search2$NREPS_EVAL*search2$NREPS*search2$TOP_N), .combine = rbind, .packages = 'spatPomp'
    ) %dopar% {
      p3 <- coef(local_ibpf[[(i-1) %/% search2$NREPS_EVAL + 1]])
      coef(h3_spat) <- p3
      apply(bpfilter(
        h3_spat, Np = search2$NP_EVAL,
        block_size = 1, parms = p3
      )@block.cond.loglik, 1, sum)
    }
  )

  # Condense unit likelihoods into model likelihood
  pfilterLikes$ll <- apply(likMat, 1, sum)

  # Group by model parameter set.
  ibpf_logLik_temp <- pfilterLikes %>%
    dplyr::group_by(which) %>%
    dplyr::summarize(logLik = logmeanexp(ll),
                     se = logmeanexp(ll, se = TRUE)[2])

  ibpf_logLik <- pfilterLikes %>%
    dplyr::group_by(which) %>%
    dplyr::summarise(starting_set = dplyr::first(starting_set)) %>%
    dplyr::left_join(
      x = ibpf_logLik_temp,
      y = .,
      by = "which"
    ) %>%
    dplyr::select(which, starting_set, logLik, se)

  search2_results <- list()
  search2_results$logLiks <- ibpf_logLik
  search2_results$params <- t(coef(local_ibpf))
  search2_results$ibpf_time <- t_ibpf_local
  search2_results$bpf_time <- t_local_bpf
  search2_results$likMat <- likMat

  if (isTRUE(search2$KEEP_TRACES)) {
    search2_results$traces <- lapply(local_ibpf, function(x) x@traces[, c('loglik', names(chol_rw_2.sd@call)[-1])])
  }

  results$search2 <- search2_results

  if (nsearches == 2) {
    return(results)
  }

  rm(
    search1, search1_results, top_n_global, t_ibpf_local, t_local_bpf,
    pfilterLikes, ibpf_logLik, chol_rw_2.sd, local_ibpf, likMat, ibpf_logLik_temp
  )

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
        cooling.fraction.50 = search3$COOLING,
        block_size = 1
        # params = r_params  Not yet implemented in spatPomp
      )
    } -> local_ibpf
  )

  rm(params)

  pfilterLikes <- data.frame(
    "ll" = 0,
    "starting_set" = rep(top_n_local, each = search3$NREPS * search3$NREPS_EVAL),
    "which" = rep(1:(search3$NREPS * search3$TOP_N), each = search3$NREPS_EVAL)
  )

  t_local_bpf <- system.time(
    likMat <- foreach(  # Get log-likelihood for each unit and set of parameters, NREPS_EVAL times each
      i = 1:(search3$NREPS_EVAL*search3$NREPS*search3$TOP_N), .combine = rbind, .packages = 'spatPomp'
    ) %dopar% {
      p3 <- coef(local_ibpf[[(i-1) %/% search3$NREPS_EVAL + 1]])
      coef(h3_spat) <- p3
      apply(bpfilter(
        h3_spat, Np = search3$NP_EVAL,
        block_size = 1, parms = p3
      )@block.cond.loglik, 1, sum)
    }
  )

  # Condense unit likelihoods into model likelihood
  pfilterLikes$ll <- apply(likMat, 1, sum)

  # Group by model parameter set.
  ibpf_logLik_temp <- pfilterLikes %>%
    dplyr::group_by(which) %>%
    dplyr::summarize(logLik = logmeanexp(ll),
                     se = logmeanexp(ll, se = TRUE)[2])

  ibpf_logLik <- pfilterLikes %>%
    dplyr::group_by(which) %>%
    dplyr::summarise(starting_set = dplyr::first(starting_set)) %>%
    dplyr::left_join(
      x = ibpf_logLik_temp,
      y = .,
      by = "which"
    ) %>%
    dplyr::select(which, starting_set, logLik, se)

  search3_results <- list()
  search3_results$logLiks <- ibpf_logLik
  search3_results$params <- t(coef(local_ibpf))
  search3_results$ibpf_time <- t_ibpf_local
  search3_results$bpf_time <- t_local_bpf
  search3_results$likMat <- likMat

  if (isTRUE(search3$KEEP_TRACES)) {
    search3_results$traces <- lapply(local_ibpf, function(x) x@traces[, c('loglik', c('loglik', names(chol_rw_3.sd@call)[-1]))])
  }

  results$search3 <- search3_results

  results
}
