#' Fit Model 3 SpatPomp
#'
#' This function fits the spatPomp version of Model 3 to cholera incidence
#' data from Oct 2010 - Jan 2019 in Haiti. The data are measured at a
#' weekly timescale. The function allows a user to fit the model in up to 3
#' stages by setting \emph{nsearches} equal to 3. Any number larger than 3
#' is ignored.
#'
#' Hyper-parameters to the IBPF algorithm can be adjusted as inputs into a list
#' for each search (e.g. search1 = list(TOP_N = 3, ...)). The possible list
#' options are:
#' \itemize{
#'    \item{TOP_N: }{Number of results from the previous search to use as a
#'    starting point for the next search. Note that this argument is ignored for
#'    the first global search, and for searches after local unit searches.}
#'    \item{NBPF: }{Number of IBPF iterations.}
#'    \item{NP: }{Number of particles to use in in the IBPF algorithm.}
#'    \item{SPAT_REGRESSION: }{The regression coefficient for shared parameters}
#'    \item{NREPS: }{For each starting point, how many replicated searches should
#'    be conducted. Here, we recommend that NREPS * TOP_N be a multiple of
#'    ncores.}
#'    \item{NP_EVAL: }{Number of particles to use in model evaluation.}
#'    \item{NREPS_EVAL: }{Number of times to replicate the evaluation using a
#'    particle filter.}
#'    \item{RW_SD: }{Specification of the rw.sd used in the IBPF algorithm.}
#'    \item{COOLING: }{Optional, defaults to 0.5. Specifies the
#'    cooling.fraction.50 for the rw.sd used in the IBPF algorithm.}
#'    \item{KEEP_TRACES: }{Optional, defaults to FALSE. Saves the traces of the
#'    IBPF output. This is a large list, so is primarily saved for debugging and
#'    testing purposes.}
#'    \item{KEEP_LIK_MAT: }{Optional, defaults to FALSE. Saves a likelihood
#'    matrix containing likelihood estimates for each unit. This is a large
#'    matrix, so it is primarily saved for debugging and testing purposes.}
#' }
#'
#' @param nsearches integer number of searches to conduct. See details below.
#' @param search1 list containing parameters used to fit the model. See details.
#' @param search2 list containing parameters used to fit the model. See details.
#' @param search3 list containing parameters used to fit the model. See details.
#' @param search_rho Boolean indicating whether or not rho should be estimated.
#' @param search_gamma Boolean indicating whether or not gamma should be estimated.
#' @param search_hur Boolean indicating whether or not coefficients related to
#'    Hurricane Matthew (2016) should be perturbed in the global search.
#' @param search_Iinit Boolean indicating whether or not units with zero case
#'    counts at time t1 should have the value \eqn{I_u(t0)} estimated.
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
fit_haiti3 <- function(
    search1 = list(
      NBPF = 5,
      NP = 50,
      SPAT_REGRESSION = 0.5,
      NREPS = 6,
      NP_EVAL = 100,
      NREPS_EVAL = 6,
      RW_SD = NULL,
      COOLING = 0.5,
      KEEP_TRACES = FALSE,
      KEEP_LIK_MAT = FALSE
    ),
    search2 = NULL,
    search3 = NULL,
    ncores = 3,
    nsearches = 1,
    search_rho = FALSE,
    search_gamma = FALSE,
    search_hur = FALSE,
    search_Iinit = FALSE
    ) {

  #
  ##
  ### Search1: Global search
  ##
  #

  if (is.null(search1$COOLING)) {
    search1$COOLING <- 0.5
  }

  if (is.null(search1$KEEP_LIK_MAT)) {
    search1$KEEP_LIK_MAT <- TRUE
  }

  if (nsearches >= 2 && is.null(search2$COOLING)) {
    search2$COOLING <- 0.5
  }

  if (nsearches >= 2 && is.null(search2$KEEP_LIK_MAT)) {
    search2$KEEP_LIK_MAT <- TRUE
  }

  if (nsearches >= 3 && is.null(search3$COOLING)) {
    search3$COOLING <- 0.5
  }

  if (nsearches >= 3 && is.null(search3$KEEP_LIK_MAT)) {
    search3$KEEP_LIK_MAT <- TRUE
  }

  # Create the model that will be fit to cholera incidence data
  h3_spat <- haiti3_spatPomp()

  # Create a list to save all of the results.
  results <- list()

  # Create vectors for the unit and shared parameters
  unit_specific_names <- c("betaB", "foi_add", "aHur", "hHur")

  shared_param_names <- c(
    "mu_B", "XthetaA", "thetaI", "lambdaR", "r", "std_W",
    "epsilon", "k"
  )

  if (search_rho) {
    shared_param_names <- c(shared_param_names, "rho")
  }

  if (search_gamma) {
    shared_param_names <- c(shared_param_names, "gamma")
  }

  est_param_names <- c(
    unit_specific_names, shared_param_names
  )

  # Add unit numbers to each parameter
  est_param_names_expanded <- paste0(rep(est_param_names, each = 10), 1:10)

  if (is.null(search1$RW_SD)) {
    stop("RW_SD must be supplied for search 1.")
  } else {
    chol_rw_1.sd <- search1$RW_SD
  }

  if (nsearches >= 2 && is.null(search2$RW_SD)) {
    stop("RW_SD must be supplied for search 2.")
  } else {
    chol_rw_2.sd <- search2$RW_SD
  }

  if (nsearches >= 3 && is.null(search3$RW_SD)) {
    stop("RW_SD must be supplied for search 3.")
  } else {
    chol_rw_3.sd <- search3$RW_SD
  }

  # Get lower bound for unit parameters (global search)
  min_val <- 1e-8
  unit_lb <- rep(c(1, 1e-9, 0, 0), each = 10)
  names(unit_lb) <- paste0(rep(unit_specific_names[1:4], each = 10), 1:10)

  # Get upper bound for unit parameters (global search)
  unit_ub <- rep(c(50, 1e-6, 0, 0), each = 10)
  names(unit_ub) <- paste0(rep(unit_specific_names[1:4], each = 10), 1:10)

  # Set unique upper-bounds for betaB, based on run_level_2 search results.
  unit_ub['betaB1'] <- 15
  unit_ub['betaB2'] <- 80
  unit_ub['betaB3'] <- 75
  unit_ub['betaB4'] <- 70
  unit_ub['betaB10'] <- 35
  unit_ub['betaB5'] <- 15
  unit_ub['betaB6'] <- 80
  unit_ub['betaB7'] <- 20
  unit_ub['betaB8'] <- 5
  unit_ub['betaB9'] <- 30
  unit_ub['foi_add10'] <- 2e-7
  unit_ub['foi_add3'] <- 2e-7
  unit_ub['foi_add4'] <- 2e-7
  unit_ub['foi_add6'] <- 5e-7
  unit_ub['foi_add9'] <- 5e-7

  if (search_hur) {
    unit_ub['aHur3'] <- 50
    unit_ub['aHur9'] <- 25
    unit_ub['hHur3'] <- 120
    unit_ub['hHur9'] <- 120

    unit_lb['aHur3'] <- 0.01
    unit_lb['aHur9'] <- 0.01
    unit_lb['hHur3'] <- 30
    unit_lb['hHur9'] <- 30
  }

  if (search_Iinit) {
    unit_lb['Iinit3'] <- 0.8 / 468301  # H3
    unit_lb['Iinit4'] <- 0.8 / 342525  # H4

    unit_ub['Iinit3']  <- 40 / 468301
    unit_ub['Iinit4']  <- 40 / 342525
  }

  # Get lower bound for shared parameters (global search)
  shared_lb <- c(
    25, 0.02, min_val, min_val,
    0.75, 0.005, 0.5, 1
  )

  if (search_rho) {
    shared_lb <- c(shared_lb, 1/8)
  }

  if (search_gamma) {
    shared_lb <- c(shared_lb, 365.25 / 5)
  }

  names(shared_lb) <- shared_param_names

  # Get upper bound for shared parameters (global search)
  shared_ub <- c(500, 0.8, 0.0001, 5, 1.75, 0.075, 0.99, 200)

  if (search_rho) {
    shared_ub <- c(shared_ub, 1/(0.2))
  }

  if (search_gamma) {
    shared_ub <- c(shared_ub, 365.25 / 2)
  }

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
    rep(fixed_params, search1$NREPS),
    byrow = TRUE, nrow = search1$NREPS
  )
  colnames(fixed_mat) <- names(all_params[!names(all_params) %in% colnames(guesses)])

  # Combine estimated and fixed parameters, and reorder based on original order.
  guesses_all <- cbind(guesses, fixed_mat)[, names(coef(h3_spat))]

  # Memory clean-up
  rm(guesses_unit, guesses_shared, fixed_mat, min_val, shared_lb,
     shared_ub, unit_lb, unit_ub,  all_params, fixed_params,
     est_param_names)
  gc()

  doParallel::registerDoParallel(ncores)
  doRNG::registerDoRNG(2198635)

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

  ibpf_params <- t(coef(Global_ibpf))

  for (i in 1:nrow(ibpf_params)) {
    for (sp in shared_param_names) {
      ibpf_params[i, paste0(sp, 1:10)] <- mean(ibpf_params[i, paste0(sp, 1:10)])
    }
  }

  t_global_bpf <- system.time(
    likMat <- foreach(  # Get log-likelihood for each unit and set of parameters, NREPS_EVAL times each
      i = 1:(search1$NREPS_EVAL*search1$NREPS), .combine = rbind, .packages = 'spatPomp'
    ) %dopar% {
      p3 <- ibpf_params[(i-1) %/% search1$NREPS_EVAL + 1, ]

      coef(h3_spat) <- p3
      apply(bpfilter(
        h3_spat, Np = search1$NP_EVAL,
        block_size = 1
      )@block.cond.loglik, 1, sum)
    }
  )

  # Condense unit likelihoods into model likelihood
  # pfilterLikes$ll <- apply(likMat, 1, sum)
  pfilterLikes$ll <- rowSums(likMat)  # Faster

  # Group by model parameter set.
  ibpf_logLik <- pfilterLikes %>%
    dplyr::group_by(which) %>%
    dplyr::summarize(logLik = logmeanexp(ll),
              se = logmeanexp(ll, se = TRUE)[2])

  search1_results <- list()
  search1_results$logLiks <- ibpf_logLik
  search1_results$params <- ibpf_params
  search1_results$ibpf_time <- t_global
  search1_results$bpf_time <- t_global_bpf

  if (isTRUE(search1$KEEP_LIK_MAT)) {
    search1_results$likMat <- likMat
  }

  if (isTRUE(search1$KEEP_TRACES)) {
    search1_results$traces <- lapply(Global_ibpf, function(x) x@traces[, c('loglik', names(chol_rw_1.sd@call)[-1])])
  }

  results$search1 <- search1_results

  if (nsearches == 1) {
    return(results)
  }

  rm(ibpf_logLik, guesses, likMat, pfilterLikes,
     Global_ibpf, chol_rw_1.sd, t_global, t_global_bpf,
     i, sp, ibpf_params)
  gc()

  #
  ##
  ### Search 2: Local Search
  ##
  #

  top_n_global <- order(-search1_results$logLiks$logLik)[1:search2$TOP_N]

  # Contains the best set of parameters, measure for the coupled model
  params <- search1_results$params[rep(top_n_global, each = search2$NREPS), ]

  # This is a set of parameters where the unit parameters come from the set
  # of parameters that are best for the unit in question, and the shared
  # parameters are averaged based on the best for each unit (i.e., if
  # there are 2 units and parameter set 5 is the best for unit 1 and parameter
  # set 3 is best for unit 2, the unit-specific parameters for unit 1 are
  # set to the unit specific parameters from search 5, and the shared parameters
  # are an (weighted?) average from shared parameters of parameters 5 and 3).
  # params_shared_average <- params

  # Save all of the likelihood values for each unit
  # pfilterLikes <- as.data.frame(search1_results$likMat)
  # pfilterLikes$which <- rep(1:search1$NREPS, each = search1$NREPS_EVAL)

  # best_unit_parm_id <- pfilterLikes %>%
  #   tidyr::pivot_longer(  # Convert unit likelihood columns into variables
  #     data = .,
  #     cols = -which,
  #     values_to = "loglik",
  #     names_to = "unit",
  #     names_prefix = "V"
  #   ) %>%
  #   dplyr::group_by(which, unit) %>%  # Which indicates parameter set
  #   dplyr::summarize(logLik = logmeanexp(loglik), # Within each parameter set and unit, estimate the log-likelihood for the unit
  #                    se = logmeanexp(loglik, se = TRUE)[2]) %>%
  #   dplyr::mutate(unit = as.integer(unit)) %>%
  #   dplyr::group_by(unit) %>%
  #   dplyr::arrange(-logLik) %>%  # Within each unit, arrange so best parameter sets and likelihood are first
  #   dplyr::slice_head(n = search2$TOP_N) %>%  # Only interested in the top TOP_N
  #   dplyr::arrange(unit, -logLik)

  # For each top-parameter set, we are going to loop through and change the
  # unit-specific parameter in each set (p) to correspond to the parameter that
  # maximized the unit likelihood. Keep in mind that we need to replicate this
  # search2$NREPS times for each set of parameters.
  #
  # Note that this doesn't guarantee maximization
  # of the model likelihood, as there is coupling between units.
  # for (p in 1:search2$TOP_N) {
  #   for (pname in unit_specific_names) {
  #     for (u in 1:10) {
  #
  #       # For unit u, which indicates the parameter set with the pth best log-likelihood.
  #       top_unit_p_id <- best_unit_parm_id[best_unit_parm_id$unit == u, 'which'][p, ] %>% as.numeric()
  #
  #       # Do the same for params_shared_average, so the unit-specific parameters
  #       # for the two datasets are the same.
  #       params_shared_average[((search2$NREPS * (p-1)) + 1):(search2$NREPS * p), paste0(pname, u)] <- search1_results$params[top_unit_p_id, paste0(pname, u)]
  #     }
  #   }
  # }

  # for the top p parameters, take the average of the shared parameters that
  # were the best for unit units.
  # for (p in 1:search2$TOP_N) {
  #   for (pname in shared_param_names) {
  #
  #     # Get the id's for the top pth parameter set for each unit (length=10)
  #     top_unit_ps <- best_unit_parm_id %>%
  #       dplyr::slice(p) %>%
  #       dplyr::pull(which)
  #
  #     # Set shared parameters to average of the top pth parameter set for each unit.
  #     params_shared_average[((search2$NREPS * (p-1)) + 1):(search2$NREPS * p), paste0(pname, 1:10)] <- mean(search1_results$params[top_unit_ps, paste0(pname, 1)])  # All are same across each unit
  #   }
  # }


  # all_params <- rbind(
  #   params,
  #   params_shared_average
  # )

  all_params <- params

  # rm(params, params_shared_average, params_unit_best, pfilterLikes, p,
  #    pname, top_unit_p_id, top_unit_ps, u, best_unit_parm_id)
  rm(params)

  doRNG::registerDoRNG(327498615)

  t_ibpf_local <- system.time(
    foreach(
      i = 1:(search2$TOP_N * search2$NREPS),
      .packages = c("spatPomp"),
      .combine = c
    ) %dopar% {
      r_params <- all_params[i, ]
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

  rm(all_params)

  pfilterLikes <- data.frame(
    "ll" = NA_real_,
    "starting_set" = rep(top_n_global, each = search2$NREPS * search2$NREPS_EVAL),
    "which" = rep(1:(search2$NREPS * search2$TOP_N), each = search2$NREPS_EVAL)
  )

  ibpf_params <- t(coef(local_ibpf))

  for (i in 1:nrow(ibpf_params)) {
    for (sp in shared_param_names) {
      ibpf_params[i, paste0(sp, 1:10)] <- mean(ibpf_params[i, paste0(sp, 1:10)])
    }
  }

  t_local_bpf <- system.time(
    likMat <- foreach(  # Get log-likelihood for each unit and set of parameters, NREPS_EVAL times each
      i = 1:(search2$NREPS_EVAL*search2$NREPS*search2$TOP_N), .combine = rbind, .packages = 'spatPomp'
    ) %dopar% {
      # p3 <- coef(local_ibpf[[(i-1) %/% search2$NREPS_EVAL + 1]])
      p3 <- ibpf_params[(i-1) %/% search2$NREPS_EVAL + 1, ]
      coef(h3_spat) <- p3
      apply(bpfilter(
        h3_spat, Np = search2$NP_EVAL,
        block_size = 1
      )@block.cond.loglik, 1, sum)
    }
  )

  # Condense unit likelihoods into model likelihood
  pfilterLikes$ll <- rowSums(likMat)

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

  # ibpf_logLik$init_method <- (ibpf_logLik$which - 1) %/% (search2$TOP_N * search2$NREPS) + 1
  search2_results <- list()
  search2_results$logLiks <- ibpf_logLik
  search2_results$params <- ibpf_params
  search2_results$ibpf_time <- t_ibpf_local
  search2_results$bpf_time <- t_local_bpf

  if (isTRUE(search2$KEEP_LIK_MAT)) {
    search2_results$likMat <- likMat
  }

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

  # Contains the best set of parameters, measure for the coupled model
  params <- search2_results$params[rep(top_n_local, each = search3$NREPS), ]

  # Creates set of parameters where the shared come from the best set for
  # the entire model, but the unit parameters are chosen to optimize the
  # unit specific likelihoods.
  # params_unit_best <- params

  # Save all of the likelihood values for each unit
  # pfilterLikes <- as.data.frame(search2_results$likMat)
  # pfilterLikes$which <- rep(1:search2$NREPS, each = search2$NREPS_EVAL)

  # best_unit_parm_id <- pfilterLikes %>%
  #   tidyr::pivot_longer(  # Convert unit likelihood columns into variables
  #     data = .,
  #     cols = -which,
  #     values_to = "loglik",
  #     names_to = "unit",
  #     names_prefix = "V"
  #   ) %>%
  #   dplyr::group_by(which, unit) %>%  # Which indicates parameter set
  #   dplyr::summarize(logLik = logmeanexp(loglik), # Within each parameter set and unit, estimate the log-likelihood for the unit
  #                    se = logmeanexp(loglik, se = TRUE)[2]) %>%
  #   dplyr::mutate(unit = as.integer(unit)) %>%
  #   dplyr::group_by(unit) %>%
  #   dplyr::arrange(-logLik) %>%  # Within each unit, arrange so best parameter sets and likelihood are first
  #   dplyr::slice_head(n = search3$TOP_N) %>%  # Only interested in the top TOP_N
  #   dplyr::arrange(unit, -logLik)

  # For each top-parameter set, we are going to loop through and change the
  # unit-specific parameter in each set (p) to correspond to the parameter that
  # maximized the unit likelihood. Keep in mind that we need to replicate this
  # search3$NREPS times for each set of parameters.
  #
  # Note that this doesn't guarantee maximization
  # of the model likelihood, as there is coupling between units.
  # for (p in 1:search3$TOP_N) {
  #   for (pname in unit_specific_names) {
  #     for (u in 1:10) {
  #
  #       # For unit u, which indicates the parameter set with the pth best log-likelihood.
  #       top_unit_p_id <- best_unit_parm_id[best_unit_parm_id$unit == u, 'which'][p, ] %>% as.numeric()
  #
  #       # set search3$NREPS rows to the best unit parameter for unit u
  #       params_unit_best[((search3$NREPS * (p-1)) + 1):(search3$NREPS * p), paste0(pname, u)] <- search2_results$params[top_unit_p_id, paste0(pname, u)]
  #     }
  #   }
  # }

  # all_params <- rbind(
  #   params,
  #   params_unit_best
  # )
  all_params <- params

  doRNG::registerDoRNG(327498615)

  t_ibpf_local <- system.time(
    foreach(
      i = 1:(search3$TOP_N * search3$NREPS),
      .packages = c("spatPomp"),
      .combine = c
    ) %dopar% {
      r_params <- all_params[i, ]
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

  rm(all_params)

  pfilterLikes <- data.frame(
    "ll" = NA_real_,
    "starting_set" = rep(top_n_local, each = search3$NREPS * search3$NREPS_EVAL),
    "which" = rep(1:(search3$NREPS * search3$TOP_N), each = search3$NREPS_EVAL)
  )

  ibpf_params <- t(coef(local_ibpf))

  for (i in 1:nrow(ibpf_params)) {
    for (sp in shared_param_names) {
      ibpf_params[i, paste0(sp, 1:10)] <- mean(ibpf_params[i, paste0(sp, 1:10)])
    }
  }

  t_local_bpf <- system.time(
    likMat <- foreach(  # Get log-likelihood for each unit and set of parameters, NREPS_EVAL times each
      i = 1:(search3$NREPS_EVAL*search3$NREPS*search3$TOP_N), .combine = rbind, .packages = 'spatPomp'
    ) %dopar% {
      # p3 <- coef(local_ibpf[[(i-1) %/% search3$NREPS_EVAL + 1]])
      p3 <- ibpf_params[(i-1) %/% search3$NREPS_EVAL + 1, ]
      coef(h3_spat) <- p3
      apply(bpfilter(
        h3_spat, Np = search3$NP_EVAL,
        block_size = 1
      )@block.cond.loglik, 1, sum)
    }
  )

  # Condense unit likelihoods into model likelihood
  # pfilterLikes$ll <- apply(likMat, 1, sum)
  pfilterLikes$ll <- rowSums(likMat)

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

  # ibpf_logLik$init_method <- (ibpf_logLik$which - 1) %/% (search3$TOP_N * search3$NREPS) + 1
  search3_results <- list()
  search3_results$logLiks <- ibpf_logLik
  search3_results$params <- ibpf_params
  search3_results$ibpf_time <- t_ibpf_local
  search3_results$bpf_time <- t_local_bpf

  if (isTRUE(search3$KEEP_LIK_MAT)) {
    search3_results$likMat <- likMat
  }

  if (isTRUE(search3$KEEP_TRACES)) {
    search3_results$traces <- lapply(local_ibpf, function(x) x@traces[, c('loglik', names(chol_rw_3.sd@call)[-1])])
  }

  results$search3 <- search3_results

  results
}
