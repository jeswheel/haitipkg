# TODO: Before actually running anything, create a "lab-log" of all of the things that we tried, that way
# if we used "banked" results, it will be possible to completely recreate these banked results. To start,
# keep track of the haitipkg version number and commit hash of the original submission (we will probably use those banked results), as well as the commit hash for the submission. Then figure out a way to track everything that we
# submit to great lakes (at least everything that we might use, I'm guessing a few low RUN_LEVELS will be
# conducted to smooth out issues before running a big job).

# TODO: Document.
fit_haiti3 <- function(
  # TODO: change search list defaults.
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
    # TODO: change search list defaults.
    search2 = list(
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
    prof_parameter,
    prof_values,
    nprof,
    start_date = "2010-11-20"
    ) {

  COOLING <- 0.5
  KEEP_LIK_MAT <- TRUE

  # Create the model that will be fit to cholera incidence data
  h3_spat <- haiti3_spatPomp(start_date = start_date)

  # Create a list to save all of the results.
  results <- list()

  # Create vectors for the unit and shared parameters
  unit_specific_names <- c("betaB", "foi_add", "aHur", "hHur")
  # unit_names_expanded <- paste0(rep(unit_specific_names, each = 10), 1:10)

  shared_param_names <- c(
    "mu_B", "XthetaA", "thetaI", "lambdaR", "r", "std_W",
    "epsilon", "k"
  )

  if (prof_parameter %in% shared_param_names) {
    prof_cols <- matrix(rep(rep(prof_values, 10), each = nprof), ncol = 10)
    colnames(prof_cols) <- paste0(prof_parameter, 1:10)
  } else {
    prof_cols <- matrix(rep(prof_values, each = nprof), ncol = 1)
    colnames(prof_cols) <- prof_parameter
  }

  S1_reps <- nprof * length(prof_values)

  est_param_names <- c(
    unit_specific_names, shared_param_names
  )

  # Add unit numbers to each parameter
  # est_param_names_expanded <- paste0(rep(est_param_names, each = 10), 1:10)

  # TODO: This parameter will be removed, just stick with what we found worked.
  # if (is.null(search1$RW_SD)) {
  #   stop("RW_SD must be supplied for search 1.")
  # } else {
  #   chol_rw_1.sd <- search1$RW_SD
  # }

  # TODO: This parameter will be removed, just stick with what we found worked.
  # if (nsearches >= 2 && is.null(search2$RW_SD)) {
  #   stop("RW_SD must be supplied for search 2.")
  # } else {
  #   chol_rw_2.sd <- search2$RW_SD
  # }

  # Unit-lower bounds
  unit_lower_bounds <- c(
    'betaB1' = 2.6,
    'betaB2' = 9,
    'betaB3' = 13,
    'betaB4' = 13,
    'betaB5' = 2.5,
    'betaB6' = 15,
    'betaB7' = 5,
    'betaB8' = 0.4,
    'betaB9' = 5,
    'betaB10' = 6,
    'foi_add1' = 5e-8,
    'foi_add2' = 1e-8,
    'foi_add3' = 1e-7,
    'foi_add4' = 5e-8,
    'foi_add5' = 2e-7,
    'foi_add6' = 2.5e-7,
    'foi_add7' = 1.1e-7,
    'foi_add8' = 1e-8,
    'foi_add9' = 1.5e-7,
    'foi_add10' = 1e-7,
    'aHur3' = 0.01,
    'aHur9' = 10,
    'hHur3' = 30,
    'hHur9' = 40,
    'Iinit3' = 0.9 / 468301,
    'Iinit4' = 4.8 / 342525
  )

  # Unit upper-bounds
  unit_upper_bounds <- c(
    'betaB1' = 8,
    'betaB2' = 35,
    'betaB3' = 50,
    'betaB4' = 50,
    'betaB5' = 10,
    'betaB6' = 50,
    'betaB7' = 20,
    'betaB8' = 2,
    'betaB9' = 20,
    'betaB10' = 20,
    'foi_add1' = 8e-7,
    'foi_add2' = 6e-7,
    'foi_add3' = 5e-7,
    'foi_add4' = 2.5e-7,
    'foi_add5' = 7.5e-7,
    'foi_add6' = 6e-7,
    'foi_add7' = 4.5e-7,
    'foi_add8' = 4e-7,
    'foi_add9' = 3.5e-7,
    'foi_add10' = 2.8e-7,
    'aHur3' = 50,
    'aHur9' = 45,
    'hHur3' = 120,
    'hHur9' = 100,
    'Iinit3' = 35 / 468301,
    'Iinit4' = 35 / 342525
  )

  shared_lower_bounds <- c(
    "mu_B" = 350,
    "XthetaA" = 0.02,
    "thetaI" = 2.5e-05,
    "lambdaR" = 0.5,
    "r" = 0.5,
    "std_W" = 0.02,
    "epsilon" = 0.65,
    "k" = 75
  )

  shared_upper_bounds <- c(
    "mu_B" = 750,
    "XthetaA" = 0.25,
    "thetaI" = 1e-04,
    "lambdaR" = 3.5,
    "r" = 1.75,
    "std_W" = 0.03,
    "epsilon" = 0.99,
    "k" = 175
  )

  est_u_names <- names(unit_lower_bounds)
  partial_u_names <- est_u_names[est_u_names != prof_parameter]
  partial_s_names <- shared_param_names[shared_param_names != prof_parameter]

  unit_lower_bounds <- unit_lower_bounds[partial_u_names]
  unit_upper_bounds <- unit_upper_bounds[partial_u_names]

  shared_lower_bounds <- shared_lower_bounds[partial_s_names]
  shared_upper_bounds <- shared_upper_bounds[partial_s_names]

  # Create data.frame with random unit parameters
  guesses_unit <- pomp::runif_design(
    lower = unit_lower_bounds,
    upper = unit_upper_bounds,
    nseq = S1_reps
  )

  guesses_shared <- pomp::runif_design(
    lower = shared_lower_bounds,
    upper = shared_upper_bounds,
    nseq = S1_reps
  )

  # Need to duplicate each of the shared parameter columns
  guesses_shared <- guesses_shared[, rep(partial_s_names, each = 10)]
  colnames(guesses_shared) <- paste0(rep(partial_s_names, each = 10), 1:10)

  # Combine the unit and shared parameters
  guesses <- cbind(guesses_shared, guesses_unit)
  guesses <- dplyr::bind_cols(prof_cols, guesses)

  # We need to add fixed parameters
  all_params <- coef(h3_spat)
  fixed_params <- all_params[!names(all_params) %in% colnames(guesses)]
  fixed_mat <- matrix(
    rep(fixed_params, S1_reps),
    byrow = TRUE, nrow = S1_reps
  )
  colnames(fixed_mat) <- names(all_params[!names(all_params) %in% colnames(guesses)])

  # Combine estimated and fixed parameters, and reorder based on original order.
  guesses_all <- cbind(guesses, fixed_mat)[, names(coef(h3_spat))]

  # Memory clean-up
  rm(guesses_unit, guesses_shared, fixed_mat, min_val, shared_lower_bound,
     shared_upper_bound, unit_lower_bound, unit_upper_bound,  all_params, fixed_params,
     est_param_names)
  gc()

  t_global <- system.time(
    foreach(
      i = 1:S1_reps,
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
      )
    } -> Global_ibpf
  )

  rm(guesses_all)

  pfilterLikes <- data.frame(
    "ll" = rep(0, search1$NREPS_EVAL*S1_reps),
    "which" = rep(1:S1_reps, each = search1$NREPS_EVAL)
  )

  ibpf_params <- t(coef(Global_ibpf))

  for (i in 1:nrow(ibpf_params)) {
    for (sp in partial_s_names) {
      ibpf_params[i, paste0(sp, 1:10)] <- mean(ibpf_params[i, paste0(sp, 1:10)])
    }
  }

  t_global_bpf <- system.time(
    likMat <- foreach(  # Get log-likelihood for each unit and set of parameters, NREPS_EVAL times each
      i = 1:(search1$NREPS_EVAL * S1_reps), .combine = rbind, .packages = 'spatPomp'
    ) %dopar% {
      p3 <- ibpf_params[(i - 1) %/% search1$NREPS_EVAL + 1, ]

      coef(h3_spat) <- p3
      apply(bpfilter(
        h3_spat, Np = search1$NP_EVAL,
        block_size = 1
      )@block.cond.loglik, 1, sum)
    }
  )

  if (prof_parameter %in% shared_param_names) {
    prof_parameter <- paste0(prof_parameter, "1")
  }

  # Condense unit likelihoods into model likelihood
  # pfilterLikes$ll <- apply(likMat, 1, sum)
  pfilterLikes$ll <- rowSums(likMat)  # Faster

  # Group by model parameter set.
  ibpf_logLik <- pfilterLikes |>
    dplyr::group_by(which) |>
    dplyr::summarize(logLik = logmeanexp(ll),
              se = logmeanexp(ll, se = TRUE)[2])

  ibpf_logLik <- bind_cols(ibpf_logLik, prof_cols[, 1, drop = FALSE])

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

  top_n_global <- search1_results$logLiks %>%
    dplyr::group_by_at(prof_parameter) %>%
    dplyr::slice_max(order_by = dplyr::desc(logLik), n = search2$TOP_N) %>%
    dplyr::pull(which)

  # top_n_global <- order(-search1_results$logLiks$logLik)[1:search2$TOP_N]

  # Contains the best set of parameters, measure for the coupled model
  all_params <- search1_results$params[rep(top_n_global, each = search2$NREPS), ]

  t_ibpf_local <- system.time(
    foreach(
      i = 1:nrow(all_params),
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
    "which" = rep(1:(search2$NREPS * search2$TOP_N * length(prof_values)), each = search2$NREPS_EVAL)
  )

  ibpf_params <- t(coef(local_ibpf))

  for (i in 1:nrow(ibpf_params)) {
    for (sp in partial_s_names) {
      ibpf_params[i, paste0(sp, 1:10)] <- mean(ibpf_params[i, paste0(sp, 1:10)])
    }
  }

  # TODO: Incorrect.
  t_local_bpf <- system.time(
    likMat <- foreach(  # Get log-likelihood for each unit and set of parameters, NREPS_EVAL times each
      i = 1:(search2$NREPS_EVAL*search2$NREPS*search2$TOP_N*length(prof_values)), .combine = rbind, .packages = 'spatPomp'
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
  ibpf_logLik_temp <- pfilterLikes |>
    dplyr::group_by(which) |>
    dplyr::summarize(logLik = logmeanexp(ll),
                     se = logmeanexp(ll, se = TRUE)[2])

  ibpf_logLik <- pfilterLikes |>
    dplyr::group_by(which) |>
    dplyr::summarise(starting_set = dplyr::first(starting_set)) |>
    dplyr::right_join(
      y = ibpf_logLik_temp,
      by = "which"
    ) |>
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

  results
}
