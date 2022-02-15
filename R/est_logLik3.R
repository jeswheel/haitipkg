#' Estimate Log-Likelihood
#'
#' This function estimates the log-likelihood for model 3.
#'
#' As of Feb 11, 2022, Model 3 can be written in 3 different forms:
#' \itemize{
#'    \item Coupled `pomp` object (`haiti3_correct()`).
#'    \item Independent `panelPomp` object (`haiti3_panel()`).
#'    \item Coupled `spatPomp` object (`haiti3_spatPomp()`).
#' }
#'
#'
#' @param version in c("panelPomp", "spatPomp", "pomp").
#' @param params Parameters that will be used in the log-likelihood estimation.
#'    Note that these parameters need to match the version of the model that
#'    will be used, otherwise a warning will be given and the default parameters
#'    for the version requested will be used.
#' @param Np Number of Particles used in the (potentially block) Particle Filter
#' @param nreps Number of repeated particle filters
#' @param ncores Number of cores used to conduct the repeated particle filters.
#'    The parallelization is done over `nreps`, so ideally `ncores` should divide
#'    `nreps`.
#'
#' @importFrom magrittr %>%
#' @importFrom foreach %do%
#' @importFrom foreach foreach
#'
#' @export
est_logLik3 <- function(version = 'panelPomp', params = NULL,
                        Np = 100, nreps = 10, ncores = 1) {

  doParallel::registerDoParallel(ncores)
  doRNG::registerDoRNG(38763911)


  if (version == 'panelPomp') {

    # Load the model
    mod <- haiti3_panel(start_time = "2010-10-23", B0 = TRUE)
    temp_params <- panelPomp::coef(mod)

    # Check to make sure that we have parameters that will work for the panelPomp
    if (is.null(params)) {  # User did not provide input parameters

      params <- temp_params

    } else if (!all(names(params) %in% names(temp_params)) && !all(names(temp_params) %in% names(params))) {  # not the correct parameters

      warning("Names of provided parameters do not match with the parameters
              of the panelPomp version of Model 3; default parameter values will
              be used instead.")

      parms <- temp_params  # Set parameters to defaults.
    }

    # Get logLikelihood Matrix
    pf3_loglik_matrix <- foreach(i=1:nreps,
                                 .combine = rbind,
                                 .packages = "panelPomp") %dopar% {
                                   panelPomp::unitlogLik(
                                     panelPomp::pfilter(
                                       mod,
                                       params = params,
                                       Np = Np
                                     )
                                   )
                                 }

    # Retrun final estimate of the log-likelihood.
    return(panel_logmeanexp(pf3_loglik_matrix, MARGIN = 2, se = TRUE))

  } else if (version == 'spatPomp') {

    mod <- haiti3_spatPomp()
    temp_params <- panelPomp::coef(mod)

    # Check to make sure that we have parameters that will work for the spatPomp
    if (is.null(params)) {  # User did not provide input parameters

      params <- temp_params

    } else if (!all(names(params) %in% names(temp_params)) && !all(names(temp_params) %in% names(params))) {  # not the correct parameters

      warning("Names of provided parameters do not match with the parameters
              of the panelPomp version of Model 3; default parameter values will
              be used instead.")

      parms <- temp_params  # Set parameters to defaults.
    }

    pomp::coef(mod) <- params

    pf <- foreach (i = 1:nreps, .combine = c, .packages = "spatPomp") %dopar% {
      mod %>% spatPomp::bpfilter(
        Np = Np,
        params = params,
        block_list = list(c(1, 2), c(5, 6), c(7), c(3, 4, 9), c(8, 10))
      )
    }

    ll <- logLik(pf)
    return(logmeanexp(ll[!is.na(ll)], se = TRUE))

  } else if (version == "pomp") {  # Coupled pomp Model.
    warning("Due to the spatial structure of the model, a better estimate of the
            log-likelihood will be obtained by running `version = 'spatPomp' for
            the coupled model.")

    # Create coupled version of the model
    mod <- haiti3_correct(dt_yrs = 1/365.25)
    temp_params <- pomp::coef(mod)

    # Check to make sure that we have parameters that will work for the spatPomp
    if (is.null(params)) {  # User did not provide input parameters

      params <- temp_params

    } else if (!all(names(params) %in% names(temp_params)) && !all(names(temp_params) %in% names(params))) {  # not the correct parameters

      warning("Names of provided parameters do not match with the parameters
              of the panelPomp version of Model 3; default parameter values will
              be used instead.")

      parms <- temp_params  # Set parameters to defaults.
    }

    # Set coupled model parameters to fitted panelPomp parameters.
    coef(mod) <- params

    pf <- foreach (i=1:nreps, .combine=c) %dopar% {
      mod_coupled %>% pfilter(Np = Np, params = params)
    }
    ll <- logLik(pf)

    return(logmeanexp(ll[!is.na(ll)], se = TRUE))
  } else {
    stop(paste0("Version = \"", version, "\" is not a valid input. Version must be in {\"spatPomp\", \"panelPomp\", \"pomp\"}."))
  }
}
