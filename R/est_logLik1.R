#' Estimate Log-Likelihood
#'
#' This function estimates the log-likelihood for model 1.
#'
#' As of Feb 20, 2022, Model 1 can be written in 4 different forms:
#' \itemize{
#'    \item Basic: aggregated, single time period `pomp` object (`haiti1()`).
#'    \item Joint: aggregated, joint time periods `pomp` object (`haiti1_joint()`).
#'    \item Disagg: disaggregated `pomp` object (`haiti1_disagg()`).
#'    \item Panel: disaggregated `panelPomp` object (`haiti1_disagg()`).
#' }
#'
#' Currently, this function supports estimation for the following forms: Basic, Joint, Panel
#'
#' @param version In c("basic", "joint", "panel").
#' @param model A `pomp` or `panelPomp` object matching the specified model version.
#' @param Np Number of Particles used in the (potentially block) Particle Filter
#' @param nreps Number of repeated particle filters
#' @param ncores Number of cores used to conduct the repeated particle filters.
#'    The parallelization is done over `nreps`, so ideally `ncores` should divide
#'    `nreps`.
#'
#' @examples
#' m1 <- haiti1_joint()
#' loglik_est <- est_logLik1(version = "joint", model = m1, Np = 1000, nreps = 10, ncores = 5)
#'
#' @importFrom magrittr %>%
#' @importFrom foreach %do%
#' @importFrom foreach foreach
#'
#' @export
est_logLik1 <- function(version = NULL, model = NULL, Np = 100, nreps = 10, ncores = 1) {
  doParallel::registerDoParallel(ncores)
  if (version == 'panel') {
    # Get logLikelihood Matrix
    logliks <- foreach(i=1:nreps, .combine = rbind, .packages = "panelPomp") %dopar% {
      panelPomp::unitlogLik(panelPomp::pfilter(model, Np = Np))
    }
    logliks <- panel_logmeanexp(logliks, MARGIN = 2, se = TRUE)
    return(logliks)

  } else if (version == 'joint' || version == 'basic' || version == 'disagg') {
    logliks <- foreach (i=1:nreps, .combine=c) %dopar% {
      model %>% pfilter(Np = Np)
    }
    logliks <- sapply(logliks, logLik)
    logliks <- logmeanexp(logliks, se = TRUE)
    return(logliks)
  } else {
    stop(paste0("Version = \"", version, "\" is not a valid input.
                Version must be in {\"basic\", \"panel\", \"joint\", \"disagg\"}."))
  }
}
