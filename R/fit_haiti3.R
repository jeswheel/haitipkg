#' Fit Haiti3
#'
#' A simple function that takes as input starting parameters, as well as
#' hyperparameters to the IBPF algorithm and fits haiti3 using the given
#' starting parameters.
#'
#' @param start_params a numeric vector of parameters for Haiti 3 model.
#' @param NBPF number of BPF iterations to use for the IBPF algorithms.
#' @param NP Number of particles for the IBPF algorithm.
#' @param SPAT_REGRESSION The regression parameter for shared parameters in the
#'    IBPF algorithm.
#' @param NP_EVAL Number of particles to use in the evaluation of parameter estimates.
#' @param NREPS_EVAL Number of replicates of the BPF likelihood evaluation.
#' @param RW_SD Random Walk SD for the IBPF algorithm.
#' @param COOLING Cooling for the IBPF algorithm
#' @param start_date Starting date used for Haiti 3.
#'
#' @return a numeric vector containing the BPF likelihood estimate of the
#'    calibrated parameters, the standard error of this estimate, followed by
#'    all of the parameters of the model (both estimated and fixed parameters.)
#' @export
#'
#' @examples
fit_haiti3 <- function(
    start_params,
    NBPF = 5,
    NP = 50,
    SPAT_REGRESSION = 0.5,
    NP_EVAL = 100,
    NREPS_EVAL = 6,
    RW_SD = NULL,
    COOLING = 0.5,
    start_date = "2010-11-20"
    ) {

  # Create the model that will be fit to cholera incidence data
  h3_spat <- haiti3_spatPomp(start_date = start_date)

  # Create vectors for the unit and shared parameters
  unit_specific_names <- c("betaB", "foi_add", "aHur", "hHur")

  shared_param_names <- c(
    "mu_B", "XthetaA", "thetaI", "lambdaR", "r", "std_W",
    "epsilon", "k"
  )

  est_param_names <- c(
    unit_specific_names, shared_param_names
  )

  # Add unit numbers to each parameter
  est_param_names_expanded <- paste0(rep(est_param_names, each = 10), 1:10)

  coef(h3_spat) <- start_params

  ibpf_out <- ibpf(
    h3_spat,
    Nbpf = NBPF,
    Np = NP,
    sharedParNames = shared_param_names,
    unitParNames = unit_specific_names,
    spat_regression = SPAT_REGRESSION,
    rw.sd = RW_SD,
    cooling.fraction.50 = COOLING,
    block_size = 1
  )

  ibpf_params <- coef(ibpf_out)

  for (sp in shared_param_names) {
    ibpf_params[paste0(sp, 1:10)] <- mean(ibpf_params[paste0(sp, 1:10)])
  }

  coef(h3_spat) <- ibpf_params
  evals <- replicate(NREPS_EVAL, logLik(bpfilter(h3_spat, Np = NP_EVAL, block_size = 1)))
  ll <- pomp::logmeanexp(evals, se = TRUE)
  names(ll) <- c("logLik", "logLik_se")

  c(ll, ibpf_params)
}
