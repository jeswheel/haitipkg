#' Fit Haiti 1
#'
#' This simple function takes as input starting values for Model 1 parameters,
#' as well as hyperparameters to the IF2 algorithm and fits Model 1 using these
#' input values.
#'
#' @param start_params A numeric vector of parameters for Model 1.
#' @param NMIF Integer number of IF2 iterations.
#' @param NP Integer number of particles for the IF2 algorithm.
#' @param NP_EVAL Integer number of particles to use in the particle filter
#'    that is used to estimate the log-likelihood of the fitted model.
#' @param NREPS_EVAL Number of repeated likelihood evaluations.
#' @param RW_SD Random walk standard deviations for the IF2 algorithm.
#' @param COOLING Cooling rate for the IF2 algorithm.
#'
#' @return A numeric vector containing the log-likelihood estimate of the model
#'    parameters, an estimate of the standard error of the estimated likelihood,
#'    and the parameter vector resulting from the IF2 algorithm.
#' @export
fit_haiti1 <- function(
    start_params, NMIF, NP,
    NP_EVAL, NREPS_EVAL, RW_SD,
    COOLING = 0.5
    ) {

  # Load the model, allow breakpoint in rho, tau, and sig_sq, but not beta or nu
  mod_1 <- haiti1_joint(
    rho_flag = FALSE,
    tau_flag = TRUE,
    sig_sq_flag = TRUE,
    beta_flag = FALSE,
    nu_flag = FALSE
  )

  coef(mod_1) <- start_params

  if2_out <- mif2(
    mod_1,
    Np = NP,
    Nmif = NMIF,
    cooling.fraction.50 = COOLING,
    rw.sd = RW_SD
  )

  if2_params <- coef(if2_out)
  coef(mod_1) <- if2_params

  evals <- replicate(
    NREPS_EVAL,
    logLik(pomp::pfilter(mod_1, Np = NP_EVAL, params = if2_params))
  )

  ll <- pomp::logmeanexp(evals, se = TRUE)
  names(ll) <- c("logLik", "logLik_se")

  c(ll, if2_params)
}
