#' Estimate Log-Likelihood
#'
#' This function estimates the log-likelihood for model 2.
#'
#' @param model Input is a spatPomp object of the model whose Log-likelihood
#' is to be estimated. The spatPomp object should already have parameter values
#' stored in it.
#'
#' @return The negative log likelihood of the model.
#'
#' @export
est_logLik2 <- function(model) {

  # Check that input is correct type
  if (class(model) != "spatPomp"){
    stop("Input is not a spatPomp model.")
  }

  # Get objective function for model
  ofun <- traj_objfun(model)

  # Get log-Likelihood
  ll <- ofun(h2@params)

  return(ll)
}
