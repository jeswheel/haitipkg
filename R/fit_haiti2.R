#' Fit model 2
#'
#' This function fits the parameters in model 2.
#'
#' @importFrom subplex subplex
#' @param model Input is a spatPomp object of the model that is to be fit.
#' @param params Vector of characters representing the parameters to be fit.
#'
#' @return Spatpomp object with the newly fitted parameters in it.
#'
#' @export
fit_haiti2 <- function(model, params=c("Mu","Beta","BetaW","v","VR")) {

  # Check that input is correct type
  if (class(model) != "spatPomp"){
    stop("Input is not a spatPomp model.")
  }

  # Get objective function for model
  ofun <- traj_objfun(model, est=params)

  # Get vector of initial param values. Using the param values in the model.
  theta <- model@params[params]

  # Fit the model
  fit <- subplex(par=log(theta), fn=ofun)

  # Update params in model with fitted params
  model@params[params] <- exp(fit$par)

  return(model)
}
