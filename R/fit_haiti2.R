#' Fit Haiti 2
#'
#' This simple function takes as input an objective function to minimize and
#' parameter initializations, and returns the corresponding maximum likelihood
#' estimate. Note that for this function to work properly, the objective
#' function needs to be created using [pomp::traj_objfun].
#'
#' @param initialization A vector of parameter values that will be used as
#'    the initialization in order to perform optimization. This initialization
#'    must match the necessary parameters in the obj_fun argument.
#' @param obj_fun The objective function to maximize, obtained as the output
#'    of [pomp::traj_objfun] of an instance of Model 2.
#' @param ... Additional arguments to be passed into the [subplex::subplex]
#'    function.
#'
#' @return A vector containing the estimated parameters and the corresponding
#'    likelihood value.
#' @export
fit_haiti2 <- function(initialization, obj_fun, ...) {
  h2_fit <- subplex::subplex(
    par = initialization,
    fn = obj_fun,
    ...
  )
}
