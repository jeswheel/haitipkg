#' Fit model 2
#'
#' This function fits the parameters to cholera incidence data in Haiti. This
#' function has two main parts:
#'   \enumerate{
#'     \item Fit Lee et al (2018) version of the model, to obtain model likelihoods.
#'     \item Fit our joint version of the model.
#'   }
#'
#' @return List containing four objects: \describe{
#'   \item{h2_params}{A numerical vector containing the parameters to the epidemic period.}
#'   \item{leeFit_ll}{numeric containing the likelihood of the Lee et al (2020a) model.}
#'   \item{leeFit_n_params}{numeric containing the number of parameters fit by Lee et al (2020a).}
#' }
#'
#' @examples
#' \dontrun{h2_fit_results <- fit_haiti2()}
#'
#' @export
fit_haiti2 <- function() {

  #### Calculate the LL of Lee et al (2020) model

  # Load the epi model, default parameters based on lee et al
  h2_epi <- haiti2()
  epi_params <- h2_epi@params

  # Create objective function to maximize
  epi_ofun_v <- traj_objfun(
    h2_epi,
    est = c("v"),
    paramnames = c("v")
  )

  # Save vector "theta" that will be used as starting values
  theta <- epi_params[c('v')]

  # Fit parameter "v"
  epi_fit_v <- subplex::subplex(par = log(theta), fn = epi_ofun_v)

  # Update parameters
  epi_params_update_v <- epi_params
  epi_params_update_v[c('v')] <- exp(epi_fit_v$par)

  # Now we do the same for the endemic period:
  h2_end <- haiti2(region = 'after')

  end_params <- h2_end@params
  end_ofun_v <- traj_objfun(
    h2_end,
    est = c("v"),
    paramnames = c("v")
  )
  theta <- end_params[c('v')]
  end_fit_v <- subplex::subplex(par = log(theta), fn = end_ofun_v)
  end_params_update_v <- end_params
  end_params_update_v[c('v')] <- exp(end_fit_v$par)

  end_ofun_v(par = end_fit_v$par) + epi_ofun_v(par = epi_fit_v$par)

  lee_epi_ll <- -epi_fit_v$value
  lee_end_ll <- -end_fit_v$value

  mod2_leeFit_ll <- lee_epi_ll + lee_end_ll
  mod2_leeFit_num_parms <- 6

  rm(
    epi_params, epi_ofun_v, theta, epi_fit_v, epi_params_update_v,
    end_params, end_ofun_v, end_fit_v, end_params_update_v, lee_epi_ll,
    lee_end_ll
  )

  gc()

  ##### Fitting our joint model ######

  # Create model and save initial values for each parameter
  h2 <- haiti2(cutoff = 10000, measure = "log")  # Default is epidemic period
  h2_params <- h2@params

  est_params <- c('Mu', 'Beta', 'BetaW', 'v', 'sigma', 'phase')
  n_fit_params <- length(est_params)

  # Get objective function to minimize
  epi_ofun <- pomp::traj_objfun(
    h2,
    est = est_params,
    params = h2_params
  )

  theta <- h2_params[c("Mu", "Beta", "BetaW", "v", "sigma")]  # parameters to fit
  # Note that phase is not log-transformed, se it's treated slightly differently
  # Fit the epi model
  h2_fit <- subplex::subplex(
    par = c(log(theta), h2_params['phase']),
    fn = epi_ofun
  )

  # Save fitted parameters
  h2_params[c("Mu","Beta","BetaW","v", "sigma")] <- exp(h2_fit$par[c("Mu", "Beta", "BetaW", "v", "sigma")])
  h2_params['phase'] <- h2_fit$par['phase']

  rm(
    h2,
    h2_fit,
    epi_ofun,
    theta
  )

  gc()

  model2_fit_results <- list()

  model2_fit_results[['h2_params']] <- h2_params
  model2_fit_results[["n_fit_params"]] <- n_fit_params
  # model2_fit_results[['end_params']] <- h2_end_params
  model2_fit_results[['leeFit_ll']] <- mod2_leeFit_ll
  model2_fit_results[['leeFit_n_params']] <- mod2_leeFit_num_parms

  return(model2_fit_results)
}
