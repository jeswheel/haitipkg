#' Fit model 2
#'
#' This function fits the parameters to cholera incidence data in Haiti. This
#' function has two main parts:
#'   \enumerate{
#'     \item Fit Lee et al (2018) version of the model, to obtain model likelihoods.
#'     \item Fit our joint version of the model.
#'   }
#'
#' @importFrom magrittr %>%
#' @return List containing four objects: \describe{
#'   \item{epi_params}{A numerical vector containing the parameters to the epidemic period.}
#'   \item{end_params}{A numerical vector containing the parameters to the endemic period.}
#'   \item{leeFit_ll}{numeric containing the likelihood of the Lee et al (2020a) model.}
#'   \item{leeFit_n_params}{numeric containing the number of parameters fit by Lee et al (2020a).}
#' }
#'
#' @examples
#' h2_fit_results <- fit_haiti2()
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

  lee_epi_ll <- -epi_fit_v$value
  lee_end_ll <- -end_fit_v$value


  mod2_leeFit_ll <- lee_epi_ll - sum(apply(log(h2_epi@data + 1), 1, sum)) + (lee_end_ll - sum(apply(log(h2_end@data + 1), 1, sum, na.rm = TRUE)))
  mod2_leeFit_num_parms <- 26

  rm(
    epi_params, epi_ofun_v, theta, epi_fit_v, epi_params_update_v,
    end_params, end_ofun_v, end_fit_v, end_params_update_v, lee_epi_ll,
    lee_end_ll
  )

  gc()

  ##### Fitting our joint model ######

  h2_epi_temp <- haiti2()

  # Create model and save initial values for each parameter
  h2_epi <- haiti2(cutoff = 10000)  # Default is epidemic period
  h2_epi_params <- h2_epi@params

  # Get objective function to minimize
  epi_ofun <- pomp::traj_objfun(
    h2_epi,
    est = c('Mu', 'Beta', 'BetaW', 'v', 'sigma', 'phase'),
    params = h2_epi_params
  )

  theta <- h2_epi_params[c("Mu", "Beta", "BetaW", "v", "sigma")]  # parameters to fit
  # Note that phase is not log-transformed, se it's treated slightly differently
  # Fit the epi model
  h2_epi_fit <- subplex::subplex(
    par = c(log(theta), h2_epi_params['phase']),
    fn = epi_ofun
  )

  # Save fitted parameters
  h2_epi_params[c("Mu","Beta","BetaW","v", "sigma")] <- exp(h2_epi_fit$par[c("Mu", "Beta", "BetaW", "v", "sigma")])
  h2_epi_params['phase'] <- h2_epi_fit$par['phase']

  # Now that the model for the epidemic phase is fit, we need to get the values of
  # the hidden states and use that as initial values for fitting the model on the
  # endemic phase:
  h2_epi_sim <- pomp::trajectory(
    h2_epi,
    times = h2_epi_temp@times,
    t0 = h2_epi_temp@t0,
    params = h2_epi_params,
    format = 'data.frame'
  )
  IVPS <- h2_epi_sim[nrow(h2_epi_sim),
                     !colnames(h2_epi_sim) %in% c('year', '.id')] %>%
    unlist()

  # Now we fit the endemic model:
  h2_end <- haiti2(region = 'after', joint = TRUE)
  h2_end_params <- h2_end@params

  # Initialize parameters based on previous result:
  h2_end_params[1:24] <- h2_epi_params[1:24]

  # Initialize hidden states:
  h2_end_params[25:length(h2_end_params)] <- IVPS

  # TODO: Uncomment to fit endemic period seperately
  # # Get objective function:
  # end_ofun <- pomp::traj_objfun(
  #   h2_end,
  #   est = c("Mu","Beta","BetaW","v","Delta", 'phase'),
  #   params = h2_end_params
  # )
  #
  # # Fit model:
  # h2_end_fit <- subplex::subplex(
  #   par = c(log(h2_end_params[c("Mu","Beta","BetaW","v","Delta")]), h2_end_params['phase']),
  #   fn = end_ofun
  # )
  #
  # h2_end_params[c("Mu","Beta","BetaW","v","Delta")] <- exp(h2_end_fit$par[c("Mu","Beta","BetaW","v","Delta")])
  # h2_end_params['phase'] <- h2_end_fit$par['phase']

  rm(
    h2_end,
    # h2_end_fit,
    h2_epi,
    h2_epi_fit, h2_epi_sim,
    # end_ofun,
    epi_ofun, theta, IVPS
  )

  gc()

  model2_fit_results <- list()

  model2_fit_results[['epi_params']] <- h2_epi_params
  model2_fit_results[['end_params']] <- h2_end_params
  model2_fit_results[['leeFit_ll']] <- mod2_leeFit_ll
  model2_fit_results[['leeFit_n_params']] <- mod2_leeFit_num_parms

  return(model2_fit_results)
}
