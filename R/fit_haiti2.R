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


  ##### Fitting Lee et al (2020a) model ######

  # Load the data:
  haiti <- haiti2_data()

  ### Fit epidemic period ###
  # Create epidemic model
  org_epi <- haiti2()

  # Save default parameters (many are fixed)
  org_params <- org_epi@params

  # Get objective function to be minimized
  org_ofun <- pomp::traj_objfun(
    org_epi, est = c("Mu", "Beta", "BetaW", "v"),
    params=org_params
  )

  # Create parameter vector of values to be fit
  theta <- org_params[c("Mu", "Beta", "BetaW", "v")]

  # Fit the parameters (log-transformed)
  org_fit <- subplex::subplex(par = log(theta), fn = org_ofun)

  # Store parameters in vector
  org_params[c("Mu", "Beta", "BetaW", "v")] <- exp(org_fit$par)


  ### Fit endemic period ###
  # Create endemic model
  org_end <- haiti2(region = "after")

  # Get default endemic parameters
  org_params2 <- org_end@params

  # Initialize values to be fit at epidemic fit values
  org_params2[c("Mu", "Beta", "BetaW", "v")] <- exp(org_fit$par)

  # Get function to be optimized for endemic model
  org_ofun2 <- pomp::traj_objfun(org_end, est = c("Mu", "v"), params = org_params2)

  # Fit endemic model parameters
  org_fit2 <- subplex::subplex(par = log(org_params2[c("Mu", "v")]), fn = org_ofun2)

  # Save endemic model parameters
  org_params2[c("Mu", "v")] <- exp(org_fit2$par)

  # Save the liklihood of the fitted model
  mod2_leeFit_ll <- -(org_fit$value + org_fit2$value)
  mod2_leeFit_num_parms <- 26  # Six parameters, 20 hidden states.

  rm(
    org_end, org_epi, org_fit, org_fit2,
    org_ofun, org_ofun2, theta
  )

  gc()

  ##### Fitting our joint model ######

  # Create model and save initial values for each parameter
  h2_epi <- haiti2()  # Default is epidemic period
  h2_epi_params <- h2_epi@params

  # Get objective function to minimize
  epi_ofun <- pomp::traj_objfun(
    h2_epi,
    est = c('Mu', 'Beta', 'BetaW', 'v', 'Delta', 'phase'),
    params = h2_epi_params
  )

  theta <- h2_epi_params[c("Mu", "Beta", "BetaW", "v", "Delta")]  # parameters to fit
  # Note that phase is not log-transformed, se it's treated slightly differently
  # Fit the epi model
  h2_epi_fit <- subplex::subplex(
    par = c(log(theta), h2_epi_params['phase']),
    fn = epi_ofun
  )

  # Save fitted parameters
  h2_epi_params[c("Mu","Beta","BetaW","v","Delta")] <- exp(h2_epi_fit$par[c("Mu","Beta","BetaW","v","Delta")])
  h2_epi_params['phase'] <- h2_epi_fit$par['phase']

  # Now that the model for the epidemic phase is fit, we need to get the values of
  # the hidden states and use that as initial values for fitting the model on the
  # endemic phase:
  h2_epi_sim <- pomp::trajectory(
    h2_epi,
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

  # Get objective function:
  end_ofun <- pomp::traj_objfun(
    h2_end,
    est = c("Mu","Beta","BetaW","v","Delta"),
    params = h2_end_params
  )

  # Fit model:
  h2_end_fit <- subplex::subplex(
    par = log(h2_end_params[c("Mu","Beta","BetaW","v","Delta")]),
    fn = end_ofun
  )

  h2_end_params[c("Mu","Beta","BetaW","v","Delta")] <- exp(h2_end_fit$par)

  rm(
    h2_end, h2_end_fit, h2_epi,
    h2_epi_fit, h2_epi_sim, haiti,
    end_ofun, epi_ofun, theta, IVPS
  )

  gc()

  model2_fit_results <- list()

  model2_fit_results[['epi_params']] <- h2_epi_params
  model2_fit_results[['end_params']] <- h2_end_params
  model2_fit_results[['leeFit_ll']] <- mod2_leeFit_ll
  model2_fit_results[['leeFit_n_params']] <- mod2_leeFit_num_parms

  return(model2_fit_results)
}
