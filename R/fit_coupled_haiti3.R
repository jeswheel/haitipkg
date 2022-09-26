#' Fit Model 3 SpatPomp
#'
#' This function fits the spatPomp version of Model 3 to cholera incidence
#' data from Oct 2010 - Jan 2019 in Haiti. The data are measured at a
#' weekly timescale.
#'
#' @import foreach
#' @import doRNG
#' @import spatPomp
#'
#' @importFrom magrittr %>%
#'
#' @export
fit_coupled_haiti3 <- function(ncores = 3) {

  #
  ##
  ### Global search
  ##
  #

  h3_spat <- haiti3_spatPomp()

  # Create vectors for the unit and shared parameters
  unit_specific_names <- c("betaB", "foi_add")
  shared_param_names <- c(
    "mu_B", "XthetaA", "thetaI", "lambdaR", "r", "std_W",
    "epsilon", "k", "sigma"
  )
  est_param_names <- c(
    unit_specific_names, shared_param_names
  )

  # Add unit numbers to each parameter
  est_param_names_expanded <- paste0(rep(est_param_names, each = 10), 1:10)

  # Create rw.sd for each parameter
  reg_rw.sd <- lapply(est_param_names_expanded, function(x) 0.02)
  names(reg_rw.sd) <- est_param_names_expanded
  chol_rw.sd <- do.call(rw.sd, reg_rw.sd)

  # Get lower bound for unit parameters (global search)
  min_val <- 1e-8
  unit_lb <- rep(c(min_val, min_val), each = 10)
  names(unit_lb) <- paste0(rep(unit_specific_names, each = 10), 1:10)

  # Get upper bound for unit parameters (global search)
  unit_ub <- rep(c(50, 1e-5), each = 10)
  names(unit_ub) <- paste0(rep(unit_specific_names, each = 10), 1:10)

  # Get lower bound for shared parameters (global search)
  shared_lb <- rep(c(5, min_val, min_val, min_val,
                     min_val, min_val, 0.25, 5, 0.01))
  names(shared_lb) <- shared_param_names

  # Get upper bound for shared parameters (global search)
  shared_ub <- c(300, 1, 2e-3, 5, 1.2, 0.15, 1, 1000, 0.5)
  names(shared_ub) <- shared_param_names

  # Create data.frame with random unit parameters
  guesses_unit <- pomp::runif_design(
    lower = unit_lb,
    upper = unit_ub,
    nseq = 20
  )

  # Create data.frame with random shared parameters
  guesses_shared <- pomp::runif_design(
    lower = shared_lb,
    upper = shared_ub,
    nseq = 20
  )

  # Need to duplicate each of the shared parameter columns
  guesses_shared <- guesses_shared[, rep(1:length(shared_param_names), each = 10)]
  colnames(guesses_shared) <- paste0(rep(shared_param_names, each = 10), 1:10)

  # Combine the unit and shared parameters
  guesses <- cbind(guesses_shared, guesses_unit)

  # We need to add fixed parameters
  all_params <- coef(h3_spat)
  fixed_params <- all_params[!names(all_params) %in% colnames(guesses)]
  fixed_mat <- matrix(rep(fixed_params, 20), byrow = TRUE, nrow = 20)
  colnames(fixed_mat) <- names(all_params[!names(all_params) %in% colnames(guesses)])
  guesses_all <- cbind(guesses, fixed_mat)[names(coef(h3_spat))]

  doParallel::registerDoParallel(ncores)
  registerDoRNG(2198635)

  foreach(
    i = 1:20,  # TODO: change to a varaible argument.
    .packages = c("spatPomp"),
    .combine = c
  ) %dopar% {
    r_params <- unlist(guesses_all[i, ])
    coef(h3_spat) <- r_params

    ibpf_out <- ibpf(
      h3_spat,
      Nbpf = 20,  # TODO: change to a variable
      Np = 2000,  # TODO: change to a variable
      sharedParNames = shared_param_names,
      unitParNames = unit_specific_names,
      spat_regression = 0.5,
      rw.sd = chol_rw.sd,
      cooling.fraction.50 = 0.5,
      block_size = 1,
      params = r_params
    )
  } -> Global_ibpf

}
