#' Build panelPomp object for Model 1
#'
#' Generate a class \sQuote{panelPomp} object for fitting to epidemic/endemic Haiti cholera data.
#' All parameters initialized at 0. These values must be overwritten.
#'
#' @param vacscen Vaccination scenario
#' @importFrom pomp Csnippet
#' @importFrom panelPomp panelPomp
#' @return An object of class \sQuote{panelPomp}.
#' @examples
#' m1 <- haiti1(vacscen = 'id0')
#' @export

haiti1_panel <- function(vacscen = 'id0') {
  ## departement list
  depts <- c('Artibonite', 'Centre', 'Grand_Anse', 'Nippes',
             'Nord', 'Nord_Est', 'Nord_Ouest', 'Ouest', 'Sud', 'Sud_Est')
  ## data by departement
  data <- haiti1_data()

  ## list of pomp objects
  pomps <- list()

  # Loop through each departement
  for (dep in depts) {
    # Using haitipkg, create and save pomp object for each departement
    pomps[[dep]] <- haiti1_dep(departement = dep, vacscen = vacscen)
  }

  ## parameters
  shared_pars <- c('mu' = ((1+22.6/1000)^(1/52.14))-1,
                   'delta' = ((1+7.5/1000)^(1/52.14))-1,
                   'gamma' = 7/2,
                   'sigma' = 7/1.4,
                   'theta0' = 0.0,
                   'alpha' = 7/2920,
                   'kappa' = 0) ## KAPPA = 0 FOR NO VACCINATIONS

  shared_pars <- c('mu' = 0, 'delta' = 0, 'gamma' = 0, 'sigma' = 0,
                   'theta0' = 0, 'alpha' = 0, 'kappa' = 0, 'rho' = 0,
                   'tau' = 0, 'beta1' = 0, 'beta2' = 0, 'beta3' = 0,
                   'beta4' = 0, 'beta5' = 0, 'beta6' = 0, 'nu' = 0, 'sig_sq' = 0) ## KAPPA = 0 FOR NO VACCINATIONS

  ## departement specific parameters --- values not determined yet
  spec_par_names <- c("rho", "tau", "beta1", "beta2", "beta3", "beta4", "beta5",
                      "beta6", "nu", "pop_0", "sig_sq",
                      "S_0","E_0","I_0","A_0","R_0",
                      "S1_0","E1_0","I1_0","A1_0","R1_0")

  spec_par_names <- c("pop_0" = 0,
                      "S_0" = 0,"E_0" = 0,"I_0" = 0,"A_0" = 0,"R_0" = 0,
                      "S1_0" = 0,"E1_0" = 0,"I1_0" = 0,"A1_0" = 0,"R1_0" = 0)
#
#   specific_pars <- data.frame(matrix(0, 21, 10))
#   colnames(specific_pars) <- depts
#   rownames(specific_pars) <- spec_par_names


  ## build panelPomp object
  m1_panel <- panelPomp(
    pomps,
    shared = shared_pars,
    specific = spec_par_names
  )

  return(m1_panel)
}
