#' Build panelPomp object for Model 1
#'
#' Generate a class \sQuote{panelPomp} object for fitting to epidemic/endemic Haiti cholera data.
#'
#' @param vacscen Vaccination scenario
#' @importFrom pomp Csnippet
#' @importFrom panelPomp panelPomp
#' @return An object of class \sQuote{panelPomp}.
#' @examples
#' m1 <- haiti1_panel(vacscen = 'id0')
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
    pomps[[dep]] <- haiti1_dep(dept = dep, vacscen = vacscen)
  }

  ## parameters --- from best fit of haiti1_joint()
  shared_pars <- unlist(MODEL1_INPUT_PARAMETERS$panel_shared_params)

  ## departement specific parameters --- values from individual mifs
  spec_par_names <- c("rho_epi", "rho_end",
                      "S_0","E_0","I_0","A_0","R_0", "pop_0",
                      "mob_c")

  pars <- unlist(MODEL1_INPUT_PARAMETERS$dep_params)

  rho_epis <- pars[grepl('rho_epi', names(pars))]
  rho_ends <- pars[grepl('rho_end', names(pars))]
  S0s <- pars[grepl('S_0', names(pars))]
  E0s <- pars[grepl('E_0', names(pars))]
  I0s <- pars[grepl('I_0', names(pars))]
  A0s <- pars[grepl('A_0', names(pars))]
  R0s <- pars[grepl('R_0', names(pars))]
  pop0s <- pars[grepl('pop_0', names(pars))]
  mob_cs <- pars[grepl('mob_c', names(pars))]

  spec_pars <- rbind(rho_epis, rho_ends, S0s, E0s, I0s, A0s, R0s, pop0s, mob_cs)
  colnames(spec_pars) <- depts
  rownames(spec_pars) <- spec_par_names

  ## build panelPomp object
  m1_panel <- panelPomp(
    pomps,
    shared = shared_pars,
    specific = spec_pars
  )

  return(m1_panel)
}
