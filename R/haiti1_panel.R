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

  ## parameters --- from best fit of haiti1_joint()
  shared_pars <- c('mu' = ((1+22.6/1000)^(1/52.14))-1,
                   'delta' = ((1+7.5/1000)^(1/52.14))-1,
                   'gamma' = 7/2,
                   'sigma' = 7/1.4,
                   'theta0' = 0.0,
                   'alpha' = 7/2920,
                   'kappa' = 0,
                   'beta1' = 4.014758,
                   'beta2' = 2.7089,
                   'beta3' = 2.742331,
                   'beta4' = 3.058927,
                   'beta5' = 3.57466,
                   'beta6' = 2.230872,
                   'nu' = 0.9976078)

  ## departement specific parameters --- values from individual mifs
  spec_par_names <- c("rho_epi", "rho_end",
                      "tau_epi", "tau_end",
                      "sig_sq_epi", "sig_sq_end",
                      "S_0","E_0","I_0","A_0","R_0", "pop_0")

  dep_params <- haiti1_dep_params
  dep_params <- dep_params[rownames(dep_params) %in% spec_par_names, ]

  spec_pars <- rbind(
    unlist(dep_params["rho_epi", ]),
    unlist(dep_params["rho_end", ]),
    unlist(dep_params["tau_epi", ]),
    unlist(dep_params["tau_end", ]),
    unlist(dep_params["sig_sq_epi", ]),
    unlist(dep_params["sig_sq_end", ]),
    unlist(dep_params["S_0", ]),
    unlist(dep_params["E_0", ]),
    unlist(dep_params["I_0", ]),
    unlist(dep_params["A_0", ]),
    unlist(dep_params["R_0", ]),
    unlist(dep_params["pop_0", ]))

  rownames(spec_pars) <- spec_par_names

  ## build panelPomp object
  m1_panel <- panelPomp(
    pomps,
    shared = shared_pars,
    specific = spec_par_names
  )

  return(m1_panel)
}
