#' Build PanelPomp representation of model 3. The shared parameters are
#' set to the best fit on data from Artibonite, and the unit specific
#' parameters are set to their best independent values using the values
#' fit in Artibonite.
#'
#' @param delta.t delta time step used in the model
#' @param departements vector of departements to include in the panelPomp.
#'    The default is all departements.
#' @param BetaB_trend Boolean indicator of whether or not a trend parameter
#'    should be added to the betaB parameter.
#' @param start_time Time of which to start the panel. In the
#'    original paper, the authors used a start time of 2014-03-01 for all
#'    departements except for Ouest, which started with a time of 2017-06-10.
#'    This means that there wasn't a national-coupled model fit to the
#'    data until 2017-06-01, which means that there was less than 2 years
#'    of data to fit the model.
#' @param B0 Boolean indicator. If TRUE, Initial value of the Bacteria
#'    compartment will be estimated; If FALSE, the equilibrium will be used.
#' @importFrom panelPomp panelPomp
#' @export

haiti3_panel <- function(delta.t = 1/365, departements = c(
  'Artibonite', 'Centre', 'Grande_Anse', 'Nippes',
  'Nord', 'Nord-Est', 'Nord-Ouest', 'Ouest', 'Sud', 'Sud-Est'),
  betaB_trend = FALSE, start_time = "2014-03-01", B0 = FALSE,
  old_mod = FALSE
) {

  # Create a list of pomp objects
  pomps <- list()

  if (old_mod) {
    # Best values fitted using Artibonite data, hand entered.
    SHARED <- c(
      'mu_B' = 111.3593, 'XthetaA' = 0.05779537, 'thetaI' = 0.0004327653,
      'lambdaR' = 0.499169, 'std_W' = 0.01400412, 'epsilon' = 0.993717,
      'k' = 516.8755, 'cas_def' = 0.9341205, 'sigma' = 0.25, 'r' = 0.6359411,
      'gammaI' = 182.625,
      'gammaA' = 182.625,
      'rhoA' = 0.1250856,
      'XrhoI' = 1,
      'Rtot_0' = 0.35, 'mu' = 0.01586626, 'alpha' = 1.461, 'cases_ext' = 1,
      't_vacc_start' = 0, 't_vacc_end' = 0, 'p1d_reg' = 0, 'r_v_year' = 0
    )
  } else {
    # Best values fitted using Artibonite data, hand entered.
    SHARED <- c(
      'mu_B' = 111.3593, 'XthetaA' = 0.05779537, 'thetaI' = 0.0004327653,
      'lambdaR' = 0.499169, 'std_W' = 0.01400412, 'epsilon' = 0.993717,
      'k' = 516.8755, 'cas_def' = 0.9341205, 'sigma' = 0.25, 'r' = 0.6359411,
      'gamma' = 182.625,
      'rho' = 0.1250856,
      'Rtot_0' = 0.35, 'mu' = 0.01586626, 'alpha' = 1.461, 'cases_ext' = 1,
      't_vacc_start' = 0, 't_vacc_end' = 0, 'p1d_reg' = 0, 'r_v_year' = 0
    )
  }

  all_deps <- c(
    'Artibonite', 'Centre', 'Grande_Anse', 'Nippes',
    'Nord', 'Nord-Est', 'Nord-Ouest', 'Ouest', 'Sud', 'Sud-Est'
  )

  specific_param_names <- c('betaB', 'foi_add', 'H', 'D')

  # Best values of individual departement fits
  SPECIFIC <- rbind(
    c(0.354169973579405, 1.38491428682277, 1.91956651003771,
      2.0582526120949, 0.403718310499447, 2.23586092140192,
      0.774580862027951, 0.079403943457212, 0.905802818696324,
      0.973571920059841),
    c(1.18214116603088e-06, 5.15295170410547e-06, 5.15701684501422e-07,
      2.34573186547564e-07, 4.89222508562004e-07, 6.05955895165197e-07,
      6.52134987324786e-07, 1.09575724078116e-07, 5.2034928078208e-07,
      6.14770089555291e-07),
    c(1727524, 746236, 468301, 342525, 1067177,
      393967, 728807, 4029705, 774976, 632601),
    c(353, 214, 245, 270, 505, 243, 347, 809, 292, 311)
  )

  if (betaB_trend) {
    SPECIFIC <- rbind(SPECIFIC, rep(0, 10))
    colnames(SPECIFIC) <- all_deps
    rownames(SPECIFIC) <- c(specific_param_names, "betaB_trend")
  } else if (B0) {
    SPECIFIC <- rbind(SPECIFIC, rep(0.024, 10))
    colnames(SPECIFIC) <- all_deps
    rownames(SPECIFIC) <- c(specific_param_names, "B0")
  } else {
    colnames(SPECIFIC) <- all_deps
    rownames(SPECIFIC) <- specific_param_names
  }

  # Loop through each departement
  for (dep in departements) {

    if (old_mod) {
      # Using haitipkg, create and save pomp object for each departement
      pomps[[dep]] <- haiti3_dep(departement = dep, delta.t = delta.t,
                                 betaB_trend = betaB_trend,
                                 start_time = start_time,
                                 B0 = B0)
    } else {
      # Using haitipkg, create and save pomp object for each departement
      pomps[[dep]] <- haiti3_dep_correct(
        departement = dep,
        delta.t = delta.t,
        start_time = start_time
      )
    }
  }

  SPECIFIC <- SPECIFIC[, departements]

  # Create the panelPomp object
  SIRB_panel <- panelPomp(
    pomps,
    shared = SHARED,
    specific = SPECIFIC
  )

  SIRB_panel
}
