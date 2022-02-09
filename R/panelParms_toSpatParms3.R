#' Panel to Spat parameters
#'
#' This function converts the parameters fit for haiti3_panel to parameters
#' that work for haiti3_spatPomp. This is done because the parameters are
#' fit using panel_pomp, but ultimately the coupling requires the use of a
#' spatPomp object in order to use the block particle filter.
#'
#' @param panel_parms Any set of parameters for a panel_pomp version of
#'  model 3. Typically, these parameters will be the parameters that resulted
#'  in the best fit for the model.
#'
#' @import pomp
#' @import spatPomp

panelParms_toSpatParms3 <- function(panel_parms) {
  mod_coupled <- haiti3_spatPomp(dt_years = 1/365.25)
  panel_parms <- panel_parms[!names(panel_parms) %in% c("t_vacc_start",
                                                        "t_vacc_end",
                                                        "p1d_reg",
                                                        "r_v_year")]

  departements <-
    c(
      'Artibonite',
      'Centre',
      'Grande_Anse',
      'Nippes',
      'Nord',
      'Nord-Est',
      'Nord-Ouest',
      'Ouest',
      'Sud',
      'Sud-Est'
    )

  params_common <- c(
    "sigma", "mu_B", "thetaI", "XthetaA", "lambdaR", "r",
    "gamma", "rho", "epsilon", "k",
    "std_W", "cas_def", "mu", "alpha", "cases_ext"
  )

  spatPompCoefs <- coef(mod_coupled)

  # Set the spatPompCoefs with departement spacific values
  for (i in 1:10) {
    dp <- departements[i]
    new_params <- best_params[grepl(paste0("\\[", dp, "\\]"), names(best_params))]
    names(new_params) <- gsub(paste0("\\[", dp, "\\]"), i, names(new_params))
    which.b0 <- grepl("^B0[[:digit:]]{1,2}$", names(new_params))
    names(new_params)[which.b0] <- paste0("Binit", i)
    spatPompCoefs[names(new_params)] <- new_params
  }

  # Set the spatPompCoefs with the shared parameter values
  for (parm in params_common) {
    spatPompCoefs[paste0(parm, 1:10)] <- best_params[parm]
  }

  spatPompCoefs

}
