#' Build pomp object for Model 1 - Disaggregated Joint Model
#'
#' Generate a class \sQuote{pomp} object for fitting to epidemic/endemic Haiti cholera data jointly with separate departments.
#'
#' @param vacscen ID code for vaccination scenario
#' @importFrom pomp Csnippet
#' @return An object of class \sQuote{pomp}
#' @examples
#' m1 <- haiti1_disagg(vacscen = "id0")
#' @export

haiti1_disagg <- function() {
  vacscen <- "id0" # do not need vaccinations
  ## get data
  dat <- haiti1_data() %>%
    dplyr::select(-1)
  colnames(dat) <- c("cases_Artibonite", "cases_Centre", "cases_Grand_Anse",
                     "cases_Nippes", "cases_Nord", "cases_Nord_Est",
                     "cases_Nord_Ouest", "cases_Ouest", "cases_Sud",
                     "cases_Sud_Est", "date_sat", "week")
  fc_set <- vac_scen(vacscen)
  ## make covariate table
  covar <- covars(tmin = 0,
                  tmax = nrow(dat) + 573, ## for 11 year forecast
                  byt = 1,
                  degree = 6,
                  nbasis = 6,
                  per = 52.14,
                  data = dat,
                  settings = fc_set)
  depts <- 10
  depts_names <- c("Artibonite",
                   "Centre",
                   "Grand_Anse",
                   "Nippes",
                   "Nord",
                   "Nord_Est",
                   "Nord_Ouest",
                   "Ouest",
                   "Sud",
                   "Sud_Est")

  ## make components pomp object building
  ## rinit
  state_names_base <- c("S", "E", "I", "A", "R")
  state_names <- c()
  for (i in depts_names) {
    state_names <- c(state_names, paste0(state_names_base, "_", i))
  }
  ivp_names <- paste0(state_names, "_0")
  denoms <- paste0("double denom_", depts_names, " = S_", depts_names,
                   "_0 + ", "E_", depts_names, "_0 + ", "I_", depts_names,
                   "_0 + ", "A_", depts_names, "_0 + ", "R_", depts_names, "_0; \n ")
  denom <- paste(denoms, collapse = "")

  frac_ivp <- paste0(" = nearbyint(frac_", depts_names, " * ", ivp_names, "); \n ")
  state_eqs <- c()
  for (i in 1:length(frac_ivp)) {
    state_eqs <- c(state_eqs, state_names[i], frac_ivp[i])
  }
  state_eqs <- paste(state_eqs, collapse = "")

  rinit_paste <- paste0("double frac_", depts_names, " = pop0_", depts_names,
                        " / (denom_", depts_names, "); \n ",
                        "incid_", depts_names, " = 0.0; \n ",
                        "foival_", depts_names, " = 0.0; \n ",
                        "Str0_", depts_names, " = 0.0; \n ",
                        "Sout_", depts_names, " = 0.0; \n ",
                        "Sin_", depts_names, " = 0.0; \n ")
  rinit_full <- paste(denom, paste(rinit_paste, collapse = ""), state_eqs, collapse = "")

  rinit <- Csnippet(
    rinit_full
  )

  ## rprocess
  # transition rates and numbers
  rates <- rep(c(2, 3, 2, 2, 2), depts)
  trans_rates <- c()
  trans_numbers <- c()
  for (i in 1:length(rates)) {
    trans_rates <- c(trans_rates, paste0("double ", state_names[i], "_rate[", rates[i], "]; \n "))
    trans_numbers <- c(trans_numbers, paste0("double ", state_names[i], "_trans[", rates[i], "]; \n "))
  }
  trans_rates <- paste(trans_rates, collapse = "")
  trans_numbers <- paste(trans_numbers, collapse = "")
  trans_numbers <- paste0(trans_numbers)


  # demonitors and time checks
  demons <- c(paste0("int pop_", depts_names, " = S_", depts_names, " + E_",
                     depts_names, " + I_", depts_names, " + A_", depts_names,
                     " + R_", depts_names, "; \n "),
              paste0("int births_", depts_names, " = rpois(mu*pop_", depts_names, "*dt); \n")) %>%
    paste(collapse = "")

  # seasonal beta and foi
  time_check <- paste0("double mybeta_", depts_names,
                       " = beta1_", depts_names, "*seas1 + beta2_",
                       depts_names, "*seas2 + beta3_", depts_names,
                       "*seas3 + beta4_", depts_names, "*seas4 + beta5_",
                       depts_names, "*seas5 + beta6_", depts_names, "*seas6; \n "
  )
  beta <- paste0(time_check, collapse = "")

  # get other_cases for each department
  all_dat <- haiti1_data() %>%
    dplyr::select(-1) %>%
    dplyr::rename(Grand_Anse = Grand.Anse,
                  Nord_Est = Nord.Est,
                  Nord_Ouest = Nord.Ouest,
                  Sud_Est = Sud.Est)
  other_cases <- all_dat %>%
    dplyr::select(week)
  all_dat <- all_dat %>%
    dplyr::ungroup() %>%
    dplyr::select(-c(week, date_sat))
  for (i in depts_names) {
    other_depts <- all_dat %>%
      dplyr::select(-c(i))
    others <- rowSums(other_depts)
    other_cases[, i] <- others
  }

  mobs <- paste0("double mobility_", depts_names, " = 0.0; \n ",
                 "oth_case = other_cases_mat[t_val][", 1:depts, "]; \n ",
                 "if(oth_case < 0) { \n mobility_", depts_names,
                 " = 0.0; \n } else { \n mobility_", depts_names, " = mob_c_",
                 depts_names, "*oth_case; } \n")
  fois <- paste0("double foi_", depts_names, " = pow(I_", depts_names, " + (1-kappa)*A_",
                 depts_names, ", nu_", depts_names, ")*mybeta_", depts_names, "/pop_",
                 depts_names, " + mobility_", depts_names, "; \n")

  other_cases_string <- foreach(
    r = iterators::iter(other_cases, by = "row"),
    .combine = c
  ) %do% {
    sprintf(" {%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f} ",
            r$week, r$Artibonite, r$Centre, r$`Grand_Anse`, r$`Nippes`,
            r$`Nord`, r$`Nord_Est`, r$`Nord_Ouest`, r$Ouest, r$Sud, r$`Sud_Est`)
  } %>%
    stringr::str_c(collapse = ", \n")

  matrix_other_cases <- stringr::str_c(
    sprintf("double other_cases_mat[%i][%i] = {\n", nrow(other_cases), ncol(other_cases)),
    other_cases_string,
    " \n };"
  )

  matrix_other_cases <- stringr::str_replace_all(matrix_other_cases, "NA", "-1")

  foi <- paste(mobs, fois, paste0("double sig_sq_", depts_names, "; \n "),
               "if (t < 233) { \n ",
               paste0("sig_sq_", depts_names, " = sig_sq_epi_", depts_names, "; \n "),
               "} else { \n ",
               paste0("sig_sq_", depts_names, " = sig_sq_end_", depts_names, "; \n "),
               "} \n ",
               paste0("\n double dgamma_", depts_names, " = rgammawn(sig_sq_", depts_names,
                      ", dt); \n foi_", depts_names, " = foi_", depts_names, "* dgamma_",
                      depts_names, "/dt; \n "), collapse = "")

  # compute rates
  s_rates <- paste0("S_", depts_names, "_rate[0] = foi_", depts_names,
                    "; \n S_", depts_names, "_rate[1] = delta; \n ")
  e_rates <- paste0("E_", depts_names, "_rate[0] = sigma*(1-theta0); \n E_",
                    depts_names, "_rate[1] = sigma*theta0; \n E_", depts_names, "_rate[2] = delta; \n ")
  i_rates <- paste0("I_", depts_names, "_rate[0] = gamma; \n I_", depts_names, "_rate[1] = delta; \n ")
  a_rates <- paste0("A_", depts_names, "_rate[0] = gamma; \n A_", depts_names, "_rate[1] = delta; \n ")
  r_rates <- paste0("R_", depts_names, "_rate[0] = alpha; \n R_", depts_names, "_rate[1] = delta; \n ")
  rates <- c(s_rates, e_rates, i_rates, a_rates, r_rates) %>%
    paste(collapse = "")

  # transition numbers
  numbers <- c(paste0("reulermultinom(2,S_", depts_names, ",&S_",
                      depts_names, "_rate[0],dt,&S_", depts_names, "_trans[0]); \n "),
               paste0("reulermultinom(3,E_", depts_names, ",&E_",
                      depts_names, "_rate[0],dt,&E_", depts_names, "_trans[0]); \n "),
               paste0("reulermultinom(2,I_", depts_names, ",&I_",
                      depts_names, "_rate[0],dt,&I_", depts_names, "_trans[0]); \n "),
               paste0("reulermultinom(2,A_", depts_names, ",&A_",
                      depts_names, "_rate[0],dt,&A_", depts_names, "_trans[0]); \n "),
               paste0("reulermultinom(2,R_", depts_names, ",&R_",
                      depts_names, "_rate[0],dt,&R_", depts_names, "_trans[0]); \n ")) %>%
    paste(collapse = "")

  # cohorts
  coh <- c(paste0("S_", depts_names, " += -S_", depts_names, "_trans[0] - S_", depts_names,
                  "_trans[1] + R_", depts_names, "_trans[0]; \n "),
           paste0("E_", depts_names, " += -E_", depts_names, "_trans[0] - E_", depts_names, "_trans[1] - E_",
                  depts_names, "_trans[2] + S_", depts_names, "_trans[0]; \n "),
           paste0("I_", depts_names, " += -I_", depts_names, "_trans[0] - I_", depts_names,
                  "_trans[1] + E_", depts_names, "_trans[0]; \n "),
           paste0("A_", depts_names, " += -A_", depts_names, "_trans[0] - A_", depts_names,
                  "_trans[1] + E_", depts_names, "_trans[0]; \n "),
           paste0("R_", depts_names, " += -R_", depts_names, "_trans[0] - R_", depts_names,
                  "_trans[1] + A_", depts_names, "_trans[0] + I_", depts_names, "_trans[0]; \n "))
  coh <- coh %>%
    paste(collapse = "")

  ## accumvars
  incids <- paste0("incid_", depts_names, " += E_", depts_names, "; \n ")
  fois <- paste0("foival_", depts_names, " += foi_", depts_names, "; \n ")
  str0 <- paste0("Str0_", depts_names, "+= S_", depts_names, "_trans[0]; \n ")
  sin <- paste0("Sin_", depts_names, " += R_", depts_names, "_trans[0] + births_", depts_names, "; \n ")
  sout <- paste0("Sout_", depts_names, " += S_", depts_names, "_trans[0] + S_", depts_names, "_trans[1]; \n ")
  last <- c(fois, str0, sin, sout, incids) %>%
    paste(collapse = "")

  rproc_paste <- c(matrix_other_cases,
                   "int t_val = (int)t; \n ",
                   "double oth_case = 0.0; \n ",
                   trans_rates, trans_numbers, demons, beta, foi, rates, numbers, coh, last) %>%
    paste(collapse = "")

  rproc <- Csnippet(
    rproc_paste
  )

  ## dmeasure
  liks <- paste0("lik = 0.0; \n ",
                 "double rho_", depts_names, " = 0.0; \n ",
                 "double tau_", depts_names, " = 0.0; \n ",
                 "if (ISNA(cases_", depts_names, ")) { \n ",
                 "lik += (give_log) ? 0 : 1; \n ",
                 "} else { \n ",
                 "rho_", depts_names, " = rho_epi_", depts_names, "; \n ",
                 "tau_", depts_names, " = tau_epi_", depts_names, "; \n ",
                 "if (t > 232) { \n ",
                 "rho_", depts_names, " = rho_end_", depts_names, "; \n ",
                 "tau_", depts_names, " = tau_end_", depts_names, "; \n",
                 "} \n ",
                 "lik += dnbinom_mu(cases_", depts_names, ", tau_", depts_names, ", rho_", depts_names, "*incid_", depts_names, ", give_log); \n ",
                 "} \n ")

  dmeas <- Csnippet(
    liks
  )

  ## rmeasure
  cases <- paste0("cases_", depts_names, " = 0.0; \n ",
                  "double rho_", depts_names, " = rho_epi_", depts_names, "; \n ",
                  "double tau_", depts_names, " = tau_epi_", depts_names, "; \n ",
                  "if (t > 232) { \n ",
                  "rho_", depts_names, " = rho_end_", depts_names, "; \n ",
                  "tau_", depts_names, " = tau_end_", depts_names, "; \n",
                  "} \n ",
                  "cases_", depts_names, " = rnbinom_mu(tau_", depts_names, ", rho_", depts_names, "*incid_", depts_names, "); \n ",
                  "if (cases_", depts_names, " > 0.0) { \n ",
                  "cases_", depts_names, " = nearbyint(cases_", depts_names, "); \n ",
                  "} else { \n",
                  "cases_", depts_names, " = 0.0; \n ",
                  "} \n ")

  rmeas <- Csnippet(
    cases
  )

  ## names
  ## state names
  state_names <- c(paste0("S_", depts_names),
                   paste0("E_", depts_names),
                   paste0("I_", depts_names),
                   paste0("A_", depts_names),
                   paste0("R_", depts_names),
                   paste0("incid_", depts_names),
                   paste0("foival_", depts_names),
                   paste0("Str0_", depts_names),
                   paste0("Sout_", depts_names),
                   paste0("Sin_", depts_names))

  ## accum vars
  accum_names <- c(paste0("incid_", depts_names),
                   paste0("foival_", depts_names),
                   paste0("Str0_", depts_names),
                   paste0("Sout_", depts_names),
                   paste0("Sin_", depts_names)) %>%
    paste(collapse = "")


  ## partrans
  param_trans <- pomp::parameter_trans(
    log = c(paste0("beta1_", depts_names),
            paste0("beta2_", depts_names),
            paste0("beta3_", depts_names),
            paste0("beta4_", depts_names),
            paste0("beta5_", depts_names),
            paste0("beta6_", depts_names),
            paste0("sig_sq_end_", depts_names),
            paste0("sig_sq_epi_", depts_names),
            paste0("tau_end_", depts_names),
            paste0("tau_epi_", depts_names),
            "sigma", "gamma", "mu", "delta", "alpha"),
    logit = c(paste0("rho_epi_", depts_names),
              paste0("rho_end_", depts_names),
              paste0("nu_", depts_names), "theta0"),
    barycentric = c(paste0("S_", depts_names, "_0"),
                    paste0("E_", depts_names, "_0"),
                    paste0("I_", depts_names, "_0"),
                    paste0("A_", depts_names, "_0"),
                    paste0("R_", depts_names, "_0"))
  )

  ## hand entered for now
  pars <- unlist(MODEL1_INPUT_PARAMETERS$dep_params)

  rho_epis <- paste0("rho_epi_", depts_names)
  rho_ends <- paste0("rho_end_", depts_names)
  tau_epis <- paste0("tau_epi_", depts_names)
  tau_ends <- paste0("tau_end_", depts_names)
  sig_sq_epis <- paste0("sig_sq_epi_", depts_names)
  sig_sq_ends <- paste0("sig_sq_end_", depts_names)
  beta1s <- paste0("beta1_", depts_names)
  beta2s <- paste0("beta2_", depts_names)
  beta3s <- paste0("beta3_", depts_names)
  beta4s <- paste0("beta4_", depts_names)
  beta5s <- paste0("beta5_", depts_names)
  beta6s <- paste0("beta6_", depts_names)
  nus <- paste0("nu_", depts_names)
  S0s <- paste0("S_", depts_names, "_0")
  E0s <- paste0("E_", depts_names, "_0")
  I0s <- paste0("I_", depts_names, "_0")
  A0s <- paste0("A_", depts_names, "_0")
  R0s <- paste0("R_", depts_names, "_0")
  pop0s <- paste0("pop0_", depts_names)
  mob_cs <- paste0("mob_c_", depts_names)

  par_names <- c(rho_epis, rho_ends, tau_epis, tau_ends, sig_sq_epis,
                 sig_sq_ends, beta1s, beta2s, beta3s, beta4s, beta5s, beta6s,
                 nus, S0s, E0s, I0s, A0s, R0s, pop0s, mob_cs,
                 "gamma", "sigma", "theta0", "alpha", "mu", "delta", "kappa")

  rho_epis <- pars[grepl('rho_epi', names(pars))]
  rho_ends <- pars[grepl('rho_end', names(pars))]
  tau_epis <- pars[grepl('tau_epi', names(pars))]
  tau_ends <- pars[grepl('tau_end', names(pars))]
  sig_sq_epis <- pars[grepl('sig_sq_epi', names(pars))]
  sig_sq_ends <- pars[grepl('sig_sq_end', names(pars))]
  beta1s <- pars[grepl('beta1', names(pars))]
  beta2s <- pars[grepl('beta2', names(pars))]
  beta3s <- pars[grepl('beta3', names(pars))]
  beta4s <- pars[grepl('beta4', names(pars))]
  beta5s <- pars[grepl('beta5', names(pars))]
  beta6s <- pars[grepl('beta6', names(pars))]
  nus <- pars[grepl('nu', names(pars))]
  S0s <- pars[grepl('S_0', names(pars))]
  E0s <- pars[grepl('E_0', names(pars))]
  I0s <- pars[grepl('I_0', names(pars))]
  A0s <- pars[grepl('A_0', names(pars))]
  R0s <- pars[grepl('R_0', names(pars))]
  pop0s <- pars[grepl('pop_0', names(pars))]
  mob_cs <- pars[grepl('mob_c', names(pars))]

  par_vals <- c(rho_epis, rho_ends, tau_epis, tau_ends, sig_sq_epis,
                sig_sq_ends, beta1s, beta2s, beta3s, beta4s, beta5s, beta6s,
                nus, S0s, E0s, I0s, A0s, R0s, pop0s, mob_cs,
                3.5, 5.0, 0.0, 0.00239726, 0.0004287149, 0.000143317, 0.0) %>%
    unlist()

  params <- par_vals
  names(params) <- par_names

  ## build pomp model
  model1 <- pomp::pomp(
    data = dat,
    times = "week",
    t0 = 0,
    rmeasure = rmeas,
    dmeasure = dmeas,
    rprocess = pomp::euler(step.fun = rproc, delta.t = 1/7),
    covar = pomp::covariate_table(covar, times = "time"),
    partrans = param_trans,
    statenames = state_names,
    paramnames = par_names,
    params = params,
    accumvars = accum_names,
    rinit = rinit
  )

  return(model1)
}
