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
  dat <- haiti1_data()
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
  denom <- ivp_names[1]
  for (i in 2:length(ivp_names)) {
    denom <- c(denom, " + ", ivp_names[i])
  }
  denom <- paste(denom, collapse = "")
  frac_ivp <- paste0(" = nearbyint(frac * ", ivp_names, "); \n ")
  state_eqs <- c()
  for (i in 1:length(frac_ivp)) {
    state_eqs <- c(state_eqs, state_names[i], frac_ivp[i])
  }
  state_eqs <- paste(state_eqs, collapse = "")

  rinit_paste <- paste0("double frac = pop_0 / (", denom, "); \n ",
                        state_eqs,
                        "incid = 0.0; \n ",
                        "foival = 0.0; \n ",
                        "Str0 = 0.0; \n ",
                        "Sout = 0.0; \n ",
                        "Sin = 0.0; \n ")
  rinit_paste <- c(rinit_paste, "incidU = 0.0; \n ", "incidV = 0.0; \n",
                     "asymV = 0.0; \n ", "newV = 0.0; \n ")

  rinit <- Csnippet(
    paste(rinit_paste, collapse = "")
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
  trans_numbers <- paste0(trans_numbers, "\n double dgamma; \n ")


  # demonitors and time checks
  demons <- c(paste0("int pop_", depts_names, " = S_", depts_names, " + E_",
                     depts_names, " + I_", depts_names, " + A_", depts_names,
                     " + R_", depts_names, "; \n "),
              paste0("int birts_", depts_names, " = rpois(mu*pop_", depts_names, "*dt); \n")) %>%
    paste(collapse = "")

  # seasonal beta and foi
  time_check <- c(
    "double mybeta = beta1*seas1 + beta2*seas2 + beta3*seas3 + beta4*seas4 + beta5*seas5 + beta6*seas6; \n "
  )
  beta <- paste0(time_check, collapse = "")

  # get other_cases for each department
  all_dat <- haiti1_data()
  other_cases <- all_dat %>%
    dplyr::select(week)
  all_dat <- all_dat %>%
    dplyr::select(-c(week, date_sat))
  for (i in depts_names) {
    other_depts <- all_dat %>%
      dplyr::select(-i)
    others <- rowSums(other_depts)
    other_cases <- cbind(other_cases, others)
  }
  colnames(other_cases) <- c("week",
                             depts_names)

  mobs <- paste0("double mobility_", depts_names, " = if(ISNA(other_cases[t, ", depts_names, "])) { \n mobility = 0.0; \n } else { \n mobility = mob_c_", depts_names, "*other_cases[t, ", depts_names, "]; } \n")
  fois <- paste0("double foi_", depts_names, " = pow(I_", depts_names, " + (1-kappa) * A_", depts_names, ", nu)*mybeta/pop_", depts_names, " + mobility_", depts_names, "; \n")
  foi <- paste0(mobs,
                fois,
                "double sig_sq; \n ",
                "if (t < 233) { \n ",
                paste0("sig_sq_", depts_names, " = sig_sq_epi_", depts_names, "; \n "),
                "} else { \n ",
                paste0("sig_sq_", depts_names, " = sig_sq_end_", depts_names, "; \n "),
                "} \n ",
                paste0("\n dgamma_", depts_names, " = rgammawn(sig_sq_", depts_names, ", dt); \n foi_", depts_names, " = foi_", depts_names, "* dgamma_", depts_names, "/dt; \n "))

  # compute rates
  s_rates <- paste0("S_", depts_names, "_rate[0] = foi_", depts_names, "; \n S_", depts_names, "_rate[1] = delta; \n ")
  e_rates <- paste0("E_", depts_names, "_rate[0] = sigma*(1-theta0); \n E_", depts_names, "_rate[1] = sigma*theta0; \n E", depts_names, "_rate[2] = delta; \n ")
  i_rates <- paste0("I_", depts_names, "_rate[0] = gamma; \n I_", depts_names, "_rate[1] = delta; \n ")
  a_rates <- paste0("A_", depts_names, "_rate[0] = gamma; \n A_", depts_names, "_rate[1] = delta; \n ")
  r_rates <- paste0("R_", depts_names, "_rate[0] = alpha; \n R_", depts_names, "_rate[1] = delta; \n ")
  rates <- c(s_rates, e_rates, i_rates, a_rates, r_rates) %>%
    paste(collapse = "")

  # transition numbers
  numbers <- c(paste0("reulermultinom(2,S_", depts_names, ",&S_", depts_names, "_rate[0],dt,&S_", depts_names, "_trans[0]); \n "),
               paste0("reulermultinom(3,E_", depts_names, ",&E_", depts_names, "_rate[0],dt,&E_", depts_names, "_trans[0]); \n "),
               paste0("reulermultinom(2,I_", depts_names, ",&I_", depts_names, "_rate[0],dt,&I_", depts_names, "_trans[0]); \n "),
               paste0("reulermultinom(2,A_", depts_names, ",&A_", depts_names, "_rate[0],dt,&A_", depts_names, "_trans[0]); \n "),
               paste0("reulermultinom(2,R_", depts_names, ",&R_", depts_names, "_rate[0],dt,&R_", depts_names, "_trans[0]); \n ")) %>%
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

  ################# unfinished past here

  # incidence
  incids <- c("incid += Etrans[0]; \n ")
  foi_val <- "foival += foi; \n "
  str0 <- "Str0 += Strans[0]; \n "
  sin <- c("Sin += Rtrans[0] + births; \n ")
  sout <- c("Sout += Strans[0] + Strans[1]; \n ")
  last <- c(foi_val, str0, sin)

  if (depts > 1) {
    incids <- c("incid += Etrans[0]", paste0(" + E", 1:depts, "trans[0]"), "; \n ")
    incid_u <- "incidU += Etrans[0]; \n "
    incid_v <- c("incidV += E1trans[0]", paste0(" + E", 2:depts, "trans[0]"), "; \n ")
    asym_v <- c("asymV += E1trans[1]", paste0(" + E", 2:depts, "trans[1]"), "; \n ")
    new_v <- c("newV += Strans[2]", paste0(" + Strans[", 3:(depts + 1), "]"),
               paste0(" + Etrans[", 3:(depts + 2), "]"),
               paste0(" + Itrans[", 2:(depts + 1), "]"),
               paste0(" + Atrans[", 2:(depts + 1), "]"),
               paste0(" + Rtrans[", 2:(depts + 1), "]"), "; \n ")
    sout <- c("Sout += Strans[0]", paste0(" + Strans[", 1:(depts + 1), "]"), "; \n ")
    last <- c(last, incid_u, incid_v, asym_v, new_v)
  }
  last <- c(last, incids, sout) %>%
    paste(collapse = "")

  if (depts > 1) {
    rproc_paste <- c(trans_rates, trans_numbers, vac_rates, demons, tchecks, beta,
                     foi, thetas, rates, numbers, coh, last) %>%
      paste(collapse = "")
  } else {
    rproc_paste <- c(trans_rates, trans_numbers, demons, beta, foi, rates, numbers,
                     coh, last) %>%
      paste(collapse = "")
  }

  rproc <- Csnippet(
    rproc_paste
  )

  ## dmeasure
  dmeas <- Csnippet("
    if (ISNA(cases)) {
      lik = (give_log) ? 0 : 1;
    } else {
      double rho = rho_epi;
      double tau = tau_epi;
      if (t > 232) {
        rho = rho_end;
        tau = tau_end;
      }
      lik = dnbinom_mu(cases, tau, rho*incid, give_log);
    }
  ")

  ## rmeasure
  rmeas <- Csnippet("
    double rho = rho_epi;
    double tau = tau_epi;
    if (t > 232) {
      rho = rho_end;
      tau = tau_end;
    }
    cases = rnbinom_mu(tau, rho*incid);
    if (cases > 0.0) {
      cases = nearbyint(cases);
    } else {
      cases = 0.0;
    }
  ")

  ## names
  if (depts > 1) {
    ## state names
    state_names <- c("S", paste0("S", 1:depts),
                     "E", paste0("E", 1:depts),
                     "I", paste0("I", 1:depts),
                     "A", paste0("A", 1:depts),
                     "R", paste0("R", 1:depts),
                     "incid", "incidU", "incidV", "asymV", "newV",
                     "foival", "Str0", "Sout", "Sin")

    ## parameter names
    param_names <- c("rho_epi", "sig_sq_epi", "tau_epi", #epidemic
                     "S_0", "E_0", "I_0", "A_0", "R_0", "pop_0", # epidemic
                     "rho_end", "sig_sq_end", "tau_end", # endemic
                     "beta1", "beta2", "beta3", "beta4", "beta5", "beta6", # shared
                     "mu", "gamma", "sigma", "theta0", "alpha", "delta", "kappa", "nu",# shared
                     paste0("S", 1:depts, "_0"),
                     paste0("E", 1:depts, "_0"),
                     paste0("I", 1:depts, "_0"),
                     paste0("A", 1:depts, "_0"),
                     paste0("R", 1:depts, "_0"))

    ## accum vars
    accum_names <- c("incid","incidU","incidV","asymV","newV",
                     "foival","Str0","Sout","Sin")
  } else {
    ## state names
    state_names <- c("S", "E", "I", "A", "R", "incid", "foival", "Str0", "Sout", "Sin")

    ## parameter names
    param_names <- c("rho_epi", "sig_sq_epi", "tau_epi", #epidemic
                     "S_0", "E_0", "I_0", "A_0", "R_0", "pop_0", # epidemic
                     "rho_end", "sig_sq_end", "tau_end", # endemic
                     "beta1", "beta2", "beta3", "beta4", "beta5", "beta6", # shared
                     "mu", "gamma", "sigma", "theta0", "alpha", "delta", "kappa", "nu") # shared

    ## accum vars
    accum_names <- c("incid", "foival","Str0","Sout","Sin")
  }

  ## partrans
  param_trans <- pomp::parameter_trans(
    log = c("beta1", "beta2", "beta3", "beta4", "beta5", "beta6", # shared
            "sigma", "gamma", "mu", "delta", "alpha",
            "sig_sq_end", "sig_sq_epi", "tau_end", "tau_epi"),
    logit = c("rho_epi", "rho_end", "nu", "theta0"),
    barycentric = c("S_0", "E_0", "I_0", "A_0", "R_0")
  )

  ## hand entered for now
  pars <- c("rho_epi" = 0.4765437,
            "rho_end" = 0.4496893,
            "tau_epi" = 688.7796,
            "tau_end" = 105.3583,
            "sig_sq_epi" = 0.1105648,
            "sig_sq_end" = 0.1677307,
            "beta1" = 4.014758,
            "beta2" = 2.7089,
            "beta3" = 2.742331,
            "beta4" = 3.058927,
            "beta5" = 3.57466,
            "beta6" = 2.230872,
            "nu" = 0.9976078,
            "gamma" = 3.5,
            "sigma" = 5.0,
            "theta0" = 0.0,
            "alpha" = 0.00239726,
            "mu" = 0.0004287149,
            "delta" = 0.000143317,
            "kappa" = 0.0,
            "S_0" = 0.9990317,
            "E_0" = 4.604823e-06,
            "I_0" = 0.000963733,
            "A_0" = 0.0,
            "R_0" = 0.0,
            "pop_0" = 10911819)

  if (depts > 1) {
    par_names <- names(pars)
    pars <- c(pars, rep(0.0, 5 * depts))
    names(pars) <- c(par_names,
                     paste0("S", 1:depts, "_0"),
                     paste0("E", 1:depts, "_0"),
                     paste0("I", 1:depts, "_0"),
                     paste0("A", 1:depts, "_0"),
                     paste0("R", 1:depts, "_0"))
  }

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
    paramnames = param_names,
    params = pars,
    accumvars = accum_names,
    rinit = rinit
  )

  return(model1)
}
