#' Build pomp object for Model 1 - Joint
#'
#' Generate a class \sQuote{pomp} object for fitting to epidemic/endemic Haiti cholera data jointly.
#'
#' @param vacscen ID code for vaccination scenario
#' @importFrom pomp Csnippet
#' @return An object of class \sQuote{pomp}
#' @examples
#' m1 <- haiti1_joint(vacscen = "id0")
#' @export

haiti1_joint <- function(vacscen = 'id0') {
  vacscen <- vacscen
  ## get data
  dat <- haiti1_agg_data()
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
  if (vacscen == 'id0') {
    depts <- 0
  } else {
    depts <- fc_set$nd
  }

  ## make components pomp object building
  ## rinit
  state_names_base <- c("S", "E", "I", "A", "R")
  state_names <- state_names_base
  if (depts > 1) {
    for (i in 1:depts) {
      state_names <- c(state_names, paste0(state_names_base, i))
    }
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
  if (depts > 1) {
    rinit_paste <- c(rinit_paste, "incidU = 0.0; \n ", "incidV = 0.0; \n",
                     "asymV = 0.0; \n ", "newV = 0.0; \n ")
  }

  rinit <- Csnippet(
    paste(rinit_paste, collapse = "")
  )

  ## rprocess
  # transition rates and numbers
  rates_base <- c((2 + depts), (3 + depts), (2 + depts), (2 + depts), (2 + depts)) ## for SEIAR
  if (depts > 1) {
    rates_other <- rep(c(2, 3, 2, 2, 2), depts)
    rates <- c(rates_base, rates_other)
  } else {
    rates <- rates_base
  }
  trans_rates <- c()
  trans_numbers <- c()
  for (i in 1:length(rates)) {
    trans_rates <- c(trans_rates, paste0("double ", state_names[i], "rate[", rates[i], "]; \n "))
    trans_numbers <- c(trans_numbers, paste0("double ", state_names[i], "trans[", rates[i], "]; \n "))
  }
  trans_rates <- paste(trans_rates, collapse = "")
  trans_numbers <- paste(trans_numbers, collapse = "")
  trans_numbers <- paste0(trans_numbers, "\n double dgamma; \n ")

  # vaccination rates
  if (depts > 1) {
    vac_rates <- paste0("double eta", 1:depts, " = 0.0; \n ")
  }

  # demonitors and time checks
  if (depts > 1) {
    demons <- c("int pop_nv = S + E + I + A + R; \n ",
                paste0("int pop_", 1:depts, " = S", 1:depts, " + E",
                       1:depts, " + I", 1:depts, " + A", 1:depts,
                       " + R", 1:depts, "; \n "),
                "int pop = pop_nv", paste0(" + pop_", 1:depts),
                "; \n int births = rpois(mu*pop*dt); \n ") %>%
      paste(collapse = "")
    tchecks <- paste0("if (vac_tcheck == ", 1:depts, ") { \n eta", 1:depts,
                      " = num_vacc/pop_nv/dt; \n } \n ") %>%
      paste(collapse = "")
  } else {
    demons <- c("int pop = S + E + I + A + R; \n ",
                "int births = rpois(mu*pop*dt); \n ") %>%
      paste(collapse = "")
  }

  # seasonal beta and foi
  time_check <- c(
    "double mybeta = beta1*seas1 + beta2*seas2 + beta3*seas3 + beta4*seas4 + beta5*seas5 + beta6*seas6; \n "
  )
  beta <- paste0(time_check, collapse = "")
  if (depts > 1) {
    foi_i <- c("I", paste0("+I", 1:depts)) %>%
      paste(collapse = "")
    foi_a <- c("A", paste0("+A", 1:depts)) %>%
      paste(collapse = "")
    foi <- paste0("double foi = pow(", foi_i, "+(1-kappa)*(", foi_a, "), nu)*mybeta/pop; \n ") %>%
      paste(collapse = "")
  } else {
    foi <- "double foi = pow(I, nu) * mybeta / pop; \n "
  }

  foi <- paste0(foi,
                "double sig_sq; \n ",
                "if (t < 233) { \n ",
                "sig_sq = sig_sq_epi; \n ",
                "} else { \n ",
                "sig_sq = sig_sq_end; \n ",
                "} \n ",
                "\n dgamma = rgammawn(sig_sq, dt); \n foi = foi * dgamma/dt; \n ")

  # theta_k
  if (vacscen != "id0") {
    thetas <- paste0("double theta", 1:depts, " = ve_d", 1:depts, "; \n ")
  }

  # compute rates
  s_rates <- c("Srate[0] = foi; \n ", "Srate[1] = delta; \n ")
  e_rates <- c("Erate[0] = sigma*(1-theta0); \n ", "Erate[1] = sigma*theta0; \n ", "Erate[2] = delta; \n ")
  i_rates <- c("Irate[0] = gamma; \n ", "Irate[1] = delta; \n ")
  a_rates <- c("Arate[0] = gamma; \n ", "Arate[1] = delta; \n ")
  r_rates <- c("Rrate[0] = alpha; \n ", "Rrate[1] = delta; \n ")
  if (depts > 1) {
    s_rates <- c(s_rates, paste0("Srate[", 2:(depts + 1), "] = eta", 1:depts, "; \n "),
                 paste0("S", 1:depts, "rate[0] = foi; \n ",
                        "S", 1:depts, "rate[1] = delta; \n "))
    e_rates <- c(e_rates, paste0("Erate[", 3:(depts + 2), "] = eta", 1:depts, "; \n "),
                 paste0("E", 1:depts, "rate[0] = sigma*(1-theta", 1:depts, "); \n ",
                        "E", 1:depts, "rate[1] = sigma*theta", 1:depts, "; \n ",
                        "E", 1:depts, "rate[2] = delta; \n "))
    i_rates <- c(i_rates, paste0("Irate[", 2:(depts + 1), "] = eta", 1:depts, "; \n "),
                 paste0("I", 1:depts, "rate[0] = gamma; \n ",
                        "I", 1:depts, "rate[1] = delta; \n "))
    a_rates <- c(a_rates, paste0("Arate[", 2:(depts + 1), "] = eta", 1:depts, "; \n "),
                 paste0("A", 1:depts, "rate[0] = gamma; \n ",
                        "A", 1:depts, "rate[1] = delta; \n "))
    r_rates <- c(r_rates, paste0("Rrate[", 2:(depts + 1), "] = eta", 1:depts, "; \n "),
                 paste0("R", 1:depts, "rate[0] = alpha; \n ",
                        "R", 1:depts, "rate[1] = delta; \n "))
  }
  rates <- c(s_rates, e_rates, i_rates, a_rates, r_rates) %>%
    paste(collapse = "")

  # transition numbers
  numbers <- paste0("reulermultinom(", depts + 2, ",S,&Srate[0],dt,&Strans[0]); \n ",
                    "reulermultinom(", depts + 3, ",E,&Erate[0],dt,&Etrans[0]); \n ",
                    "reulermultinom(", depts + 2, ",I,&Irate[0],dt,&Itrans[0]); \n ",
                    "reulermultinom(", depts + 2, ",A,&Arate[0],dt,&Atrans[0]); \n ",
                    "reulermultinom(", depts + 2, ",R,&Rrate[0],dt,&Rtrans[0]); \n ")
  if (depts > 1) {
    numbers <- c(numbers,
                 paste0("reulermultinom(2,S", 1:depts, ",&S", 1:depts, "rate[0],dt,&S", 1:depts, "trans[0]); \n "),
                 paste0("reulermultinom(3,E", 1:depts, ",&E", 1:depts, "rate[0],dt,&E", 1:depts, "trans[0]); \n "),
                 paste0("reulermultinom(2,I", 1:depts, ",&I", 1:depts, "rate[0],dt,&I", 1:depts, "trans[0]); \n "),
                 paste0("reulermultinom(2,A", 1:depts, ",&A", 1:depts, "rate[0],dt,&A", 1:depts, "trans[0]); \n "),
                 paste0("reulermultinom(2,R", 1:depts, ",&R", 1:depts, "rate[0],dt,&R", 1:depts, "trans[0]); \n ")) %>%
      paste(collapse = "")
  }

  # cohorts
  if (depts > 1) {
    coh_unvac <- c("S += -Strans[0]", paste0(" - Strans[", 1:(depts + 1), "]"), " + Rtrans[0] + births; \n ",
                   "E += -Etrans[0]", paste0(" - Etrans[", 1:(depts + 2), "]"), " + Strans[0]; \n ",
                   "I += -Itrans[0]", paste0(" - Itrans[", 1:(depts + 1), "]"), " + Etrans[0]; \n ",
                   "A += -Atrans[0]", paste0(" - Atrans[", 1:(depts + 1), "]"), " + Etrans[1]; \n ",
                   "R += -Rtrans[0]", paste0(" - Rtrans[", 1:(depts + 1), "]"), " + Itrans[0] + Atrans[0]; \n ")
    coh_vacs <- c(paste0("S", 1:depts, " += -S", 1:depts, "trans[0] - S", 1:depts,
                         "trans[1] + R", 1:depts, "trans[0] + Strans[", 2:(depts + 1), "]; \n "),
                  paste0("E", 1:depts, " += -E", 1:depts, "trans[0] - E", 1:depts, "trans[1] - E",
                         1:depts, "trans[2] + S", 1:depts, "trans[0] + Etrans[", 3:(depts + 2), "]; \n "),
                  paste0("I", 1:depts, " += -I", 1:depts, "trans[0] - I", 1:depts,
                         "trans[1] + E", 1:depts, "trans[0] + Itrans[", 2:(depts + 1), "]; \n "),
                  paste0("A", 1:depts, " += -A", 1:depts, "trans[0] - A", 1:depts,
                         "trans[1] + E", 1:depts, "trans[0] + Atrans[", 2:(depts + 1), "]; \n "),
                  paste0("R", 1:depts, " += -R", 1:depts, "trans[0] - R", 1:depts,
                         "trans[1] + A", 1:depts, "trans[0] + I", 1:depts, "trans[0] + Rtrans[", 2:(depts + 1), "]; \n "))
    coh <- c(coh_unvac, coh_vacs) %>%
      paste(collapse = "")
  } else {
    coh <- c("S += -Strans[0]", paste0(" - Strans[", 1:(depts + 1), "]"), " + Rtrans[0] + births; \n ",
             "E += -Etrans[0]", paste0(" - Etrans[", 1:(depts + 2), "]"), " + Strans[0]; \n ",
             "I += -Itrans[0]", paste0(" - Itrans[", 1:(depts + 1), "]"), " + Etrans[0]; \n ",
             "A += -Atrans[0]", paste0(" - Atrans[", 1:(depts + 1), "]"), " + Etrans[1]; \n ",
             "R += -Rtrans[0]", paste0(" - Rtrans[", 1:(depts + 1), "]"), " + Itrans[0] + Atrans[0]; \n ") %>%
      paste(collapse = "")
  }

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

  pars <- unlist(MODEL1_INPUT_PARAMETERS$joint_pars)

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
