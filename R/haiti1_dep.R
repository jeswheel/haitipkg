#' Build pomp object for a single department in Model 1
#'
#' Generate a class \sQuote{pomp} object for fitting to epidemic/endemic Haiti cholera data.
#' Replace " " in departement names with "_".
#'
#' @param departement String for which departement to use for the pomp object
#' @param vacscen String for what vaccination campaign to do
#' @importFrom pomp Csnippet
#' @return An object of class \sQuote{pomp}.
#' @examples
#' m1 <- haiti1_dep('Artibonite', 'id0')
#' @export

haiti1_dep <- function(dept = 'Artibonite', vacscen = 'id0') {
    ## get data
    dat <- haiti1_data() %>%
      dplyr::select(-1) %>%
      dplyr::rename(Grand_Anse = Grand.Anse,
                    Nord_Est = Nord.Est,
                    Nord_Ouest = Nord.Ouest,
                    Sud_Est = Sud.Est)
    dat <- dat %>%
      dplyr::select(c(dept, "week"))
    colnames(dat) <- c("cases", "week")
    other_cases <- haiti1_data()  %>%
      dplyr::select(-1) %>%
      dplyr::rename(Grand_Anse = Grand.Anse,
                    Nord_Est = Nord.Est,
                    Nord_Ouest = Nord.Ouest,
                    Sud_Est = Sud.Est) %>%
      dplyr::select(-dept)
    other_cases$other_cases = rowSums(other_cases[ , c(1,9)])
    other_cases <- other_cases %>%
      dplyr::select(c("week", "other_cases"))
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
    covar <- dplyr::left_join(covar, other_cases, by = c("time" = "week"))

    ## determine if dept is included in vaccination campaign
    vac <- FALSE
    if (vacscen == 'id0') { ## no vaccination
      vac <- FALSE
    } else if (vacscen == 'id1') { ## fast national
      vac <- TRUE
    } else if (vacscen == 'id2') { ## 2-dept
      if (dept == 'Centre' | dept == 'Artibonite') {
        vac <- TRUE
      }
    } else if (vacscen == 'id3') { ## slow national
      vac <- TRUE
    } else if (vacscen == 'id4') { ## 3-dept
      if (dept == 'Centre' | dept == 'Artibonite' | dept == 'Ouest') {
        vac <- TRUE
      }
    } else { ## fast national high-coverage
      vac <- TRUE
    }

    dept_order <- switch(dept,
                         "Centre" = 1,
                         "Artibonite" = 2,
                         "Ouest" = 3,
                         "Nord_Ouest" = 4,
                         "Nord" = 5,
                         "Sud" = 6,
                         "Nippes" = 7,
                         "Nord_Est" = 8,
                         "Sud_Est" = 9,
                         "Grand_Anse" = 10)

    ## make components pomp object building
    ## rinit
    state_names_base <- c("S", "E", "I", "A", "R")
    state_names <- state_names_base
    if (vac) {
      state_names <- c(state_names, paste0(state_names_base, 1))
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
    if (vac) {
      rinit_paste <- c(rinit_paste, "incidU = 0.0; \n ", "incidV = 0.0; \n",
                       "asymV = 0.0; \n ", "newV = 0.0; \n ")
    }

    rinit <- Csnippet(
      paste(rinit_paste, collapse = "")
    )

    ## rprocess
    # transition rates and numbers
    if (vac) {
      rates <- rep(c(2, 3, 2, 2, 2), 2)
    } else {
      rates <- c(2, 3, 2, 2, 2)
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
    if (vac) {
      vac_rates <- paste0("double eta = 0.0; \n ")
    }

    # demonitors and time checks
    if (vac) {
      demons <- c("int pop_nv = S + E + I + A + R; \n ",
                  "int pop_1 = S1 + E1 + I1 + A1 + R1; \n",
                  "int pop = pop_nv + pop_1",
                  "; \n int births = rpois(mu*pop*dt); \n ") %>%
        paste(collapse = "")
      tchecks <- c(paste0("if (vac_tcheck == ", dept_order, ") { \n eta",
                   " = num_vacc/pop_nv/dt; \n } \n ")) %>%
        paste(collapse = "")
    } else {
      demons <- c("int pop = S + E + I + A + R; \n ",
                  "int births = rpois(mu*pop*dt); \n ") %>%
        paste(collapse = "")
    }

    # seasonal beta and foi
    beta <- "double mybeta = beta1*seas1 + beta2*seas2 + beta3*seas3 + beta4*seas4 + beta5*seas5 + beta6*seas6; \n "
    mob <- "double mobility;
    if(ISNA(other_cases)) {
      mobility = 0.0;
    } else {
      mobility = mob_c * other_cases;
    } \n
    "
    if (vac) {
      foi_i <- c("I + I1") %>%
        paste(collapse = "")
      foi_a <- c("A + A1") %>%
        paste(collapse = "")
      foi <- paste0("double foi = pow(", foi_i, "+(1-kappa)*(", foi_a, "), nu)*mybeta/pop + mobility; \n ") %>%
        paste(collapse = "")
    } else {
      foi <- "double foi = pow(I, nu) * mybeta / pop + mobility; \n "
    }

    foi <- paste0(mob,
                  foi,
                  "double sig_sq; \n ",
                  "if (t < 233) { \n ",
                  "sig_sq = sig_sq_epi; \n ",
                  "} else { \n ",
                  "sig_sq = sig_sq_end; \n ",
                  "} \n ",
                  "\n dgamma = rgammawn(sig_sq, dt); \n foi = foi * dgamma/dt; \n ")

    # theta_k
    if (vac) {
      thetas <- paste0("double theta1 = ve_d", dept_order, "; \n ")
    }

    # compute rates
    s_rates <- c("Srate[0] = foi; \n ", "Srate[1] = delta; \n ")
    e_rates <- c("Erate[0] = sigma*(1-theta0); \n ", "Erate[1] = sigma*theta0; \n ", "Erate[2] = delta; \n ")
    i_rates <- c("Irate[0] = gamma; \n ", "Irate[1] = delta; \n ")
    a_rates <- c("Arate[0] = gamma; \n ", "Arate[1] = delta; \n ")
    r_rates <- c("Rrate[0] = alpha; \n ", "Rrate[1] = delta; \n ")
    if (vac) {
      s_rates <- c(s_rates, paste0("Srate[2] = eta; \n "),
                   paste0("S1rate[0] = foi; \n ",
                          "S1rate[1] = delta; \n "))
      e_rates <- c(e_rates, paste0("Erate[3] = eta; \n "),
                   paste0("E1rate[0] = sigma*(1-theta1); \n ",
                          "E1rate[1] = sigma*theta1; \n ",
                          "E1rate[2] = delta; \n "))
      i_rates <- c(i_rates, paste0("Irate[2] = eta; \n "),
                   paste0("I1rate[0] = gamma; \n ",
                          "I1rate[1] = delta; \n "))
      a_rates <- c(a_rates, paste0("Arate[2] = eta; \n "),
                   paste0("A1rate[0] = gamma; \n ",
                          "A1rate[1] = delta; \n "))
      r_rates <- c(r_rates, paste0("Rrate[2] = eta; \n "),
                   paste0("R1rate[0] = alpha; \n ",
                          "R1rate[1] = delta; \n "))
    }
    rates <- c(s_rates, e_rates, i_rates, a_rates, r_rates) %>%
      paste(collapse = "")

    # transition numbers
    if (vac) {
      numbers <- paste0("reulermultinom(3,S,&Srate[0],dt,&Strans[0]); \n ",
                        "reulermultinom(4,E,&Erate[0],dt,&Etrans[0]); \n ",
                        "reulermultinom(3,I,&Irate[0],dt,&Itrans[0]); \n ",
                        "reulermultinom(3,A,&Arate[0],dt,&Atrans[0]); \n ",
                        "reulermultinom(3,R,&Rrate[0],dt,&Rtrans[0]); \n ")
      numbers <- c(numbers,
                   paste0("reulermultinom(2,S1,&S1rate[0],dt,&S1trans[0]); \n "),
                   paste0("reulermultinom(3,E1,&E1rate[0],dt,&E1trans[0]); \n "),
                   paste0("reulermultinom(2,I1,&I1rate[0],dt,&I1trans[0]); \n "),
                   paste0("reulermultinom(2,A1,&A1rate[0],dt,&A1trans[0]); \n "),
                   paste0("reulermultinom(2,R1,&R1rate[0],dt,&R1trans[0]); \n ")) %>%
        paste(collapse = "")
    } else {
      numbers <- paste0("reulermultinom(2,S,&Srate[0],dt,&Strans[0]); \n ",
                        "reulermultinom(3,E,&Erate[0],dt,&Etrans[0]); \n ",
                        "reulermultinom(2,I,&Irate[0],dt,&Itrans[0]); \n ",
                        "reulermultinom(2,A,&Arate[0],dt,&Atrans[0]); \n ",
                        "reulermultinom(2,R,&Rrate[0],dt,&Rtrans[0]); \n ")
    }

    # cohorts
    if (vac) {
      coh_unvac <- c("S += -Strans[0]", paste0(" - Strans[", 1:2, "]"), " + Rtrans[0] + births; \n ",
                     "E += -Etrans[0]", paste0(" - Etrans[", 1:3, "]"), " + Strans[0]; \n ",
                     "I += -Itrans[0]", paste0(" - Itrans[", 1:2, "]"), " + Etrans[0]; \n ",
                     "A += -Atrans[0]", paste0(" - Atrans[", 1:2, "]"), " + Etrans[1]; \n ",
                     "R += -Rtrans[0]", paste0(" - Rtrans[", 1:2, "]"), " + Itrans[0] + Atrans[0]; \n ")
      coh_vacs <- c(paste0("S1 += -S1trans[0] - S1trans[1] + R1trans[0] + Strans[2]; \n "),
                    paste0("E1 += -E1trans[0] - E1trans[1] - E1trans[2] + S1trans[0] + Etrans[3]; \n "),
                    paste0("I1 += -I1trans[0] - I1trans[1] + E1trans[0] + Itrans[2]; \n "),
                    paste0("A1 += -A1trans[0] - A1trans[1] + E1trans[0] + Atrans[2]; \n "),
                    paste0("R1 += -R1trans[0] - R1trans[1] + A1trans[0] + I1trans[0] + Rtrans[2]; \n "))
      coh <- c(coh_unvac, coh_vacs) %>%
        paste(collapse = "")
    } else {
      coh <- c("S += -Strans[0]", paste0(" - Strans[", 1:1, "]"), " + Rtrans[0] + births; \n ",
               "E += -Etrans[0]", paste0(" - Etrans[", 1:2, "]"), " + Strans[0]; \n ",
               "I += -Itrans[0]", paste0(" - Itrans[", 1:1, "]"), " + Etrans[0]; \n ",
               "A += -Atrans[0]", paste0(" - Atrans[", 1:1, "]"), " + Etrans[1]; \n ",
               "R += -Rtrans[0]", paste0(" - Rtrans[", 1:1, "]"), " + Itrans[0] + Atrans[0]; \n ") %>%
        paste(collapse = "")
    }

    # incidence
    incids <- c("incid += Etrans[0]; \n ")
    foi_val <- "foival += foi; \n "
    str0 <- "Str0 += Strans[0]; \n "
    sin <- c("Sin += Rtrans[0] + births; \n ")
    sout <- c("Sout += Strans[0] + Strans[1]; \n ")
    last <- c(foi_val, str0, sin)

    if (vac) {
      incids <- "incid += Etrans[0] + E1trans[0]; \n "
      incid_u <- "incidU += Etrans[0]; \n "
      incid_v <- "incidV += E1trans[0]; \n "
      asym_v <- "asymV += E1trans[1]; \n "
      new_v <- "newV += Strans[2] + Etrans[3] + Itrans[2] + Atrans[2] + Rtrans[2]; \n "
      sout <- c("Sout += Strans[0]", paste0(" + Strans[", 1:2, "]"), "; \n ")
      last <- c(last, incid_u, incid_v, asym_v, new_v)
    }
    last <- c(last, incids, sout) %>%
      paste(collapse = "")

    if (vac) {
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
    if (vac) {
      ## state names
      state_names <- c("S", "S1",
                       "E", "E1",
                       "I", "I1",
                       "A", "A1",
                       "R", "R1",
                       "incid", "incidU", "incidV", "asymV", "newV",
                       "foival", "Str0", "Sout", "Sin")

      ## parameter names
      param_names <- c("rho_epi", "sig_sq_epi", "tau_epi", # epidemic
                       "rho_end", "sig_sq_end", "tau_end", # endemic
                       "mob_c", # mobility term intensity
                       "beta1", "beta2", "beta3", "beta4", "beta5", "beta6",
                       "gamma", "sigma", "theta0", "alpha", "mu", "delta",
                       "nu", "kappa", "pop_0",
                       "S_0","E_0","I_0","A_0","R_0",
                       "S1_0","E1_0","I1_0","A1_0","R1_0")

      ## accum vars
      accum_names <- c("incid","incidU","incidV","asymV","newV",
                       "foival","Str0","Sout","Sin")
    } else {
      ## state names
      state_names <- c("S", "E", "I", "A", "R", "incid", "foival", "Str0", "Sout", "Sin")

      ## parameter names
      param_names <- c("rho_epi", "sig_sq_epi", "tau_epi", # epidemic
                       "rho_end", "sig_sq_end", "tau_end", # endemic
                       "mob_c", # mobility term intensity
                       "beta1", "beta2", "beta3", "beta4", "beta5", "beta6",
                       "gamma", "sigma", "theta0", "alpha", "mu", "delta",
                       "nu", "kappa", "pop_0",
                       "S_0", "E_0", "I_0", "A_0", "R_0")

      ## accum vars
      accum_names <- c("incid", "foival","Str0","Sout","Sin")
    }

    ## partrans
    param_trans <- pomp::parameter_trans(
      log = c("beta1", "beta2", "beta3", "beta4", "beta5", "beta6",
              "sigma", "gamma", "mu", "delta", "alpha",
              "sig_sq_end", "sig_sq_epi", "tau_end", "tau_epi",
              "mob_c"),
      logit = c("rho_epi", "rho_end", "nu", "theta0"),
      barycentric = c("S_0", "E_0", "I_0", "A_0", "R_0")
    )

    params <- MODEL1_INPUT_PARAMETERS$dep_params[[dept]] %>%
      unlist()

    if (vac) {
      par_names <- names(params)
      params <- c(params, rep(0.0, 5))
      names(params) <- c(par_names, "S1_0", "E1_0", "I1_0", "A1_0", "R1_0")
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
      params = params,
      statenames = state_names,
      paramnames = param_names,
      accumvars = accum_names,
      rinit = rinit
    )

    return(model1)
  }

