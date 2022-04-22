#' Fit Model 1
#'
#' This function fits the specified version of model 1.
#'
#' As of Feb 20, 2022, Model 1 can be written in 4 different forms:
#' \itemize{
#'    \item Basic: aggregated, single time period `pomp` object (`haiti1()`).
#'    \item Joint: aggregated, joint time periods `pomp` object (`haiti1_joint()`).
#'    \item Disagg: disaggregated `pomp` object (`haiti1_disagg()`).
#'    \item Panel: disaggregated `panelPomp` object (`haiti1_disagg()`).
#' }
#'
#' The Basic version is the closest to the model described and used by Lee et al.
#' in that it only has an additional overdispersion term. Basic-Full is the same
#' as Basic but fitting through all 430 weeks of data (not split between endemic
#' and epidemic). We currently provide methods for fitting both the Basic and Joint forms.
#'
#' @param version In c("basic", "basic-full", "joint").
#' @param run_level in c(1, 2, 3)
#'
#' @examples
#' m1 <- haiti1_joint()
#' fit_m1 <- fit_haiti1(version = "joint", run_level = 2)
#'
#' @importFrom magrittr %>%
#' @importFrom foreach %do%
#' @importFrom foreach foreach
#'
#' @export
fit_haiti1 <- function(version = "basic", run_level = 1) {
  # Set Run Level
  num_parts <- switch(run_level, 1e3, 1e3, 5e3)  # Number of particles
  num_iters <- switch(run_level, 10, 100, 200)  # Number of MIF iterations
  num_profs <- switch(run_level,  1,  3,  5)  # Number of profiles

  ## Fixed parameters
  gamma <- 7/2 ## 2 day infectious period
  sigma <- 7/1.4 ## 1.4 day latent period
  theta0 <- 0 ## 0% going to asymptomatic in base model
  alpha <- 7/2920 ## 8 year mean duration of natural immunity
  ## 22.6/1000 average annual birth rate, adjusted for compounding by week
  mu <- ((1+22.6/1000)^(1/52.14))-1
  ## 7.5/1000 average annual death rate, adjusted for compounding by week
  delta <- ((1+7.5/1000)^(1/52.14))-1

  ## beta parameter settings
  blo <- 1E-9; bup <- 10 ## uniform beta settings
  ## median R0 among init values: median(exp(rnorm(1000, log(bmn), bse)))*2/7
  bmn <- 4.5; bse <- 0.5
  ## nu parameter settings
  nlo <- 0.95; nup <- 1
  ## sigma_se parameter settings
  sigsqlo <- 1E-9; sigsqup <- 5

  if (version == 'basic' || version == 'basic-full') {
    ## make models
    if (version == 'basic') { ## separate periods
      epi_mod <- haiti1(period = "epidemic")
      end_mod <- haiti1(period = "endemic")
      pop_haiti_epi <- MODEL1_INPUT_PARAMETERS$adj_pars_epi[22] %>% unlist()
      pop_haiti_end <- MODEL1_INPUT_PARAMETERS$adj_pars_end[22] %>% unlist()
    } else { ## basic-full
      full_mod <- haiti1(period = "1-430")
      pop_haiti_epi <- MODEL1_INPUT_PARAMETERS$adj_pars_epi[22] %>% unlist()
    }

    haiti.dat <- haiti1_agg_data() ## get case data (country-wide)

    ## starting states
    E0 <- 10/pop_haiti_epi ## rpois(nsamps, 10)/pop
    I0 <- 10/pop_haiti_epi ## rpois(nsamps, 10)/pop
    A0 <- 0.0
    R0 <- 0.0
    S0 <- 1-R0-I0-E0-A0

    ## make guesses
    rhosamps <- pomp::profile_design(rho = seq(1E-8, 1, length = 20),
                                     upper = c(tau = 20, beta1 = bup, beta2 = bup,
                                               beta3 = bup, beta4 = bup, beta5 = bup,
                                               beta6 = bup, nu = nup, sig_sq = sigsqup),
                                     lower = c(tau = 1, beta1 = blo, beta2 = blo,
                                               beta3 = blo, beta4 = blo, beta5 = blo,
                                               beta6 = blo, nu = nlo, sig_sq = sigsqlo),
                                     nprof = num_profs)
    tausamps <- pomp::profile_design(tau = seq(1, 20, length = 20),
                                     upper = c(rho = 1, beta1 = bup, beta2 = bup,
                                               beta3 = bup, beta4 = bup, beta5 = bup,
                                               beta6 = bup, nu = nup, sig_sq = sigsqup),
                                     lower = c(rho = 1E-8, beta1 = blo, beta2 = blo,
                                               beta3 = blo, beta4 = blo, beta5 = blo,
                                               beta6 = blo, nu = nlo, sig_sq = sigsqlo),
                                     nprof = num_profs)
    betasamps <- pomp::profile_design(beta1 = seq(blo, bup, length = 20),
                                      upper = c(rho = 1, tau = 20, beta2 = bup,
                                                beta3 = bup, beta4 = bup, beta5 = bup,
                                                beta6=bup, nu = nup, sig_sq = sigsqup),
                                      lower = c(rho = 1E-8, tau = 1, beta2 = blo,
                                                beta3=blo, beta4 = blo, beta5 = blo,
                                                beta6=blo, nu = nlo, sig_sq = sigsqlo),
                                      nprof = num_profs)
    nusamps <- pomp::profile_design(nu = seq(nlo, nup, length = 20),
                                    upper = c(rho = 1, tau = 20, beta1 = bup, beta2 = bup,
                                              beta3 = bup, beta4 = bup, beta5 = bup,
                                              beta6 = bup, sig_sq = sigsqup),
                                    lower = c(rho = 1E-8, tau = 1, beta1 = blo,
                                              beta2 = blo, beta3 = blo, beta4 = blo,
                                              beta5 = blo, beta6 = blo, sig_sq = sigsqlo),
                                    nprof = num_profs)
    sigsqsamps <- pomp::profile_design(sig_sq = seq(sigsqlo, sigsqup, length = 20),
                                       upper = c(rho = 1, tau = 20, beta1 = bup, beta2 = bup,
                                                 beta3 = bup, beta4 = bup, beta5 = bup,
                                                 beta6 = bup),
                                       lower = c(rho = 1E-8, tau = 1, beta1 = blo,
                                                 beta2 = blo, beta3 = blo, beta4 = blo,
                                                 beta5 = blo, beta6 = blo),
                                       nprof = num_profs)
    starts <- dplyr::bind_rows(rhosamps, tausamps, betasamps, nusamps, sigsqsamps) %>%
      dplyr::mutate(parid = seq_along(rho)) %>%
      dplyr::mutate(theta0 = theta0,mu = mu, delta = delta, nu = nu, sigma = sigma, alpha = alpha,
             gamma = gamma, S_0 = S0, E_0 = E0, I_0 = I0, A_0 = A0, R_0 = R0, pop_0 = pop_haiti_epi) %>%
      dplyr::select(parid, rho, tau, beta1, beta2, beta3, beta4, beta5,
             beta6, nu, gamma, sigma, theta0, alpha, mu, delta, sig_sq,
             S_0, E_0, I_0, A_0, R_0, pop_0)

    est_params_epi <- c(paste0("beta", 1:6), "rho", "tau", "nu", "sig_sq", "E_0", "I_0")
    rw_sds_epi <- pomp::rw.sd(beta1 = 0.02, beta2 = 0.02, beta3 = 0.02,
                              beta4 = 0.02, beta5 = 0.02, beta6 = 0.02,
                              tau = 0.02, rho = 0.02, nu = 0.02,
                              sig_sq = 0.02, E_0 = ivp(0.2), I_0 = ivp(0.2))

    if (version == 'basic') { ## separate periods
      est_params_end <- c(paste0("beta", 1:6), "rho", "tau", "nu", "sig_sq")
      rw_sds_end <- pomp::rw.sd(beta1 = 0.02, beta2 = 0.02, beta3 = 0.02,
                                beta4 = 0.02, beta5 = 0.02, beta6 = 0.02,
                                tau = 0.02, rho = 0.02, nu = 0.02, sig_sq = 0.02)
      full_mod <- epi_mod
    }

    foreach(start = iterators::iter(starts, by = 'row'),
           .combine = rbind, .inorder = FALSE,
           .packages = c("pomp", "magrittr"),
           .errorhandling = c('remove'),
           .export = c("haiti.dat", "num_betas", "pop.haiti",
                       "covartab", "settings"),
           .noexport = c(),
           .verbose = TRUE) %dopar%
    {
         po <- full_mod
         allpars <- c("rho","tau","beta1","beta2","beta3","beta4","beta5","beta6",
                      "gamma","sigma","theta0","alpha","mu","delta","nu", "sig_sq",
                      "S_0","E_0","I_0","A_0","R_0", "pop_0")
         start_coef <- unlist(start[which(names(start) %in% allpars)])
         names(start_coef) <- c("rho","tau","beta1","beta2","beta3","beta4","beta5","beta6",
                                "nu", "gamma", "sigma", "theta0", "alpha", "mu", "delta", "sig_sq",
                                "S_0", "E_0", "I_0", "A_0", "R_0", "pop_0")
         coef(po) <- start_coef

         ## if no Exposed/Infectious, make it small but nonzero
         if (coef(po, "E_0") == 0.0 ) {
           coef(po, "E_0") <- 1e-9
         }
         if (coef(po, "I_0") == 0.0) {
           coef(po, "I_0") <- 1e-9
         }

         ## perform iterated filtering
         mf.mod_epi <- pomp::mif2(po, Nmif = num_iters,
                                  rw.sd = rw_sds_epi,
                                  Np = num_parts,
                                  cooling.type = "hyperbolic",
                                  cooling.fraction.50 = 0.5,
                                  verbose = FALSE)

         ## get likelihood estimate
         ll_epi <- est_logLik1(version = "basic", model = mf.mod_epi,
                               Np = num_parts, nreps = num_iters)

         if (version == 'basic') { ## separate periods
           sims <- pomp::simulate(po, params = po@params,
                                  nsim = 25, format = "data.frame") %>%
             dplyr::mutate(pop = S + E + I + A + R) %>%
             dplyr::select(week, S, E, I, A, R, incid, pop, cases)

           states <- sims %>%
             dplyr::group_by(week) %>%
             dplyr::summarise(S_med = median(S), E_med = median(E), I_med = median(I),
                              A_med = median(A), R_med = median(R),
                              incid_med = median(incid), pop_med = median(pop),
                              cases_med = median(cases),
                              S_mean = mean(S), E_mean = mean(E), I_mean = mean(I),
                              A_mean = mean(A), R_mean = mean(R),
                              incid_mean = mean(incid), pop_mean = mean(pop),
                              cases_mean = mean(cases),
                              S_lo = quantile(S, probs = c(.025)),
                              E_lo = quantile(E, probs = c(.025)),
                              I_lo = quantile(I, probs = c(.025)),
                              A_lo = quantile(A, probs = c(.025)),
                              R_lo = quantile(R, probs = c(.025)),
                              incid_lo = quantile(incid, probs = c(.025)),
                              pop_lo = quantile(pop, probs = c(.025)),
                              cases_lo = quantile(cases, probs = c(.025)),
                              S_hi = quantile(S, probs = c(.975)),
                              E_hi = quantile(E, probs = c(.975)),
                              I_hi = quantile(I, probs = c(.975)),
                              A_hi = quantile(A, probs = c(.975)),
                              R_hi = quantile(R, probs = c(.975)),
                              incid_hi = quantile(incid, probs = c(.975)),
                              pop_hi = quantile(pop, probs = c(.975)),
                              cases_hi = quantile(cases, probs = c(.975))) %>%
             dplyr::mutate(date = lubridate::ymd("2010-10-16") + lubridate::weeks(week))

           end_starts <- states %>%
             dplyr::filter(week == max(week)) %>%
             dplyr::mutate(S_0 = S_med/pop_med, E_0 = E_med/pop_med,
                           I_0 = I_med/pop_med, A_0 = A_med/pop_med,
                           incid_0 = incid_med) %>%
             dplyr::mutate(R_0 = 1-(S_0+E_0+I_0+A_0), pop_0 = pop_med) %>%
             dplyr::select(S_0, E_0, I_0, A_0, R_0, incid_0, pop_0)

           epi_coef <- coef(po) %>%
             t() %>%
             data.frame()
           epi_coef <- epi_coef %>%
             dplyr::select("rho","tau","beta1","beta2","beta3","beta4","beta5","beta6",
                           "gamma","sigma","theta0","alpha","mu","delta","nu", "sig_sq")
           end_coef <- c(epi_coef, end_starts)

           po <- end_mod
           timezero(po) <- po@times[1] - 1
           coef(po) <- unlist(end_coef)

           ## if no Exposed/Infectious, make it small but nonzero
           if (coef(po, "E_0") == 0.0 ) {
             coef(po, "E_0") <- 1e-9
           }
           if (coef(po, "I_0") == 0.0) {
             coef(po, "I_0") <- 1e-9
           }

           ## perform iterated filtering
           mf.mod_end <- pomp::mif2(po, Nmif = num_iters,
                                    rw.sd = rw_sds_end,
                                    Np = num_parts,
                                    cooling.type = "hyperbolic",
                                    cooling.fraction.50 = 0.5,
                                    verbose = FALSE)

           ## get likelihood estimate
           ll_end <- est_logLik1(version = "basic", model = mf.mod_end,
                                 Np = num_parts, nreps = num_iters)

           ## record parameter estimates
           dummy <- data.frame(as.list(coef(mf.mod_epi)),
                               as.list(coef(mf.mod_end)),
                               loglik_epi = ll_epi[1],
                               loglik.se_epi = ll_epi[2],
                               loglik_end = ll_end[1],
                               loglik.se_end = ll_end[2])

           rm(mf.mod_epi, mf.mod_end, ll_epi, ll_end, sims)
         } else {
           dummy <- data.frame(as.list(coef(mf.mod_epi)),
                               loglik_full = ll_epi[1],
                               loglik.se_full = ll_epi[2])
           rm(mf.mod_epi, ll_epi)
         }
         gc()
         dummy
    } ->  fits

    if (version == 'basic') { ## separate periods
      epi_best <- fits %>%
        dplyr::arrange(-loglik_epi) %>%
        dplyr::slice_head(n = 1) %>%
        dplyr::pull(which)

      end_best <- fits %>%
        dplyr::arrange(-loglik_end) %>%
        dplyr::slice_head(n = 1) %>%
        dplyr::pull(which)

      res <- list("epidemic_fit" = epi_best,
                  "endemic_fit" = end_best)
    } else {
      best <- fits %>%
        dplyr::arrange(-loglik_full) %>%
        dplyr::slice_head(n = 1) %>%
        dplyr::pull(which)

      res <- list("full_fit" = best)
    }
  } else if (version == 'joint') {
    ## make model
    mod <- haiti1_joint(rho_flag = T, tau_flag = T, sig_sq_flag = T, beta_flag = F, nu_flag = F)
    haiti.dat <- haiti1_agg_data() ## get case data (country-wide)
    pop_haiti <- MODEL1_INPUT_PARAMETERS$adj_pars_epi[22] %>% unlist()

    ## starting states
    E0 <- 10/pop_haiti ## rpois(nsamps, 10)/pop
    I0 <- 10/pop_haiti ## rpois(nsamps, 10)/pop
    A0 <- 0.0
    R0 <- 0.0
    S0 <- 1-R0-I0-E0-A0

    rhosamps <- pomp::profile_design(rho = seq(1E-8, 1, length = 30),
                               upper = c(tau = 20, beta1 = bup, beta2 = bup,
                                         beta3 = bup, beta4 = bup, beta5 = bup,
                                         beta6 = bup, nu = nup, sig_sq = sigsqup),
                               lower = c(tau = 1, beta1 = blo, beta2 = blo,
                                         beta3 = blo, beta4 = blo, beta5 = blo,
                                         beta6 = blo, nu = nlo, sig_sq = sigsqlo),
                               nprof = num_profs)
    tausamps <- pomp::profile_design(tau = seq(1, 20, length = 30),
                               upper = c(rho = 1, beta1 = bup, beta2 = bup,
                                         beta3 = bup, beta4 = bup, beta5 = bup,
                                         beta6 = bup, nu = nup, sig_sq = sigsqup),
                               lower = c(rho = 1E-8, beta1 = blo, beta2 = blo,
                                         beta3 = blo, beta4 = blo, beta5 = blo,
                                         beta6 = blo, nu = nlo, sig_sq = sigsqlo),
                               nprof = num_profs)
    betasamps <- pomp::profile_design(beta1 = seq(blo, bup, length = 30),
                                upper = c(rho = 1, tau = 20, beta2 = bup,
                                          beta3 = bup, beta4 = bup, beta5 = bup,
                                          beta6=bup, nu = nup, sig_sq = sigsqup),
                                lower = c(rho = 1E-8, tau = 1, beta2 = blo,
                                          beta3=blo, beta4 = blo, beta5 = blo,
                                          beta6=blo, nu = nlo, sig_sq = sigsqlo),
                                nprof = num_profs)
    nusamps <- pomp::profile_design(nu = seq(nlo, nup, length = 30),
                              upper = c(rho = 1, tau = 20, beta1 = bup, beta2 = bup,
                                        beta3 = bup, beta4 = bup, beta5 = bup,
                                        beta6 = bup, sig_sq = sigsqup),
                              lower = c(rho = 1E-8, tau = 1, beta1 = blo,
                                        beta2 = blo, beta3 = blo, beta4 = blo,
                                        beta5 = blo, beta6 = blo, sig_sq = sigsqlo),
                              nprof = num_profs)
    sigsqsamps <- pomp::profile_design(sig_sq = seq(sigsqlo, sigsqup, length = 30),
                                       upper = c(rho = 1, tau = 20, beta1 = bup, beta2 = bup,
                                                 beta3 = bup, beta4 = bup, beta5 = bup,
                                                 beta6 = bup, nu = nup),
                                       lower = c(rho = 1E-8, tau = 1, beta1 = blo,
                                                 beta2 = blo, beta3 = blo, beta4 = blo,
                                                 beta5 = blo, beta6 = blo, nu = nlo),
                                       nprof = num_profs)
    starts <- dplyr::bind_rows(rhosamps, tausamps, betasamps, nusamps, sigsqsamps) %>%
      dplyr::mutate(parid = seq_along(rho)) %>%
      dplyr::mutate(theta0 = theta0, mu = mu, delta = delta, nu = nu,
                    sigma = sigma, alpha = alpha, kappa = 0, gamma = gamma,
                    S_0 = S0, E_0 = E0, I_0 = I0, A_0 = A0, R_0 = R0, pop_0 = pop_haiti) %>%
      dplyr::select(parid, rho, tau, beta1, beta2, beta3, beta4, beta5,
                    beta6, nu, gamma, sigma, theta0, alpha, mu, delta, sig_sq,
                    kappa, S_0, E_0, I_0, A_0, R_0, pop_0)
    add_starts <- starts %>%
      dplyr::select(rho, tau, sig_sq) %>%
      dplyr::rename(rho_end = rho, sig_sq_end = sig_sq, tau_end = tau)
    starts <- cbind(starts, add_starts) %>%
      dplyr::rename(rho_epi = rho, tau_epi = tau, sig_sq_epi = sig_sq)

    est_params <- c(paste0("beta", 1:6),
                    "rho_epi", "rho_end", "tau_epi",
                    "tau_end", "nu", "sig_sq_epi", "sig_sq_end", "E_0", "I_0")
    rw_sds_joint <- pomp::rw.sd(beta1 = 0.02, beta2 = 0.02, beta3 = 0.02,
                                beta4 = 0.02, beta5 = 0.02, beta6 = 0.02,
                                rho_epi = ifelse(time > 232, 0.0, 0.02),
                                rho_end = ifelse(time <= 232, 0.0, 0.02),
                                tau_epi = ifelse(time > 232, 0.0, 0.02),
                                tau_end = ifelse(time <= 232, 0.0, 0.02),
                                sig_sq_epi = ifelse(time > 232, 0.0, 0.02),
                                sig_sq_end = ifelse(time <= 232, 0.0, 0.02),
                                nu = 0.02,
                                E_0 = ifelse(time > 232, 0.0, pomp::ivp(0.2)),
                                I_0 = ifelse(time > 232, 0.0, pomp::ivp(0.2)))

    foreach(start = iter(starts, by = 'row'),
            .combine = rbind, .inorder = FALSE,
            .packages = c("pomp", "magrittr"),
            .errorhandling = c('remove'),
            .export = c("haiti.dat", "num_betas", "pop.haiti",
                        "covartab", "settings"),
            .noexport = c(),
            .verbose = TRUE) %dopar%
      {
        po <- mod
        timezero(po) <- 0
        allpars <- c("rho_epi", "rho_end",
                     "beta1", "beta2", "beta3",
                     "beta4", "beta5", "beta6",
                     "sig_sq_epi", "sig_sq_end", "tau_epi", "tau_end",
                     "gamma","sigma","theta0","alpha","mu","delta","nu", "kappa",
                     "S_0","E_0","I_0","A_0","R_0", "pop_0")
        coef(po) <- unlist(start[which(names(start) %in% allpars)])

        ## if no Exposed/Infectious, make it small but nonzero
        if (coef(po, "E_0") == 0.0 ) {
          coef(po, "E_0") <- 1e-9
        }
        if (coef(po, "I_0") == 0.0) {
          coef(po, "I_0") <- 1e-9
        }

        ## perform iterated filtering
        mf.mod <- pomp::mif2(po, Nmif = num_iters,
                             rw.sd = rw_sds_joint,
                             Np = num_parts,
                             cooling.type = "hyperbolic",
                             cooling.fraction.50 = 0.5,
                             verbose = FALSE)

        ll <- est_logLik1(version = "joint", model = mf.mod,
                          Np = num_parts, nreps = num_iters)

        ## record parameter estimates
        dummy <- data.frame(as.list(coef(mf.mod)),
                            loglik = ll[1],
                            loglik.se = ll[2])

        rm(mf.mod, ll)
        gc()
        dummy
      } -> fits

    res <- fits %>%
      dplyr::arrange(-loglik) %>%
      dplyr::slice_head(n = 1) %>%
      dplyr::pull(which)

  } else {
    stop(paste0("Version = \"", version, "\" is not a valid input.
                Version must be in {\"basic\", \"joint\"}."))
  }
  return(res)
}
