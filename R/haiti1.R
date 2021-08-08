#' Build pomp object for Model 1
#'
#' Generate a class \sQuote{pomp} object for fitting to epidemic/endemic Haiti cholera data.
#'
#' @importFrom pomp Csnippet
#' @return An object of class \sQuote{pomp}.
#' @examples
#' m1 <- haiti1()
#' @export

haiti1 <- function() {
  ## make components pomp object building
  rinit <- Csnippet("
    double pop = pop_0;
    double frac = pop / (S_0 + E_0 + I_0 + A_0 + R_0);
    S = nearbyint(frac * S_0);
    E = nearbyint(frac * E_0);
    I = nearbyint(frac * I_0);
    A = nearbyint(frac * A_0);
    R = nearbyint(frac * R_0);
    incid = 0;
  ")

  rproc <- Csnippet("
    // transition rates
    double S_rate[2], E_rate[3], I_rate[2], A_rate[2], R_rate[2];
    // transition numbers
    double S_trans[2], E_trans[3], I_trans[2], A_trans[2], R_trans[2];

    // seasonality
    double beta = beta1*seas1 + beta2*seas2 +
                  beta3*seas3 + beta4*seas4 +
                  beta5*seas5 + beta6*seas6;

    // population and births
    int pop = S + E + I + A + R;
    int births = rpois(mu * pop * dt);

    // force of infection
    double foi = pow(I, nu) * (beta / pop);

    // transition rate calculations
    S_rate[0] = foi;  // S -> E
    S_rate[1] = delta; // S -> death

    E_rate[0] = sigma * (1 - theta0); // E -> I
    E_rate[1] = sigma * theta0; // E -> A
    E_rate[2] = delta; // E -> death

    I_rate[0] = gamma; // I -> R
    I_rate[1] = delta; // I -> death

    A_rate[0] = gamma; // A -> R
    A_rate[1] = delta; // A -> death

    R_rate[0] = alpha; // R -> S
    R_rate[1] = delta; // R -> death

    // transition numbers
    reulermultinom(2, S, &S_rate[0], dt, &S_trans[0]);
    reulermultinom(3, E, &E_rate[0], dt, &E_trans[0]);
    reulermultinom(2, I, &I_rate[0], dt, &I_trans[0]);
    reulermultinom(2, A, &A_rate[0], dt, &A_trans[0]);
    reulermultinom(2, R, &R_rate[0], dt, &R_trans[0]);

    // states
    // S += - (S to E) - (S to death) + (R to S) + births
    S += -S_trans[0] - S_trans[1] + R_trans[0] + births;
    // E += - (E to I) - (E to A) - (E to death) + (S to E)
    E += -E_trans[0] - E_trans[1] - E_trans[2] + S_trans[0];
    // I += - (I to R) - (I to death) + (E to I)
    I += -I_trans[0] - I_trans[1] + E_trans[0];
    // A += - (A to R) - (A to death) + (E to R)
    A += -A_trans[0] - A_trans[1] + E_trans[1];
    // R += - (R to S) - (R to death) + (I to R) + (A to R)
    R += -R_trans[0] - R_trans[1] + I_trans[0] + A_trans[0];

    // incidence is cumulative entries into I state
    incid += E_trans[0];
  ")

  skel <- Csnippet('
    // transition rates
    double Srate[2];
    double Erate[3];
    double Irate[2];
    double Arate[2];
    double Rrate[2];

    // transition terms
    double Strans[2];
    double Etrans[3];
    double Itrans[2];
    double Atrans[2];
    double Rtrans[2];

    // some population demonitors
    int pop = S + E + I + A + R;
    double births = mu*pop;

    // make seasonal beta term for current time
    double mybeta = beta1*seas1 + beta2*seas2 + beta3*seas3 +
                    beta4*seas4 + beta5*seas5 + beta6*seas6;
    double foi = pow(I, nu) * mybeta / pop;

    //compute the rates for all the transitions
    Srate[0]= foi;  // S -> E
    Srate[1]= delta;

    Erate[0]= sigma*(1-theta0); // E -> I
    Erate[1]= sigma*theta0; // E -> A
    Erate[2]= delta;

    Irate[0]= gamma; // I -> R
    Irate[1]= delta;

    Arate[0]= gamma; // A -> R
    Arate[1]= delta;

    Rrate[0]= alpha; // R -> S exponential waning immunity from natural infection
    Rrate[1]= delta;

    // compute the transition terms
    for (int i = 0; i < 2; i++) {
      Strans[i] = Srate[i] * S;
    }
    for (int i = 0; i < 3; i++) {
      Etrans[i] = Erate[i] * E;
    }
    for (int i = 0; i < 2; i++) {
      Itrans[i] = Irate[i] * I;
    }
    for (int i = 0; i < 2; i++) {
      Atrans[i] = Arate[i] * A;
    }
    for (int i = 0; i < 2; i++) {
      Rtrans[i] = Rrate[i] * R;
    }

    // balance the equations
    DS = -Strans[0] - Strans[1] + Rtrans[0] + births;
    DE = -Etrans[0] - Etrans[1] - Etrans[2] + Strans[0];
    DI = -Itrans[0] - Itrans[1] + Etrans[0];
    DA = -Atrans[0] - Atrans[1] + Etrans[1];
    DR = -Rtrans[0] - Rtrans[1] + Itrans[0] + Atrans[0];
    Dincid = Etrans[0]; // incidence is cumulative entries into I state
  ')

  rmeas <- Csnippet("
    cases = rnbinom_mu(tau, rho*incid);
    if (cases > 0.0) {
      cases = nearbyint(cases);
    }
    else {
      cases = 0.0;
    }
  ")

  dmeas <- Csnippet("
    lik = dnbinom_mu(cases, tau, rho*incid, give_log);
  ")

  state_names <- c("S", "E", "I", "A", "R", "incid")

  param_names <- c("rho", "tau", "beta1", "beta2", "beta3",
                   "beta4", "beta5", "beta6", "gamma", "sigma",
                   "theta0", "alpha", "mu", "delta", "nu",
                   "S_0", "E_0", "I_0", "A_0", "R_0", "pop_0")

  param_trans <- pomp::parameter_trans(
    log = c("beta1", "beta2", "beta3", "beta4", "beta5", "beta6",
            "tau", "sigma", "gamma", "mu", "delta", "alpha"),
    logit = c("rho", "nu", "theta0"),
    barycentric = c("S_0", "E_0", "I_0", "A_0", "R_0")
  )

  ## get data
  dat <- haiti1_agg_data()

  ## make covariate table
  covar <- covars(tmin = 0,
                  tmax = nrow(dat) + 573, ## for 11 year forecast
                  byt = 1,
                  degree = 6,
                  nbasis = 6,
                  per = 52.14)

  ## build pomp model
  model1 <- pomp::pomp(
      data = dat,
      times = "week",
      t0 = 0,
      dmeasure = dmeas,
      rmeasure = rmeas,
      rprocess = pomp::euler(step.fun = rproc, delta.t = 1/7),
      skeleton = pomp::vectorfield(skel),
      covar = pomp::covariate_table(covar, times = "time"),
      partrans = param_trans,
      statenames = state_names,
      paramnames = param_names,
      accumvars = c("incid"),
      rinit = rinit
    )

  return(model1)
}
