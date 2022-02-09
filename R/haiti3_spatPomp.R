#' Build pomp object for model 3.
#'
#' Generate a \sQuote{spatPomp} object for fitting to Haiti cholera data.
#'    This model is a stochastic compartmental model applied at the
#'    level of the ten Haitian departments. It is the stochastic
#'    translation of a deterministic SIRB model based on Ordinary
#'    Differential Equations (ODEs), and has been implemented as a
#'    discrete-state model based on a Partially-Observed Markov Process
#'    (POMP), simulating the stochastic transitions between compartments
#'    as discrete events.
#'
#'    The model subdivides the population of each department into
#'    compartments counting the number of individuals at the different
#'    stages of the disease: Susceptible individuals (S),
#'    Infected symptomatic (I), infected Asymptomatic (A), and
#'    Recoverd (R). The main feature of this model is that it contains
#'    an environmental compartment describing the bacterial concentration (B)
#'    in the local environment, which is used to estimate the force of
#'    infection.
#'
#'    This model removes the redundancy in the Recovered compartments found
#'    in the original code.
#'
#'    This model was developed by Lemaitre, Joseph, et. al at the Laboratory
#'    of Ecohydrology, Ecole Polytechnique Federale de Lausanne (CH).
#'
#'
#'
#' @importFrom magrittr %>%
#' @importFrom foreach %do%
#' @importFrom foreach foreach
#' @import pomp
#' @import spatPomp
#' @seealso \code{\link{haiti2}} and \code{\link{haiti1}} for other models used
#' to fit the cholera epidemic in Haiti.
#' @return \code{\link[spatPomp]{spatPomp}} representation of model 3 described in \href{https://www.sciencedirect.com/science/article/pii/S2214109X20303107}{Lee, Elizabeth et. al.} and it's accompanying \href{https://ars.els-cdn.com/content/image/1-s2.0-S2214109X20303107-mmc3.pdf}{Supplemental Material}.
#'
#' @examples
#' mod3 <- haiti3_spatPomp()
#' @export

haiti3_spatPomp <- function(dt_years = 1/365.25) {

  # Create vector of departement names
  departements = c(
    'Artibonite', 'Centre', 'Grande_Anse',
    'Nippes', 'Nord', 'Nord_Est', 'Nord_Ouest',
    'Ouest', 'Sud', 'Sud_Est'
  )

  # First Define some helper functions:
  dateToYears <- function(date, origin = as.Date("2014-01-01"), yr_offset = 2014) {
    # This function converts a date to a decimal representation
    #
    # ex: "1976-03-01" -> 1976.163

    julian(date, origin = origin) / 365.25 + yr_offset
  }

  yearsToDate <- function(year_frac, origin = as.Date("2014-01-01"), yr_offset = 2014.0) {
    # This function is the inverse function of dateToYears; it takes
    # a decimal representation of a date and converts it into a Date.
    #
    # ex: 1976.163 -> "1976-03-01"

    as.Date((year_frac - yr_offset) * 365.25, origin = origin)
  }

  yearsToDateTime <- function(year_frac, origin = as.Date("2014-01-01"), yr_offset = 2014.0) {
    # Same as the function above, but a DateTime object rather than a Date
    # object.
    #
    # ex: 1976.163 -> "1976-03-01"
    as.POSIXct((year_frac - yr_offset) * 365.25 * 3600 * 24, origin = origin)
  }

  # List all state-names in pomp object:
  unit_state_names <- c(
    "S", "I", "A", "R_one", "R_two", "R_three",
    "VSd", "VR1d", "VR2d", "VR3d",
    "VSdd", "VR1dd", "VR2dd", "VR3dd",
    "VSd_alt", "VR1d_alt", "VR2d_alt", "VR3d_alt",
    "VSdd_alt", "VR1dd_alt",
    "VR2dd_alt", "VR3dd_alt", "C", "B"
  )

  all_state_names <- paste0(rep(unit_state_names, each = 10), 1:10)

  # All parameters that are common to each departement (disease specific parameters)
  params_common <- c(
    "sigma", "mu_B", "thetaI", "XthetaA", "lambdaR", "r",
    "gamma", "rho", "epsilon", "k",
    "std_W", "cas_def", "mu", "alpha", "cases_ext"
  )

  # Parameters that are unique to each department:
  params_diff <- c(
    "foi_add", "betaB", "H", "D", "t_vacc_start",
    "t_vacc_end", "p1d_reg", "r_v_year", "t_vacc_start_alt",
    "t_vacc_end_alt", "p1d_reg_alt", "r_v_year_alt", "Binit"
  )

  all_params_names <- c(params_common, params_diff)
  all_unit_params_names <- paste0(rep(all_params_names, each = 10), 1:10)
  all_unit_params <- rep(0, length(all_unit_params_names))
  names(all_unit_params) <- all_unit_params_names

  # Load the input parameters
  # t_start <- dateToYears(as.Date(MODEL3_INPUT_PARAMETERS$t_start))
  t_start <- dateToYears(as.Date("2010-10-23"))
  t_end   <- dateToYears(as.Date(haitipkg:::MODEL3_INPUT_PARAMETERS$t_end))

  # all_state_names <- c('DosesAll', 'CasesAll')
  # all_param_names <- params_common

  all_matrix_cases_at_t_start.string <- ""

  # haitiCholera is a saved data.frame in the package
  MODEL3_CASES <- haitiCholera %>%
    dplyr::rename(
      date = date_saturday, Grande_Anse = Grand.Anse,
      Nord_Est = Nord.Est, Nord_Ouest = Nord.Ouest,
      Sud_Est = Sud.Est
    ) %>%
    dplyr::mutate(date = as.Date(date)) %>%
    dplyr::select(-report)

  all_cases <- MODEL3_CASES %>%
    dplyr::mutate(
      date = as.Date(date, format = '%Y-%m-%d'),
      time = dateToYears(date)
    ) %>%
    tidyr::pivot_longer(
      data = .,
      cols = 2:11,
      names_to = "departement",
      values_to = "cases"
    ) %>%
    dplyr::arrange(departement, time) %>%
    dplyr::select(-date)

  std_rain <- function(x) {
    # This function simply standardizes the rain for us.
    x / max(x)
  }

  all_rain <- haitiRainfall %>%
    dplyr::summarize(
      date = date, dplyr::across(Artibonite:`Sud-Est`, std_rain)
    ) %>%
    dplyr::mutate(
      time = dateToYears(date)
    ) %>%
    dplyr::filter(time > t_start - 0.01 & time < (t_end + 0.01))


  colnames(all_rain) <- c(
    "date",
    paste0(
      'rain_std', c(
        'Artibonite', 'Centre', 'Grande_Anse',
        'Nippes', 'Nord', 'Nord_Est', 'Nord_Ouest',
        'Ouest', 'Sud', 'Sud_Est'
      )
    ),
    'time'
  )

  all_rain <- all_rain %>%
    tidyr::pivot_longer(
      data = .,
      cols = 2:11,
      names_to = "departement",
      values_to = "rain_std",
      names_prefix = "rain_std"
    ) %>%
    dplyr::arrange(departement, time) %>%
    dplyr::select(-date)


  all_cases_at_t_start.string <- ""
  for (i in 1:10) {  # Loop through data to get starting observations in each dep.

    # Select the departement
    dp <- departements[i]
    dep_cases <- all_cases %>%
      dplyr::filter(departement == dp)

    # Look at all observations prior to t_start
    cases_at_t_start <- dep_cases %>% dplyr::filter(time <= t_start)

    # Loop through each of the rows in the remaining data, and write the row
    # as an array in C.
    cases_at_t_start.string <- foreach::foreach(
      r = iterators::iter(cases_at_t_start, by = "row"),
      .combine = c
    ) %do% {
      sprintf("{%f, %f}", r$time, r$cases)
    } %>%
      stringr::str_c(collapse = ", \n") %>%
      paste0(" {", ., '}')  # Wrap all rows in {}, so single object for each dep

    # special cases for starting the array, ending it, or just in the middle.
    if (i == 1) {  # Start
      cases_at_t_start.string <- paste0(
        sprintf(
          "double cases_at_t_start[10][%i][%i] = {\n",
          nrow(cases_at_t_start),
          2
        ), cases_at_t_start.string, ',\n'
      )
    } else if (i == 10) {  # End
      cases_at_t_start.string <- paste0(
        cases_at_t_start.string, "\n};"
      )
    } else {  # Middle
      cases_at_t_start.string <- paste0(
        cases_at_t_start.string, ",\n"
      )
    }

    # all_matrix_cases_at_t_start.string <- stringr::str_c(all_matrix_cases_at_t_start.string, matrix_cases_at_t_start.string)
    all_cases_at_t_start.string <- paste0(all_cases_at_t_start.string, cases_at_t_start.string)
  }

  MODEL3_INPUT_PARAMETERS <- haitipkg:::MODEL3_INPUT_PARAMETERS

  populations <- unlist(purrr::flatten(MODEL3_INPUT_PARAMETERS["population"]))
  densities <- unlist(purrr::flatten(MODEL3_INPUT_PARAMETERS["density"]))
  p1d_alt_year <- unlist(purrr::flatten(MODEL3_INPUT_PARAMETERS['p1d_alt_year']))
  nb_doses_alt_year <- unlist(purrr::flatten(MODEL3_INPUT_PARAMETERS['nb_doses_alt_year']))
  t_vacc_start_alt  <- unlist(purrr::flatten(MODEL3_INPUT_PARAMETERS["t_vacc_start_alt"]))
  t_vacc_end_alt <- unlist(purrr::flatten(MODEL3_INPUT_PARAMETERS["t_vacc_end_alt"]))

  names(populations) <- gsub("-", "_", names(populations))
  names(densities) <- gsub("-", "_", names(densities))
  names(p1d_alt_year) <- gsub("-", "_", names(p1d_alt_year))
  names(nb_doses_alt_year) <- gsub("-", "_", names(nb_doses_alt_year))
  names(t_vacc_start_alt) <- gsub("-", "_", names(t_vacc_start_alt))
  names(t_vacc_end_alt) <- gsub("-", "_", names(t_vacc_end_alt))

  for (i in 1:10) {
    dp <- departements[i]

    all_unit_params[paste0('H', i)] <- populations[dp]
    all_unit_params[paste0('D', i)] <- densities[dp]

    # Vaccination information:
    t_vacc_start_alt_dep = dateToYears(as.Date(t_vacc_start_alt[dp]))
    t_vacc_end_alt_dep   = dateToYears(as.Date(t_vacc_end_alt[dp]))
    r_v_alt_year_dep = nb_doses_alt_year[dp] / (t_vacc_end_alt_dep - t_vacc_start_alt_dep)
    p1d_alt_dep = p1d_alt_year[dp]

    all_unit_params[paste0("t_vacc_start_alt", i)] = t_vacc_start_alt_dep
    all_unit_params[paste0("t_vacc_end_alt", i)] = t_vacc_end_alt_dep
    all_unit_params[paste0("p1d_reg_alt", i)] = p1d_alt_dep
    all_unit_params[paste0("r_v_year_alt", i)] = r_v_alt_year_dep
  }

  initalizeStatesTemplate =   "
  double *S = &S1;
  double *I = &I1;
  double *A = &A1;
  double *R_one = &R_one1;
  double *R_two = &R_two1;
  double *R_three = &R_three1;
  double *VSd = &VSd1;
  double *VR1d = &VR1d1;
  double *VR2d = &VR2d1;
  double *VR3d = &VR3d1;
  double *VSdd = &VSdd1;
  double *VR1dd = &VR1dd1;
  double *VR2dd = &VR2dd1;
  double *VR3dd = &VR3dd1;
  double *VSd_alt = &VSd_alt1;
  double *VR1d_alt = &VR1d_alt1;
  double *VR2d_alt = &VR2d_alt1;
  double *VR3d_alt = &VR3d_alt1;
  double *VSdd_alt = &VSdd_alt1;
  double *VR1dd_alt = &VR1dd_alt1;
  double *VR2dd_alt = &VR2dd_alt1;
  double *VR3dd_alt = &VR3dd_alt1;
  double *C = &C1;
  double *B = &B1;
  const double *thetaI = &thetaI1;
  const double *XthetaA = &XthetaA1;
  const double *mu = &mu1;
  const double *sigma = &sigma1;
  const double *epsilon = &epsilon1;
  const double *alpha = &alpha1;
  const double *gamma = &gamma1;
  const double *rho = &rho1;
  const double *Binit = &Binit1;
  const double *r = &r1;
  const double *mu_B = &mu_B1;
  const double *H = &H1;
  const double *D = &D1;
  const double *lambdaR = &lambdaR1;


  double R0[2] = {0,0};

  for (int u = 0; u < U; u++) {

    double thetaA = thetaI[u] * XthetaA[u];

    A[u] = nearbyint((1-sigma[u])/sigma[u] * 1/epsilon[u] * cases_at_t_start[u][n_cases_start-1][1]/7 * 365 /(mu[u] + gamma[u]));
    I[u] = nearbyint(1/epsilon[u] * cases_at_t_start[u][n_cases_start-1][1]/7 * 365 /(mu[u] + alpha[u] + gamma[u]));  // Steady state

    // Reset R0 vector
    R0[0] = 0;
    R0[1] = 0;

    for(int i = 0; i < n_cases_start; i++){
      R0[0] +=                   cases_at_t_start[u][i][1]/epsilon[u]  * exp((cases_at_t_start[u][i][0] - t_start)  * (rho[u] + mu[u])); /* because t_i in past so t_ - t_0 negative */
      R0[1] += (1-sigma[u])/sigma[u] * cases_at_t_start[u][i][1]/epsilon[u]  * exp((cases_at_t_start[u][i][0] - t_start)  * (rho[u] + mu[u]));
    }

    R_one[u] = nearbyint((R0[0] + R0[1]) / 3);
    R_two[u] = nearbyint((R0[0] + R0[1]) / 3);
    R_three[u] = nearbyint((R0[0] + R0[1]) / 3);

    if (A[u] + I[u] + R_one[u] + R_two[u] + R_three[u] >= H[u]) {
      double R_tot = H[u] - A[u] - I[u] - 100.0;
      if (R_tot <= 0)
      {
      I[u]     = nearbyint(H[u] - 100);
      A[u]     = nearbyint(0);
      R_tot = nearbyint(0);
      }
      R_one[u] = nearbyint(R_tot / 3);
      R_two[u] = nearbyint(R_tot / 3);
      R_three[u] = nearbyint(R_tot / 3);
    }

    S[u]   = nearbyint(H[u] - A[u] - I[u] - R_one[u] - R_two[u] - R_three[u]);
    B[u]   = (I[u] * thetaI[u]/mu_B[u] + A[u] * thetaA/mu_B[u]) * D[u] * (1 + lambdaR[u] * pow(Binit[u], r[u]));
    C[u]   = 0;

    VSd[u] = 0;
    VR1d[u] = 0;
    VR2d[u] = 0;
    VR3d[u] = 0;

    VSdd[u] = 0;
    VR1dd[u] = 0;
    VR2dd[u] = 0;
    VR3dd[u] = 0;

    VSd_alt[u] = 0;
    VR1d_alt[u] = 0;
    VR2d_alt[u] = 0;
    VR3d_alt[u] = 0;

    VSdd_alt[u] = 0;
    VR1dd_alt[u] = 0;
    VR2dd_alt[u] = 0;
    VR3dd_alt[u] = 0;

  }
"

  initalizeStates <- pomp::Csnippet(initalizeStatesTemplate)

  ###
  ### rmeas
  ###

  rmeasTemplate <- "
  const double *C = &C1;
  double *cases = &cases1;
  const double *epsilon = &epsilon1;
  const double *k = &k1;
  const double *cas_def = &cas_def1;
  int u;

  for (u = 0; u < U; u++) {
    if (t > 2018) {
      cases[u] = rnbinom_mu(k[u], epsilon[u] * C[u] * cas_def[u]);
    } else {
      cases[u] = rnbinom_mu(k[u], epsilon[u] * C[u]);
    }
  }
  "

  rmeas <- pomp::Csnippet(rmeasTemplate)

  ###
  ### dmeas
  ###

  dmeasTemplate <- "
int u;
const double *epsilon = &epsilon1;
const double *k = &k1;
const double *cas_def = &cas_def1;
const double *cases = &cases1;
const double *C = &C1;

lik = 0;
for (u = 0; u < U; u++) {
  if (ISNA(cases[u])) {
      lik += (give_log) ? 0 : 1;
  } else {
    if (t > 2018) {
      lik += dnbinom_mu(cases[u], k[u], epsilon[u] * C[u] * cas_def[u], give_log);
    } else {
      lik += dnbinom_mu(cases[u], k[u], epsilon[u] * C[u], give_log);
      }
    }
}
"

dmeas <- pomp::Csnippet(dmeasTemplate)


# TODO: Below assumes that k, epsilon, and cas_def constant across units.
unit_dmeasTemplate <- "
  double k = k1;
  double epsilon = epsilon1;
  double cas_def = cas_def1;

  if (ISNA(cases)) {
      lik = (give_log) ? 0 : 1;
  } else {
    if (t > 2018) {
      lik = dnbinom_mu(cases, k, epsilon * C * cas_def, give_log);
    } else {
      lik = dnbinom_mu(cases, k, epsilon * C, give_log);
      }
    }
"

unit_dmeas <- pomp::Csnippet(unit_dmeasTemplate)

###
### rproc
###

# Every department gets a copy.
rprocTemplate <- "

double *S = &S1;
double *I = &I1;
double *A = &A1;
double *R_one = &R_one1;
double *R_two = &R_two1;
double *R_three = &R_three1;
double *VSd = &VSd1;
double *VR1d = &VR1d1;
double *VR2d = &VR2d1;
double *VR3d = &VR3d1;
double *VSdd = &VSdd1;
double *VR1dd = &VR1dd1;
double *VR2dd = &VR2dd1;
double *VR3dd = &VR3dd1;
double *VSd_alt = &VSd_alt1;
double *VR1d_alt = &VR1d_alt1;
double *VR2d_alt = &VR2d_alt1;
double *VR3d_alt = &VR3d_alt1;
double *VSdd_alt = &VSdd_alt1;
double *VR1dd_alt = &VR1dd_alt1;
double *VR2dd_alt = &VR2dd_alt1;
double *VR3dd_alt = &VR3dd_alt1;
double *C = &C1;
double *B = &B1;
const double *rain_std = &rain_std1;

// getting all non-constant parameters used in the model
const double *betaB = &betaB1;
const double *foi_add = &foi_add1;
const double *H = &H1;
const double *D = &D1;
const double *t_vacc_start = &t_vacc_start1;
const double *t_vacc_end = &t_vacc_end1;
const double *p1d_reg = &p1d_reg1;
const double *r_v_year = &r_v_year1;
const double *t_vacc_start_alt = &t_vacc_start_alt1;
const double *t_vacc_end_alt = &t_vacc_end_alt1;
const double *p1d_reg_alt = &p1d_reg_alt1;
const double *r_v_year_alt = &r_v_year_alt1;

// Below I'm assuming that all these parameters are constant across the units.
const double thetaI = thetaI1;
const double XthetaA = XthetaA1;
const double mu = mu1;
const double sigma = sigma1;
const double epsilon = epsilon1;
const double alpha = alpha1;
const double gamma = gamma1;
const double rho = rho1;
// const double *Binit = &Binit1;  IVP
const double r = r1;
const double mu_B = mu_B1;
const double *cases_ext = &cases_ext1;
const double lambdaR = lambdaR1;
const double std_W = std_W1;
const double k = k1;
const double cas_def = cas_def1;

double foi, foi_stoc;   // force of infection and its stochastic version
double dw;              // extra-demographic stochasticity on foi
double dB;              // deterministic forward time difference of bacteria in the environment
double RK1, RK2, RK3, RK4;  // coefficients of  the Runge-Kutta method
double rate[48];        // vector of all rates in model
double dN[60];          // vector of transitions between classes during integration timestep
double mobility;
double p1d, pdd;
double r_v_wdn = 0.0;   // rate of vaccination: 0 if out of time window, r_v if not
int previous_vacc_campaign; // flag that indicate if we are on the first or second campain
double t_eff, t_eff_alt;

// Define rates that are constant for all units (departements):
// S compartment
rate[4] = mu;           // S -> natural death

// I compartment
rate[5] = mu;           // natural deaths
rate[6] = alpha;        // cholera-induced deaths
rate[7] = gamma;        // I -> R

// A compartment
rate[8] = mu;           // natural death
rate[9] = gamma;        // A -> R_one

// R_one
rate[12] = 3 * rho;     // loss of natural immunity; R_one -> R_two
rate[13] = mu;          // natural death R_one -> death

// R_two
rate[16] = 3 * rho;     // loss of natural immunity; R_two -> R_three
rate[17] = mu;          // natural death R_two -> death

// R_three
rate[20] = 3 * rho;     // loss of natural immunity; R_three -> S
rate[21] = mu;          // natural death R_three -> death

// VSd
rate[26] = mu;          // VSd -> death

// VR1d
rate[27] = mu;          // natural death:    VR1d -> death
rate[28] = 3 * rho;     // loss of immunity: VR1d -> VR2d

// VR2d
rate[29] = mu;          // natural death:    VR2d -> death
rate[30] = 3 * rho;     // loss of immunity: VR2d -> VR3d

// VR3d
rate[31] = mu;          // natural death:    VR3d -> death
rate[32] = 3 * rho;     // loss of immunity: VR3d -> VSd

// VSdd
rate[35] = mu;          // natural death

// VR1dd
rate[36] = mu;          // natural death:    VR1dd -> death
rate[37] = 3 * rho;     // loss of immunity: VR1dd -> VR2dd

// VR2dd
rate[38] = mu;          // natural death:    VR2dd -> death
rate[39] = 3 * rho;     // loss of immunity: VR2dd -> VR3dd

// VR3dd
rate[40] = mu;          // natural death:    VR3dd -> death
rate[41] = 3 * rho;     // loss of immunity: VR3dd -> VSdd

// VSd_alt
rate[44] = mu;          // natural death

// VSdd_alt
rate[47] = mu;          // natural death

// Loop through each unit (departement)
for (int u = 0; u < U; u++) {
  int scenario =  cases_ext[u];
  double thetaA = thetaI * XthetaA;

  previous_vacc_campaign = TRUE;
  r_v_wdn = 0;
  dw = 0;

  mobility = I[0] + I[1] + I[2] + I[3] + I[4] + I[5] + I[6] + I[7] + I[8] + I[9] +
             A[0] + A[1] + A[2] + A[3] + A[4] + A[5] + A[6] + A[7] + A[8] + A[9] -
             (I[u] + A[u]);

  // force of infection
  foi = betaB[u] * (B[u] / (1 + B[u])) + foi_add[u] * mobility;

  if(std_W > 0.0) {
    dw = rgammawn(std_W, dt);   // white noise (extra-demographic stochasticity)
    foi_stoc = foi * dw/dt;     // apply stochasticity
  } else {
    foi_stoc = foi;
  }


  if (t <= (t_vacc_end_alt[u] + dt)){
	  previous_vacc_campaign = TRUE;

	  if (t >= t_vacc_start_alt[u] && t <= (t_vacc_end_alt[u] + dt)) {
	    r_v_wdn = (r_v_year_alt[u] / (S[u] + A[u] + R_one[u] + R_two[u] + R_three[u]));
	  }

	  p1d = p1d_reg_alt[u];
  } else {
	  previous_vacc_campaign = FALSE;

	  if (t >= t_vacc_start[u] && t <= (t_vacc_end[u] + dt)) {
    	r_v_wdn = (r_v_year[u] / (S[u] + A[u] + R_one[u] + R_two[u] + R_three[u]));
	  }
	  p1d = p1d_reg[u];
  }

  pdd = 1 - p1d;

  // time in the vacc_eff referential. We assume different timing for 1d and 2d
  t_eff =     t - (t_vacc_start[u] + (t_vacc_end[u] - t_vacc_start[u])/2);
  t_eff_alt = t - (t_vacc_start_alt[u] + (t_vacc_end_alt[u] - t_vacc_start_alt[u])/2);


  // define transition rates for each type of event (i.e what multplies the thing)

  // S compartment
  rate[0] = sigma * foi_stoc;         // infections
  rate[1] = (1 - sigma) * foi_stoc;   // asymptomatic infections
  rate[2] = p1d * r_v_wdn;            // S -> VSd
  rate[3] = pdd * r_v_wdn;            // S -> VSdd

  // A compartment
  rate[10] = p1d * r_v_wdn;   // A -> VR1d
  rate[11] = pdd * r_v_wdn;   // A -> VR1dd

  // R_one
  rate[14] = p1d * r_v_wdn;  // R_one -> VR1d
  rate[15] = pdd * r_v_wdn;  // R_one -> VR1dd

  // R_two
  rate[18] = p1d * r_v_wdn;  // R_two -> VR2d
  rate[19] = pdd * r_v_wdn;  // R_two -> VR2dd

  // R_three
  rate[22] = p1d * r_v_wdn;  // R_three -> VR3d
  rate[23] = pdd * r_v_wdn;  // R_three -> VR3dd

  // VSd
  rate[24] = sigma * (1 - eff_v_1d(t_eff, scenario)) * foi_stoc;       // VSd -> I
  rate[25] = (1 - sigma) * (1 - eff_v_1d(t_eff, scenario)) * foi_stoc; // VSd -> A

  // VSdd
  rate[33] = sigma * (1 - eff_v_2d(t_eff, scenario)) * foi_stoc;       // VSdd -> I
  rate[34] = (1 - sigma) * (1 - eff_v_2d(t_eff, scenario)) * foi_stoc; // VSdd -> A

  /* For previous vacc campagain */

  // VSd_alt
  rate[42] = sigma       * (1 - eff_v_1d(t_eff_alt, scenario)) * foi_stoc; // VSd -> I
  rate[43] = (1 - sigma) * (1 - eff_v_1d(t_eff_alt, scenario)) * foi_stoc; // VSd -> A

  // VSdd_alt
  rate[45] = sigma       * (1 - eff_v_2d(t_eff_alt, scenario)) * foi_stoc; // VSdd -> I
  rate[46] = (1 - sigma) * (1 - eff_v_2d(t_eff_alt, scenario)) * foi_stoc; // VSdd -> A

  // All other alts already defined.

  // simulate all transitions
  reulermultinom(5,  S[u],     &rate[0],  dt, &dN[0]);
  reulermultinom(3,  I[u],     &rate[5],  dt, &dN[5]);
  reulermultinom(4,  A[u],     &rate[8],  dt, &dN[8]);
  reulermultinom(4, R_one[u],    &rate[12], dt, &dN[12]);
  reulermultinom(4, R_two[u],    &rate[16], dt, &dN[16]);
  reulermultinom(4, R_three[u],    &rate[20], dt, &dN[20]);

  // Vaccinated 1 dose
  reulermultinom(3,  VSd[u],   &rate[24], dt, &dN[24]);
  reulermultinom(2, VR1d[u],   &rate[27], dt, &dN[27]);
  reulermultinom(2, VR2d[u],   &rate[29], dt, &dN[29]);
  reulermultinom(2, VR3d[u],   &rate[31], dt, &dN[31]);

  // Vaccinated 2 doses
  reulermultinom(3,  VSdd[u],   &rate[33], dt, &dN[33]);
  reulermultinom(2, VR1dd[u],   &rate[36], dt, &dN[36]);
  reulermultinom(2, VR2dd[u],   &rate[38], dt, &dN[38]);
  reulermultinom(2, VR3dd[u],   &rate[40], dt, &dN[40]);

  /* For the previous vaccination campain */

  // Alt Vaccinated 1 dose
  reulermultinom(3,  VSd_alt[u],   &rate[42], dt, &dN[42]);
  reulermultinom(2, VR1d_alt[u],   &rate[27], dt, &dN[45]);
  reulermultinom(2, VR2d_alt[u],   &rate[29], dt, &dN[47]);
  reulermultinom(2, VR3d_alt[u],   &rate[31], dt, &dN[49]);

  // Alt Vaccinated 2 doses
  reulermultinom(3,  VSdd_alt[u],   &rate[45], dt, &dN[51]);
  reulermultinom(2, VR1dd_alt[u],   &rate[27], dt, &dN[54]);
  reulermultinom(2, VR2dd_alt[u],   &rate[29], dt, &dN[56]);
  reulermultinom(2, VR3dd_alt[u],   &rate[31], dt, &dN[58]);

  // bacteria as continous state variable
  // implement Runge-Kutta integration assuming S, I, R, V* stay constant during dt
  RK1 = dt * fB(I[u], A[u], B[u], mu_B, thetaI, thetaA, lambdaR, rain_std[u], r, D[u]);
  RK2 = dt * fB(I[u], A[u], B[u], mu_B, thetaI, thetaA, lambdaR, rain_std[u], r, D[u]);
  RK3 = dt * fB(I[u], A[u], B[u], mu_B, thetaI, thetaA, lambdaR, rain_std[u], r, D[u]);
  RK4 = dt * fB(I[u], A[u], B[u], mu_B, thetaI, thetaA, lambdaR, rain_std[u], r, D[u]);

  // bacteria increment
  dB = (RK1 + 2*RK2 + 2*RK3 + RK4) / 6.0;

  // Update States
  I[u]  += dN[0] + dN[24] + dN[33] + dN[42] + dN[51] - dN[7] - dN[6] - dN[5];            // S -> I, VSd -> I, VSdd -> I, VSd_alt -> I, VSdd_alt -> I, I -> R, I-> death, I -> death
  A[u]  += dN[1] + dN[25] + dN[34] + dN[43] + dN[52] - dN[9] - dN[10] - dN[11] - dN[8];  // S -> A, VSd -> A, VSdd -> A, VSd_alt -> A, VSdd_alt -> A, A -> R_one, A -> VR1d, A -> VR1dd, natural death.
  R_one[u] += dN[7] + dN[9] - dN[12] - dN[14] - dN[15] - dN[13];                            // I-> R_one, A -> R_one, R_one -> R_two, R_one -> VR1d, R_one -> VR1dd, R_one -> death
  R_two[u] += dN[12] - dN[16] - dN[18] - dN[19] - dN[17];                                   // R_one -> R_two, R_two -> R_three, R_two -> VR2d, R_two -> VR2dd, R_two -> death
  R_three[u] += dN[16] - dN[20] - dN[22] - dN[23] - dN[21];                                   // R_two -> R_three, R_three -> S , R_three -> VR3d, R_three -> VR3dd, R_three -> death


  if (previous_vacc_campaign){
  	VSd_alt[u]    +=  dN[2];            // S -> VSd_alt
  	VR1d_alt[u]   +=  dN[14] + dN[10];  // R_one -> VR1d_alt, A -> VR1d_alt
  	VR2d_alt[u]   +=  dN[18];           // R_two -> VR2d_alt
  	VR3d_alt[u]   +=  dN[22];           // R_three -> VR3d_alt

  	VSdd_alt[u]   += dN[3];             // S -> VSdd_alt
  	VR1dd_alt[u]  += dN[15] + dN[11];   // R_one -> VR1dd_alt, A -> VR1dd_alt
  	VR2dd_alt[u]  += dN[19];            // R_two -> VR2dd_alt
  	VR3dd_alt[u]  += dN[23];            // R_three -> VR3dd_alt

  } else {
  	VSd[u]    += dN[2];                 // S -> VSd
  	VR1d[u]   += dN[14] + dN[10];       // R_one -> VR1d, A -> VR1d
  	VR2d[u]   += dN[18];                // R_two -> VR2d
  	VR3d[u]   += dN[22];                // R_three -> VR3d

  	VSdd[u]   += dN[3];                 // S -> VSdd
  	VR1dd[u]  += dN[15] + dN[11];       // R_one -> VR1dd, A -> VR1dd
  	VR2dd[u]  += dN[19];                // R_two -> VR2dd
  	VR3dd[u]  += dN[23];                // R_three -> VR3dd
  }

  VSd[u]   += dN[32] - dN[24] - dN[25] - dN[26]; // VR3d -> VSd, VSd -> I, VSd -> A, VSd -> death
  VR1d[u]  += - dN[27] - dN[28];                 // VR1d -> death, VR1d -> VR2d
  VR2d[u]  += dN[28] - dN[29] - dN[30];          // VR1d -> VR2d, VR2d -> death, VR2d -> VR3d
  VR3d[u]  += dN[30] - dN[31] - dN[32];          // VR2d -> VR3d, VR3d -> death, VR3d -> VSd

  VSdd[u]   +=  dN[41] - dN[33] - dN[34] - dN[35]; // VR3dd -> VSdd, VSdd -> I, VSdd -> A, VSdd -> death
  VR1dd[u]  += - dN[36] - dN[37];                  // VR1dd -> death, VR1dd -> VR2dd
  VR2dd[u]  +=  dN[37] - dN[38] - dN[39];          // VR1dd -> VR2dd, VR2dd -> death, VR2dd -> VR3dd
  VR3dd[u]  +=  dN[39] - dN[40] - dN[41];          // VR2dd -> VR3dd, VR3dd -> death, VR3dd -> VSdd

  // previous vacccination campain
  VSd_alt[u]   += dN[50] - dN[42] - dN[43] - dN[44]; // VR3d_alt -> VSd_alt, VSd_alt -> I, VSd_alt -> A, VSd_alt -> death
  VR1d_alt[u]  += - dN[45] - dN[46];                 // death, VR1d_alt -> VR2d_alt
  VR2d_alt[u]  += dN[46] - dN[47] - dN[48];          // VR1d_alt -> VR2d_alt, death, VR2d_alt -> VR3d_alt
  VR3d_alt[u]  += dN[48] - dN[49] - dN[50];          // VR2d_alt -> VR3d_alt, death, VR3d_alt -> VSd_alt

  VSdd_alt[u]   +=  dN[59] - dN[51] - dN[52] - dN[53]; // VR3dd_alt -> VSdd_alt, VSdd_alt -> I, VSdd_alt -> A, VSdd_alt -> death
  VR1dd_alt[u]  += - dN[54] - dN[55];
  VR2dd_alt[u]  +=  dN[55] - dN[56] - dN[57] ;
  VR3dd_alt[u]  +=  dN[57] - dN[58] - dN[59];

  C[u]   +=  dN[0] + dN[24] + dN[33] + dN[42] + dN[51]; // S -> I, VSd -> I, VSdd -> I, VSd_alt -> I, VSdd_alt -> I
  B[u]   += (((dB) < -B[u]) ? (-B[u] + 1.0e-3) : (dB)); // condition to ensure B>0

  // susceptibles so as to match total population
  S[u] = nearbyint(H[u] - I[u] - A[u] - R_one[u] - R_two[u] - R_three[u] -
  	VSd[u] - VR1d[u] - VR2d[u] - VR3d[u] -
  	VSdd[u] - VR1dd[u] -VR2dd[u] -VR3dd[u] -
  	VSd_alt[u] - VR1d_alt[u] - VR2d_alt[u] - VR3d_alt[u] -
  	VSdd_alt[u] - VR1dd_alt[u] - VR2dd_alt[u] - VR3dd_alt[u]);


  // if (!previous_vacc_campaign) {
  //	DosesAll  += dN[2] + dN[10] + dN[14] + dN[18] + dN[22] + 2 * (dN[3] + dN[11] + dN[15] + dN[19] + dN[23]);
  // }

  // CasesAll  +=  dN[0] + dN[24] + dN[33] + dN[42] + dN[51];

}  // end of u-loop
  "

final_rproc_c <- pomp::Csnippet(rprocTemplate)

# C function to compute the time-derivative of bacterial concentration OK
derivativeBacteria.c <- " double fB(int I, int A, double B,
double mu_B, double thetaI, double XthetaA, double lambdaR, double rain, double r, double D) {
  double thetaA = thetaI * XthetaA;
  double dB;
  dB = -mu_B * B +  (1 + lambdaR * pow(rain, r)) * D * (thetaI * (double) I + thetaA * (double) A);
  return(dB);
};
  "

eff_v.c <- "

 double eff_v_2d(double t_since_vacc, int scenario) {
  double eff_v_2d = 0.0;
  switch(scenario){
  case 1:
      if      (t_since_vacc <=   1./12) eff_v_2d =  0.76              ;
      else if (t_since_vacc <=   2./12) eff_v_2d =  0.753527484533759;
      else if (t_since_vacc <=   3./12) eff_v_2d =  0.746961516042262;
      else if (t_since_vacc <=   4./12) eff_v_2d =  0.740300745209625;
      else if (t_since_vacc <=   5./12) eff_v_2d =  0.733543803237948;
      else if (t_since_vacc <=   6./12) eff_v_2d =  0.726689301566023;
      else if (t_since_vacc <=   7./12) eff_v_2d =  0.719735831583988;
      else if (t_since_vacc <=   8./12) eff_v_2d =  0.712681964343848;
      else if (t_since_vacc <=   9./12) eff_v_2d =  0.705526250265831;
      else if (t_since_vacc <=  10./12) eff_v_2d =  0.698267218840495;
      else if (t_since_vacc <=  11./12) eff_v_2d =  0.690903378326535;
      else if (t_since_vacc <=  12./12) eff_v_2d =  0.683433215444227;
      else if (t_since_vacc <=  13./12) eff_v_2d =  0.675855195064453;
      else if (t_since_vacc <=  14./12) eff_v_2d =  0.668167759893222;
      else if (t_since_vacc <=  15./12) eff_v_2d =  0.660369330151649;
      else if (t_since_vacc <=  16./12) eff_v_2d =  0.652458303251305;
      else if (t_since_vacc <=  17./12) eff_v_2d =  0.644433053464886;
      else if (t_since_vacc <=  18./12) eff_v_2d =  0.636291931592122;
      else if (t_since_vacc <=  19./12) eff_v_2d =  0.628033264620864;
      else if (t_since_vacc <=  20./12) eff_v_2d =  0.619655355383277;
      else if (t_since_vacc <=  21./12) eff_v_2d =  0.61115648220707 ;
      else if (t_since_vacc <=  22./12) eff_v_2d =  0.602534898561692;
      else if (t_since_vacc <=  23./12) eff_v_2d =  0.593788832699414;
      else if (t_since_vacc <=  24./12) eff_v_2d =  0.584916487291234;
      else if (t_since_vacc <=  25./12) eff_v_2d =  0.575916039057525;
      else if (t_since_vacc <=  26./12) eff_v_2d =  0.566785638393345;
      else if (t_since_vacc <=  27./12) eff_v_2d =  0.557523408988343;
      else if (t_since_vacc <=  28./12) eff_v_2d =  0.548127447441173;
      else if (t_since_vacc <=  29./12) eff_v_2d =  0.538595822868346;
      else if (t_since_vacc <=  30./12) eff_v_2d =  0.528926576507423;
      else if (t_since_vacc <=  31./12) eff_v_2d =  0.519117721314497;
      else if (t_since_vacc <=  32./12) eff_v_2d =  0.509167241555842;
      else if (t_since_vacc <=  33./12) eff_v_2d =  0.499073092393685;
      else if (t_since_vacc <=  34./12) eff_v_2d =  0.488833199465984;
      else if (t_since_vacc <=  35./12) eff_v_2d =  0.478445458460148;
      else if (t_since_vacc <=  36./12) eff_v_2d =  0.467907734680592;
      else if (t_since_vacc <=  37./12) eff_v_2d =  0.457217862610059;
      else if (t_since_vacc <=  38./12) eff_v_2d =  0.446373645464601;
      else if (t_since_vacc <=  39./12) eff_v_2d =  0.435372854742138;
      else if (t_since_vacc <=  40./12) eff_v_2d =  0.424213229764494;
      else if (t_since_vacc <=  41./12) eff_v_2d =  0.412892477212831;
      else if (t_since_vacc <=  42./12) eff_v_2d =  0.401408270656362;
      else if (t_since_vacc <=  43./12) eff_v_2d =  0.38975825007427  ;
      else if (t_since_vacc <=  44./12) eff_v_2d =  0.377940021370718;
      else if (t_since_vacc <=  45./12) eff_v_2d =  0.365951155882864;
      else if (t_since_vacc <=  46./12) eff_v_2d =  0.353789189881759;
      else if (t_since_vacc <=  47./12) eff_v_2d =  0.341451624066056;
      else if (t_since_vacc <=  48./12) eff_v_2d =  0.328935923048392;
      else if (t_since_vacc <=  49./12) eff_v_2d =  0.316239514834368;
      else if (t_since_vacc <=  50./12) eff_v_2d =  0.303359790293999;
      else if (t_since_vacc <=  51./12) eff_v_2d =  0.290294102625533;
      else if (t_since_vacc <=  52./12) eff_v_2d =  0.27703976681153  ;
      else if (t_since_vacc <=  53./12) eff_v_2d =  0.263594059067087;
      else if (t_since_vacc <=  54./12) eff_v_2d =  0.249954216280098;
      else if (t_since_vacc <=  55./12) eff_v_2d =  0.236117435443426;
      else if (t_since_vacc <=  56./12) eff_v_2d =  0.222080873078887;
      else if (t_since_vacc <=  57./12) eff_v_2d =  0.207841644652907;
      else if (t_since_vacc <=  58./12) eff_v_2d =  0.193396823983748;
      else if (t_since_vacc <=  59./12) eff_v_2d =  0.178743442640173;
      else if (t_since_vacc <=  60./12) eff_v_2d = 0.163878489331427;
      else if (t_since_vacc <=  61./12) eff_v_2d =  0.148798909288418;
   break;
    case 2:
      eff_v_2d = 0.76;
    break;
    case 3:
      if      (t_since_vacc <=   1./12) eff_v_2d = 0.62899141753524;
      else if (t_since_vacc <=   2./12) eff_v_2d = 0.62363463243243;
      else if (t_since_vacc <=   3./12) eff_v_2d = 0.61820050371012;
      else if (t_since_vacc <=   4./12) eff_v_2d = 0.61268791464710;
      else if (t_since_vacc <=   5./12) eff_v_2d = 0.60709573239845;
      else if (t_since_vacc <=   6./12) eff_v_2d = 0.60142280776277;
      else if (t_since_vacc <=   7./12) eff_v_2d = 0.59566797494594;
      else if (t_since_vacc <=   8./12) eff_v_2d = 0.58983005132162;
      else if (t_since_vacc <=   9./12) eff_v_2d = 0.58390783718819;
      else if (t_since_vacc <=  10./12) eff_v_2d = 0.57790011552220;
      else if (t_since_vacc <=  11./12) eff_v_2d = 0.57180565172828;
      else if (t_since_vacc <=  12./12) eff_v_2d = 0.56562319338543;
      else if (t_since_vacc <=  13./12) eff_v_2d = 0.55935146998966;
      else if (t_since_vacc <=  14./12) eff_v_2d = 0.55298919269287;
      else if (t_since_vacc <=  15./12) eff_v_2d = 0.54653505403800;
      else if (t_since_vacc <=  16./12) eff_v_2d = 0.53998772769036;
      else if (t_since_vacc <=  17./12) eff_v_2d = 0.53334586816505;
      else if (t_since_vacc <=  18./12) eff_v_2d = 0.52660811055048;
      else if (t_since_vacc <=  19./12) eff_v_2d = 0.51977307022784;
      else if (t_since_vacc <=  20./12) eff_v_2d = 0.51283934258662;
      else if (t_since_vacc <=  21./12) eff_v_2d = 0.50580550273589;
      else if (t_since_vacc <=  22./12) eff_v_2d = 0.49867010521154;
      else if (t_since_vacc <=  23./12) eff_v_2d = 0.49143168367921;
      else if (t_since_vacc <=  24./12) eff_v_2d = 0.48408875063295;
      else if (t_since_vacc <=  25./12) eff_v_2d = 0.47663979708957;
      else if (t_since_vacc <=  26./12) eff_v_2d = 0.46908329227848;
      else if (t_since_vacc <=  27./12) eff_v_2d = 0.46141768332718;
      else if (t_since_vacc <=  28./12) eff_v_2d = 0.45364139494210;
      else if (t_since_vacc <=  29./12) eff_v_2d = 0.44575282908489;
      else if (t_since_vacc <=  30./12) eff_v_2d = 0.43775036464403;
      else if (t_since_vacc <=  31./12) eff_v_2d = 0.42963235710167;
      else if (t_since_vacc <=  32./12) eff_v_2d = 0.42139713819568;
      else if (t_since_vacc <=  33./12) eff_v_2d = 0.41304301557684;
      else if (t_since_vacc <=  34./12) eff_v_2d = 0.40456827246104;
      else if (t_since_vacc <=  35./12) eff_v_2d = 0.39597116727650;
      else if (t_since_vacc <=  36./12) eff_v_2d = 0.38724993330585;
      else if (t_since_vacc <=  37./12) eff_v_2d = 0.37840277832307;
      else if (t_since_vacc <=  38./12) eff_v_2d = 0.36942788422520;
      else if (t_since_vacc <=  39./12) eff_v_2d = 0.36032340665871;
      else if (t_since_vacc <=  40./12) eff_v_2d = 0.35108747464049;
      else if (t_since_vacc <=  41./12) eff_v_2d = 0.34171819017333;
      else if (t_since_vacc <=  42./12) eff_v_2d = 0.33221362785594;
      else if (t_since_vacc <=  43./12) eff_v_2d = 0.32257183448719;
      else if (t_since_vacc <=  44./12) eff_v_2d = 0.31279082866482;
      else if (t_since_vacc <=  45./12) eff_v_2d = 0.30286860037818;
      else if (t_since_vacc <=  46./12) eff_v_2d = 0.29280311059522;
      else if (t_since_vacc <=  47./12) eff_v_2d = 0.28259229084344;
      else if (t_since_vacc <=  48./12) eff_v_2d = 0.27223404278483;
      else if (t_since_vacc <=  49./12) eff_v_2d = 0.26172623778464;
      else if (t_since_vacc <=  50./12) eff_v_2d = 0.25106671647396;
      else if (t_since_vacc <=  51./12) eff_v_2d = 0.24025328830599;
      else if (t_since_vacc <=  52./12) eff_v_2d = 0.22928373110581;
      else if (t_since_vacc <=  53./12) eff_v_2d = 0.21815579061378;
      else if (t_since_vacc <=  54./12) eff_v_2d = 0.20686718002227;
      else if (t_since_vacc <=  55./12) eff_v_2d = 0.19541557950571;
      else if (t_since_vacc <=  56./12) eff_v_2d = 0.18379863574388;
      else if (t_since_vacc <=  57./12) eff_v_2d = 0.17201396143827;
      else if (t_since_vacc <=  58./12) eff_v_2d = 0.16005913482151;
      else if (t_since_vacc <=  59./12) eff_v_2d = 0.14793169915969;
      else if (t_since_vacc <=  60./12) eff_v_2d = 0.13562916224751;
      else if (t_since_vacc <=  61./12) eff_v_2d = 0.12314899589607;
    break;
  }
  if (t_since_vacc >  61./12)
      eff_v_2d =  0.0;
  return eff_v_2d * (1-(1-0.4688)*0.11) ; /* 11 % of U5 person have VE as 0.4688*VE_adult */
 } double eff_v_1d(double t_since_vacc, int scenario) {
  if (t_since_vacc < 1)
                  return eff_v_2d(t_since_vacc, scenario);
                  else
                  return 0;
                  };
  "

zeronameUnit = paste0(c("C"), 1:10)

pt <- pomp::parameter_trans(
  log = paste0(rep(c(
    "mu_B", "thetaI", "lambdaR", "r", "std_W", "k",
    "betaB", "foi_add"
  ), each = 10), 1:10),
  logit = paste0(rep(c(
    "XthetaA",
    "epsilon",
    "cas_def", "Binit"
  ), each = 10), 1:10)
)

# These params are constant for all departements, so use regex to set all values at once.
all_unit_params[grepl("^sigma[[:digit:]]{1,2}$", names(all_unit_params))] <- 0.25
all_unit_params[grepl("^gamma[[:digit:]]{1,2}$", names(all_unit_params))] <- 182.625
all_unit_params[grepl("^rho[[:digit:]]{1,2}$", names(all_unit_params))] <- 1 / (365 * 8) * 365.25
all_unit_params[grepl("^cases_ext[[:digit:]]{1,2}$", names(all_unit_params))] <- 1

all_unit_params[grepl("^mu[[:digit:]]{1,2}$", names(all_unit_params))] <- 0.01586625546
all_unit_params[grepl("^alpha[[:digit:]]{1,2}$", names(all_unit_params))] <- 1.461
all_unit_params[grepl("^mu_B[[:digit:]]{1,2}$", names(all_unit_params))] <- 133.19716102404308
all_unit_params[grepl("^XthetaA[[:digit:]]{1,2}$", names(all_unit_params))] <- 0.0436160721505241
all_unit_params[grepl("^thetaI[[:digit:]]{1,2}$", names(all_unit_params))] <- 3.4476623459780395e-4
all_unit_params[grepl("^lambdaR[[:digit:]]{1,2}$", names(all_unit_params))] <- 0.2774237712085347
all_unit_params[grepl("^r[[:digit:]]{1,2}$", names(all_unit_params))] <- 0.31360358752214235
all_unit_params[grepl("^std_W[[:digit:]]{1,2}$", names(all_unit_params))] <- 0.008172280355938182
all_unit_params[grepl("^epsilon[[:digit:]]{1,2}$", names(all_unit_params))] <- 0.9750270707877388
all_unit_params[grepl("^k[[:digit:]]{1,2}$", names(all_unit_params))] <- 101.2215999283583
all_unit_params[grepl("^cas_def[[:digit:]]{1,2}$", names(all_unit_params))] <- 0.10
all_unit_params[grepl("^Binit[[:digit:]]{1,2}$", names(all_unit_params))] <- 0.24

# These parameters are different for each departement, so these need to be set seperately.

old_params <- c()

old_params["betaBArtibonite"] =   0.516191
old_params["betaBSud_Est"] =      1.384372
old_params["betaBNippes"] =       2.999928
old_params["betaBNord_Est"] =     3.248645
old_params["betaBOuest"] =        0.090937
old_params["betaBCentre"] =       1.977686
old_params["betaBNord"] =         0.589541
old_params["betaBSud"] =          1.305966
old_params["betaBNord_Ouest"] =   1.141691
old_params["betaBGrande_Anse"] =  2.823539
old_params["foi_addArtibonite"] =    1.530994e-06
old_params["foi_addSud_Est"] =     6.105491e-07
old_params["foi_addNippes"] =       3.056857e-07
old_params["foi_addNord_Est"] =     8.209611e-07
old_params["foi_addOuest"] =       1.070717e-06
old_params["foi_addCentre"] =      0.0000106504579266415
old_params["foi_addNord"] =         5.319736e-07
old_params["foi_addSud"] =          1.030357e-06
old_params["foi_addNord_Ouest"] =   5.855759e-07
old_params["foi_addGrande_Anse"] =  8.762740e-07

for (i in 1:10) {
  dp <- departements[i]
  all_unit_params[paste0("betaB", i)] <- old_params[paste0('betaB', dp)]
  all_unit_params[paste0("foi_add", i)] <- old_params[paste0('foi_add', dp)]
}


# all_rain <- all_rain %>% dplyr::select(time, dplyr::starts_with('rain_std'))


sirb_cholera <- spatPomp::spatPomp(
  data = as.data.frame(all_cases),
  units = "departement",
  times = "time",
  t0 = t_start - dt_years,
  unit_statenames = unit_state_names,
  covar = as.data.frame(all_rain),
  rprocess = euler(step.fun = final_rproc_c, delta.t = dt_years),
  unit_accumvars = c("C"),
  paramnames = names(all_unit_params),
  partrans = pt,
  # global C definitions
  globals = Csnippet(
    stringr::str_c(
      sprintf("int n_cases_start = %i;",  nrow(cases_at_t_start)),
      sprintf("double t_start = %f;",  t_start),
      sprintf("double t_end = %f;",  t_end),
      derivativeBacteria.c,
      all_cases_at_t_start.string,
      eff_v.c,
      sep = " "
    )
  ),
  rinit = initalizeStates,
  dmeasure = dmeas,
  dunit_measure = unit_dmeas,
  rmeasure = rmeas,
  params = all_unit_params
)

coef(sirb_cholera) <- all_unit_params

return(sirb_cholera)
}

