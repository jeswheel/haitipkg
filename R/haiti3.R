#' Build pomp object for model 3.
#'
#' Generate a \sQuote{pomp} object for fitting to Haiti cholera data.
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
#'    This model was developed by Lemaitre, Joseph, et. al at the Laboratory
#'    of Ecohydrology, Ecole Polytechnique Federale de Lausanne (CH).
#'
#'
#'
#' @importFrom magrittr %>%
#' @importFrom foreach %do%
#' @importFrom foreach foreach
#' @seealso \code{\link{haiti2}} and \code{\link{haiti1}} for other models used
#' to fit the cholera epidemic in Haiti.
#' @return \code{\link[pomp]{pomp}} representation of model 3 described in \href{https://www.sciencedirect.com/science/article/pii/S2214109X20303107}{Lee, Elizabeth et. al.} and it's accompanying \href{https://ars.els-cdn.com/content/image/1-s2.0-S2214109X20303107-mmc3.pdf}{Supplemental Material}.
#'
#' @examples
#' mod3 <- haiti3()
#' @export


haiti3 <- function() {

  # Create vector of departement names
  departements = c(
    'Artibonite', 'Centre', 'Grande_Anse',
    'Nippes', 'Nord', 'Nord-Est', 'Nord-Ouest',
    'Ouest', 'Sud', 'Sud-Est'
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
  state_names <- c(
    "S", "I", "A", "RI1", "RI2", "RI3", "RA1", "RA2", "RA3",
    "VSd", "VRI1d", "VRI2d", "VRI3d", "VRA1d", "VRA2d", "VRA3d",
    "VSdd", "VRI1dd", "VRI2dd", "VRI3dd", "VRA1dd", "VRA2dd",
    "VRA3dd", "VSd_alt", "VRI1d_alt", "VRI2d_alt", "VRI3d_alt",
    "VRA1d_alt", "VRA2d_alt", "VRA3d_alt", "VSdd_alt", "VRI1dd_alt",
    "VRI2dd_alt", "VRI3dd_alt", "VRA1dd_alt", "VRA2dd_alt",
    "VRA3dd_alt", "B", "C", "W"
  )

  # All parameters that are common to each departement (disease specific parameters)
  params_common <- c(
    "sigma", "mu_B", "thetaI", "XthetaA", "lambdaR", "r",
    "gammaI", "gammaA", "rhoA", "XrhoI", "epsilon", "k",
    "std_W", "cas_def", "Rtot_0", "mu", "alpha", "cases_ext"
  )

  # Parameters that are unique to each department:
  params_diff <- c(
    "foi_add", "betaB", "H", "D","t_vacc_start",
    "t_vacc_end", "p1d_reg", "r_v_year", "t_vacc_start_alt",
    "t_vacc_end_alt", "p1d_reg_alt", "r_v_year_alt"
  )

  # Loads the input parameters
  # load('R/sysdata.rda')  # These are loaded by the package automatically
  t_start <- dateToYears(as.Date(MODEL3_INPUT_PARAMETERS$t_start))
  t_end   <- dateToYears(as.Date(MODEL3_INPUT_PARAMETERS$t_end))

  all_state_names <- c('IncidenceAll', 'DosesAll', 'CasesAll')
  all_param_names <- params_common

  all_matrix_cases_at_t_start.string <- ""
  all_matrix_cases_other.string <- ""

  # haitiCholera is a saved data.frame in the package
  MODEL3_CASES <- haitiCholera %>%
    dplyr::rename(
      date = date_saturday, Grande_Anse = Grand.Anse,
      `Nord-Est` = Nord.Est, `Nord-Ouest` = Nord.Ouest,
      `Sud-Est` = Sud.Est
    ) %>%
    dplyr::mutate(date = as.Date(date)) %>%
    dplyr::select(-report)

  all_cases <- MODEL3_CASES %>%
    dplyr::mutate(
      date = as.Date(date, format = '%Y-%m-%d'),
      time = dateToYears(date)
    )

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

  # Loop through each of the departements and:
  # (1) Create a dataset with each departement case count
  for (dp in departements) {

    cases <- MODEL3_CASES %>%
      tidyr::gather(dep, cases, -date) %>%
      dplyr::filter(dep == dp) %>%
      dplyr::mutate(
        date = as.Date(date, format = "%Y-%m-%d"),
        time = dateToYears(date)
      )

    all_cases <- cbind(all_cases, placeholder = cases$cases)
    names(all_cases)[names(all_cases) == "placeholder"] <- paste0('cases', gsub('-','_',dp))

    cases_other_dept <- MODEL3_CASES  %>%
      tidyr::gather(dep, cases,-date) %>%
      dplyr::filter(dep != dp) %>%
      dplyr::mutate(
        date = as.Date(date, format = "%Y-%m-%d"),
        time = dateToYears(date)
      )

    cases_other_dept <- aggregate(
      cases_other_dept$cases,
      by = list(Category = cases_other_dept$time),
      FUN = sum,
      na.rm = TRUE,
      na.action = NULL
    ) %>%
      dplyr::mutate(time = Category) %>%
      dplyr::mutate(cases = x)


    cases_at_t_start <- cases %>% dplyr::filter(dateToYears(date) <= t_start)
    cases_at_t_start.string <- foreach::foreach(r = iterators::iter(cases_at_t_start, by = "row"),
                                                .combine = c) %do% {
                                                  sprintf(" {%f, %f} ", r$time, r$cases)
                                                } %>%
      stringr::str_c(collapse = ", \n")

    matrix_cases_at_t_start.string <- stringr::str_c(
      sprintf(
        "double cases_at_t_start%s[%i][%i] = {\n",
        gsub('-', '_', dp),
        nrow(cases_at_t_start),
        2
      ),
      cases_at_t_start.string,
      " \n };"
    )

    # Cases from other departement as a mobility rational
    cases_other.string <-
      foreach::foreach(r = iterators::iter(cases_other_dept, by = "row"),
                       .combine = c) %do% {
                         sprintf(" {%f, %f} ", r$time, r$cases)
                       } %>%
      stringr::str_c(collapse = ", \n")

    matrix_cases_other.string <-
      stringr::str_c(sprintf("double cases_other%s[%i][%i] = {\n", gsub('-', '_', dp), nrow(cases_other_dept), 2),
                     cases_other.string,
                     " \n };")

    all_matrix_cases_at_t_start.string <- stringr::str_c(all_matrix_cases_at_t_start.string, matrix_cases_at_t_start.string)
    all_matrix_cases_other.string <- stringr::str_c(all_matrix_cases_other.string, matrix_cases_other.string)

    all_state_names <- append(all_state_names, lapply(state_names, paste0, gsub('-', '_',dp)))
    all_param_names <- append(all_param_names, lapply(params_diff, paste0, gsub('-', '_',dp)))

    all_state_names <- unlist(all_state_names)
    all_param_names <- unlist(all_param_names)

    all_params <- purrr::set_names(seq_along(all_param_names) * 0, all_param_names)

  }

  for (dp in departements) {
    populations <- unlist(purrr::flatten(MODEL3_INPUT_PARAMETERS["population"]))
    densities <- unlist(purrr::flatten(MODEL3_INPUT_PARAMETERS["density"]))
    all_params[paste0('H', gsub('-', '_', dp))] <- populations[dp]
    all_params[paste0('D', gsub('-', '_', dp))] <- densities[dp]
    p1d_alt_year <- unlist(purrr::flatten(MODEL3_INPUT_PARAMETERS['p1d_alt_year']))
    nb_doses_alt_year <- unlist(purrr::flatten(MODEL3_INPUT_PARAMETERS['nb_doses_alt_year']))
    t_vacc_start_alt  <- unlist(purrr::flatten(MODEL3_INPUT_PARAMETERS["t_vacc_start_alt"]))
    t_vacc_end_alt <- unlist(purrr::flatten(MODEL3_INPUT_PARAMETERS["t_vacc_end_alt"]))
    # Vaccination information:
    t_vacc_start_alt = dateToYears(as.Date(t_vacc_start_alt[dp]))
    t_vacc_end_alt   = dateToYears(as.Date(t_vacc_end_alt[dp]))
    r_v_alt_year = nb_doses_alt_year[dp] / (t_vacc_end_alt - t_vacc_start_alt)
    p1d_alt = p1d_alt_year[dp]

    all_params[paste0("t_vacc_start_alt", gsub('-','_',dp))] = t_vacc_start_alt
    all_params[paste0("t_vacc_end_alt",gsub('-','_',dp))] = t_vacc_end_alt
    all_params[paste0("p1d_reg_alt", gsub('-','_',dp))] = p1d_alt
    all_params[paste0("r_v_year_alt", gsub('-','_',dp))] = r_v_alt_year
  }

  initalizeStatesTemplate =   "
A%s     = nearbyint((1-sigma)/sigma  * 1/epsilon * cases_at_t_start%s[n_cases_start-1][1]/7 * 365 /(mu+gammaA));
I%s     = nearbyint(1/epsilon * cases_at_t_start%s[n_cases_start-1][1]/7 * 365 /(mu+alpha+gammaI))  ;  // Steady state, DP says its correct.

R0[0] = 0;
R0[1] = 0;
B_acc = 0;

for(int i = 0; i < n_cases_start; i++){
R0[0] +=                   cases_at_t_start%s[i][1]/epsilon  * exp((cases_at_t_start%s[i][0] - t_start)  * (rhoI+mu)); /* because t_i in past so t_ - t_0 negative */
R0[1] += (1-sigma)/sigma * cases_at_t_start%s[i][1]/epsilon  * exp((cases_at_t_start%s[i][0] - t_start)  * (rhoA+mu));
B_acc += (thetaA * (1-sigma)/sigma * cases_at_t_start%s[i][1]/epsilon + thetaI * cases_at_t_start%s[i][1]/epsilon) *
(1 + lambdaR * pow(0.024, r)) * D%s * exp((cases_at_t_start%s[i][0] - t_start)  * mu_B);

}

B%s = B_acc;
RI1%s   = nearbyint(R0[0]/3);
RI2%s   = nearbyint(R0[0]/3);
RI3%s   = nearbyint(R0[0]/3);
RA1%s   = nearbyint(R0[1]/3);
RA2%s   = nearbyint(R0[1]/3);
RA3%s   = nearbyint(R0[1]/3);
if (A%s + I%s + RI1%s + RI2%s + RI3%s + RA1%s + RA2%s + RA3%s >= H%s)
{
  double R_tot = H%s - A%s - I%s - 100.0;
  if (R_tot <= 0)
  {
  I%s     = nearbyint(H%s - 100);
  A%s     = nearbyint(0);
  R_tot = nearbyint(0);
  }
  RI1%s   = nearbyint(sigma * R_tot/3.0);
  RI2%s   = nearbyint(sigma * R_tot/3.0);
  RI3%s   = nearbyint(sigma * R_tot/3.0);
  RA1%s   = nearbyint((1-sigma) * R_tot/3.0);
  RA2%s   = nearbyint((1-sigma) * R_tot/3.0);
  RA3%s   = nearbyint((1-sigma) * R_tot/3.0);
}
S%s   = nearbyint(H%s - A%s - I%s - RI1%s - RI2%s - RI3%s - RA1%s - RA2%s - RA3%s);
B%s   = (I%s * thetaI/mu_B + A%s * thetaA/mu_B) * D%s * (1 + lambdaR * pow(0.024, r)); // TODO custom initial conditions equivalent to the 'forcing' in the continous model
C%s   = 0;
W%s   = 0;
VSd%s = 0;
VRI1d%s = 0;
VRI2d%s = 0;
VRI3d%s = 0;
VRA1d%s = 0;
VRA2d%s = 0;
VRA3d%s = 0;
VSdd%s = 0;
VRI1dd%s = 0;
VRI2dd%s = 0;
VRI3dd%s = 0;
VRA1dd%s = 0;
VRA2dd%s = 0;
VRA3dd%s = 0;
VSd_alt%s = 0;
VRI1d_alt%s = 0;
VRI2d_alt%s = 0;
VRI3d_alt%s = 0;
VRA1d_alt%s = 0;
VRA2d_alt%s = 0;
VRA3d_alt%s = 0;
VSdd_alt%s = 0;
VRI1dd_alt%s = 0;
VRI2dd_alt%s = 0;
VRI3dd_alt%s = 0;
VRA1dd_alt%s = 0;
VRA2dd_alt%s = 0;
VRA3dd_alt%s= 0;


"

  initalizeStatesAll <- "double R0[2] = {0,0};
IncidenceAll = 0;
DosesAll = 0;
CasesAll = 0;
double B_acc = 0;
double rhoI = rhoA * XrhoI;
double thetaA = thetaI * XthetaA;"

  for (dp in departements) {
    initalizeStatesAll = paste0(initalizeStatesAll, gsub('%s', gsub('-', '_', dp), initalizeStatesTemplate))
  }

  initalizeStates <- pomp::Csnippet(initalizeStatesAll)


# Build pomp object -------------------------------------------------------??
  rmeasTemplate <- "
  double mean_cases%s = epsilon * C%s;
  if (t > 2018)
  mean_cases%s = mean_cases%s * cas_def;
  cases%s = rnbinom_mu(k, mean_cases%s);
  "

  rmeasAll = ""
  for (dp in departements) {
    rmeasAll = paste0(rmeasAll, gsub('%s', gsub('-', '_', dp), rmeasTemplate))
  }

  rmeas <- pomp::Csnippet(rmeasAll)


## NegBinomial density (if k -> inf then becomes Poisson)
  dmeasTemplate <- "
    double mean_cases%s = epsilon * C%s;
    if (t > 2018)
       mean_cases%s = mean_cases%s * cas_def;
    if (ISNA(cases%s)) {
       lik += (give_log) ? 0 : 1;
    } else {
       if (S%s < 10000) {
          lik += (give_log) ? -99999 : 1.0e-18;
       } else {
           lik += dnbinom_mu(cases%s, k, mean_cases%s, give_log) ;
       }
    }
"

  dmeasAll = "lik = 0;"
  for (dp in departements) {
    dmeasAll = paste0(dmeasAll, gsub('%s', gsub('-', '_', dp), dmeasTemplate))
  }

  dmeas <- pomp::Csnippet(dmeasAll)

  # sirb_file <- 'https://raw.githubusercontent.com/jcblemai/haiti-mass-ocv-campaign/master/scripts/sirb_model_all_dept.c'
  # sirb_file_init <- 'https://raw.githubusercontent.com/jcblemai/haiti-mass-ocv-campaign/master/scripts/sirb_model_all_dept_init.c'

  sirb_file <- "
  previous_vacc_campaign = TRUE ;
r_v_wdn = 0.0;
p1d = 0;
mobility = 0;
pdd = 0;
t_eff =  0;
t_eff_alt = 0;
dw = 0;





/* Compute mobility term */
/*
if (t < t_end) {
      for(int i = 0; i < n_cases_start - 1; i++){
           if (t >= cases_other%s[i][0] && t <= cases_other%s[i+1][0])
               mobility = cases_other%s[i][1];
      }
      if (t > cases_other%s[n_cases_start-1][0])
          mobility = cases_other%s[n_cases_start-1][1];
}*/
//else {
    mobility =  (IArtibonite + ICentre + IGrande_Anse + INippes + INord + INord_Est + INord_Ouest + IOuest + ISud + ISud_Est - I%s) * epsilon * 3.5;
    /* in week for comparison (3.5 because theta is 2) */
    //if (t > 2018)
        //mobility = mobility * cas_def;
//}

  // force of infection
foi = betaB%s * (B%s / (1 + B%s)) + foi_add%s * mobility;
if(std_W > 0.0)
{
    dw = rgammawn(std_W, dt);   // white noise (extra-demographic stochasticity)
    foi_stoc = foi * dw/dt;      // apply stochasticity
} else
{
    foi_stoc = foi;
}

if (t <= (t_vacc_end_alt%s + dt)){
	previous_vacc_campaign = TRUE;
	if (t >= t_vacc_start_alt%s && t <= (t_vacc_end_alt%s + dt)) {
    	r_v_wdn = (r_v_year_alt%s / (S%s + A%s + RI1%s + RI2%s + RI3%s + RA1%s + RA2%s + RA3%s));
	}
	p1d = p1d_reg_alt%s;
} else {
	previous_vacc_campaign = FALSE;
	if (t >= t_vacc_start%s && t <= (t_vacc_end%s + dt)) {
    	r_v_wdn = (r_v_year%s / (S%s + A%s + RI1%s + RI2%s + RI3%s + RA1%s + RA2%s + RA3%s));
	}
	p1d = p1d_reg%s;
}
pdd = 1 - p1d;



// time in the vacc_eff referential. We assume different timing for 1d and 2d
t_eff =     t - (t_vacc_start%s + (t_vacc_end%s - t_vacc_start%s)/2);
t_eff_alt = t - (t_vacc_start_alt%s + (t_vacc_end_alt%s - t_vacc_start_alt%s)/2);


// define transition rates for each type of event (i.e what multplies the thing)
// S compartment
rate[0] = sigma * foi_stoc;   // infections
rate[1] = (1 - sigma) * foi_stoc;   // asymptomatic infections
rate[2] = p1d * r_v_wdn;
rate[3] = pdd * r_v_wdn;
// I compartment
rate[4] = mu;           // natural deaths
rate[5] = alpha;        // cholera-induced deaths
rate[6] = gammaI;       // recovery from infection
// A compartment
rate[7] = mu;           // natural death
rate[8] = gammaA;       // symptoms development
rate[9] = p1d * r_v_wdn;
rate[10] = pdd * r_v_wdn;
// RI1,2,3 compartment
rate[11] = 3*rhoI;        // loss of natural immunity
rate[12] = mu;            // natural death
rate[13] = p1d * r_v_wdn;
rate[14] = pdd * r_v_wdn;
// RA1,2,3 compartment
rate[15] = 3*rhoA;        // loss of natural immunity
rate[16] = mu;            // natural death
rate[17] = p1d * r_v_wdn;
rate[18] = pdd * r_v_wdn;
// V1d_S compartments
rate[19] = sigma       * (1 - eff_v_1d(t_eff, scenario)) * foi_stoc; // symptomatic infections
rate[20] = (1 - sigma) * (1 - eff_v_1d(t_eff, scenario)) * foi_stoc; // asymptomatic infections
rate[21] = mu;          // natural death
// V1d_RI1,2,3 compartment and V2d_RI1,2,3 compartment
rate[22] = mu;          // natural death
rate[23] = 3*rhoI;
// V1d_RA1,2,3 compartment and V2d_RA1,2,3 compartment
rate[24] = mu;          // natural death
rate[25] = 3*rhoA;
// V2d_S compartments
rate[26] = sigma * (1 - eff_v_2d(t_eff, scenario)) * foi_stoc; // symptomatic infections
rate[27] = (1 - sigma) * (1 - eff_v_2d(t_eff, scenario)) * foi_stoc; // asymptomatic infections
rate[28] = mu;          // natural death

/* For previous vacc campagain */
// V1d_S compartments
rate[29] = sigma       * (1 - eff_v_1d(t_eff_alt, scenario)) * foi_stoc; // symptomatic infections
rate[30] = (1 - sigma) * (1 - eff_v_1d(t_eff_alt, scenario)) * foi_stoc; // asymptomatic infections
rate[31] = mu;          // natural death
// V1d_RI1,2,3 compartment and V2d_RI1,2,3 compartment
rate[32] = mu;          // natural death
rate[33] = 3*rhoI;
// V1d_RA1,2,3 compartment and V2d_RA1,2,3 compartment
rate[34] = mu;          // natural death
rate[35] = 3*rhoA;
// V2d_S compartments
rate[36] = sigma *       (1 - eff_v_2d(t_eff_alt, scenario)) * foi_stoc; // symptomatic infections
rate[37] = (1 - sigma) * (1 - eff_v_2d(t_eff_alt, scenario)) * foi_stoc; // asymptomatic infections
rate[38] = mu;          // natural death



// simulate all transitions
/* Probably we can reuse the rates (because const in C function)
but the dN should be different */
reulermultinom(4, S%s,     &rate[0],  dt, &dN[0]);
reulermultinom(3, I%s,     &rate[4],  dt, &dN[4]);
reulermultinom(4, A%s,     &rate[7],  dt, &dN[7]);
reulermultinom(4, RI1%s,   &rate[11], dt, &dN[11]);
reulermultinom(4, RI2%s,   &rate[11], dt, &dN[15]);
reulermultinom(4, RI3%s,   &rate[11], dt, &dN[19]);
reulermultinom(4, RA1%s,   &rate[15], dt, &dN[23]);
reulermultinom(4, RA2%s,   &rate[15], dt, &dN[27]);
reulermultinom(4, RA3%s,   &rate[15], dt, &dN[31]);
/* Vaccinated 1 dose */
reulermultinom(3, VSd%s,   &rate[19], dt, &dN[35]);
reulermultinom(2, VRI1d%s, &rate[22], dt, &dN[38]);
reulermultinom(2, VRI2d%s, &rate[22], dt, &dN[40]);
reulermultinom(2, VRI3d%s, &rate[22], dt, &dN[42]);
reulermultinom(2, VRA1d%s, &rate[24], dt, &dN[44]);
reulermultinom(2, VRA2d%s, &rate[24], dt, &dN[46]);
reulermultinom(2, VRA3d%s, &rate[24], dt, &dN[48]);
/* Vaccinated 2 doses */
reulermultinom(3, VSdd%s,  &rate[26], dt, &dN[50]);
reulermultinom(2, VRI1dd%s,&rate[22], dt, &dN[53]);
reulermultinom(2, VRI2dd%s,&rate[22], dt, &dN[55]);
reulermultinom(2, VRI3dd%s,&rate[22], dt, &dN[57]);
reulermultinom(2, VRA1dd%s,&rate[24], dt, &dN[59]);
reulermultinom(2, VRA2dd%s,&rate[24], dt, &dN[61]);
reulermultinom(2, VRA3dd%s,&rate[24], dt, &dN[63]);
/* For the previous vaccination campain */
/* Vaccinated 1 dose */
reulermultinom(3, VSd_alt%s,   &rate[29], dt, &dN[65]);
reulermultinom(2, VRI1d_alt%s, &rate[32], dt, &dN[68]);
reulermultinom(2, VRI2d_alt%s, &rate[32], dt, &dN[70]);
reulermultinom(2, VRI3d_alt%s, &rate[32], dt, &dN[72]);
reulermultinom(2, VRA1d_alt%s, &rate[34], dt, &dN[74]);
reulermultinom(2, VRA2d_alt%s, &rate[34], dt, &dN[76]);
reulermultinom(2, VRA3d_alt%s, &rate[34], dt, &dN[78]);
/* Vaccinated 2 doses */
reulermultinom(3, VSdd_alt%s,  &rate[36], dt, &dN[80]);
reulermultinom(2, VRI1dd_alt%s,&rate[32], dt, &dN[83]);
reulermultinom(2, VRI2dd_alt%s,&rate[32], dt, &dN[85]);
reulermultinom(2, VRI3dd_alt%s,&rate[32], dt, &dN[87]);
reulermultinom(2, VRA1dd_alt%s,&rate[34], dt, &dN[89]);
reulermultinom(2, VRA2dd_alt%s,&rate[34], dt, &dN[91]);
reulermultinom(2, VRA3dd_alt%s,&rate[34], dt, &dN[93]);

// bacteria as continous state variable
// implement Runge-Kutta integration assuming S, I, R, V* stay constant during dt
k1 = dt * fB(I%s, A%s, B%s, mu_B, thetaI, thetaA, lambdaR, rain_std%s, r, D%s);
k2 = dt * fB(I%s, A%s, B%s, mu_B, thetaI, thetaA, lambdaR, rain_std%s, r, D%s);
k3 = dt * fB(I%s, A%s, B%s, mu_B, thetaI, thetaA, lambdaR, rain_std%s, r, D%s);
k4 = dt * fB(I%s, A%s, B%s, mu_B, thetaI, thetaA, lambdaR, rain_std%s, r, D%s);
// bacteria increment
dB = (k1 + 2*k2 + 2*k3 + k4) / 6.0;


I%s   += dN[0] + dN[35] + dN[50] + dN[35+30] + dN[50+30] - dN[4] - dN[5] - dN[6];
A%s   += dN[1] + dN[36] + dN[51] + dN[36+30] + dN[51+30] - dN[7] - dN[8] - dN[9] - dN[10];
RI1%s += dN[6] -  dN[11] - dN[12] - dN[13] - dN[14] ;
RI2%s += dN[11] - dN[15] - dN[16] - dN[17] - dN[18];
RI3%s += dN[15] - dN[19] - dN[20] - dN[21] - dN[22];
RA1%s += dN[8]  - dN[23] - dN[24] - dN[25] - dN[26];
RA2%s += dN[23] - dN[27] - dN[28] - dN[29] - dN[30];
RA3%s += dN[27] - dN[31] - dN[32] - dN[33] - dN[34];

if (previous_vacc_campaign){
	VSd_alt%s    += dN[2];
	VRI1d_alt%s  += dN[13];
	VRI2d_alt%s  += dN[17];
	VRI3d_alt%s  += dN[21];
	VRA1d_alt%s  += dN[9] + dN[25] ;
	VRA2d_alt%s  += dN[29];
	VRA3d_alt%s  += dN[33];

	VSdd_alt%s    += dN[3];
	VRI1dd_alt%s  += dN[14];
	VRI2dd_alt%s  += dN[18];
	VRI3dd_alt%s  += dN[22];
	VRA1dd_alt%s  += dN[10] + dN[26];
	VRA2dd_alt%s  += dN[30];
	VRA3dd_alt%s  += dN[34];
} else {
	VSd%s    += dN[2];
	VRI1d%s  += dN[13];
	VRI2d%s  += dN[17];
	VRI3d%s  += dN[21];
	VRA1d%s  += dN[9] + dN[25] ;
	VRA2d%s  += dN[29];
	VRA3d%s  += dN[33];

	VSdd%s    += dN[3];
	VRI1dd%s  += dN[14];
	VRI2dd%s  += dN[18];
	VRI3dd%s  += dN[22];
	VRA1dd%s  += dN[10] + dN[26];
	VRA2dd%s  += dN[30];
	VRA3dd%s  += dN[34];

}

VSd%s    += dN[43] + dN[49] - dN[35] - dN[36] - dN[37];
VRI1d%s  += - dN[38] - dN[39];
VRI2d%s  += dN[39] - dN[40] - dN[41];
VRI3d%s  += dN[41] - dN[42] - dN[43];
VRA1d%s  += - dN[44] - dN[45];
VRA2d%s  += dN[45] - dN[46] - dN[47];
VRA3d%s  += dN[47] - dN[48] - dN[49];

VSdd%s    +=  dN[58] + dN[64] - dN[50] - dN[51] - dN[52];
VRI1dd%s  += - dN[53] - dN[54];
VRI2dd%s  +=  dN[54] - dN[55] - dN[56] ;
VRI3dd%s  +=  dN[56] - dN[57] - dN[58];
VRA1dd%s  += - dN[59] - dN[60];
VRA2dd%s  +=  dN[60] - dN[61] - dN[62];
VRA3dd%s  +=  dN[62] - dN[63] - dN[64];

/* *previous* vacccination campain */

VSd_alt%s    += dN[43+30] + dN[49+30] - dN[35+30] - dN[36+30] - dN[37+30];
VRI1d_alt%s  += - dN[38+30] - dN[39+30];
VRI2d_alt%s  += dN[39+30] - dN[40+30] - dN[41+30];
VRI3d_alt%s  += dN[41+30] - dN[42+30] - dN[43+30];
VRA1d_alt%s  += - dN[44+30] - dN[45+30];
VRA2d_alt%s  += dN[45+30] - dN[46+30] - dN[47+30];
VRA3d_alt%s  += dN[47+30] - dN[48+30] - dN[49+30];

VSdd_alt%s    +=  dN[58+30] + dN[64+30] - dN[50+30] - dN[51+30] - dN[52+30];
VRI1dd_alt%s  += - dN[53+30] - dN[54+30];
VRI2dd_alt%s  +=  dN[54+30] - dN[55+30] - dN[56+30] ;
VRI3dd_alt%s  +=  dN[56+30] - dN[57+30] - dN[58+30];
VRA1dd_alt%s  += - dN[59+30] - dN[60+30];
VRA2dd_alt%s  +=  dN[60+30] - dN[61+30] - dN[62+30];
VRA3dd_alt%s  +=  dN[62+30] - dN[63+30] - dN[64+30];


C%s   +=  dN[0] + dN[35] + dN[50] + dN[35+30] + dN[50+30];
W%s   +=  (dw - dt)/std_W;  // standardized i.i.d. white noise
B%s += (((dB) < -B%s) ? (-B%s + 1.0e-3) : (dB)); // condition to ensure B>0

// susceptibles so as to match total population
S%s = nearbyint(H%s - I%s - A%s - RI1%s - RI2%s - RI3%s - RA1%s - RA2%s - RA3%s -
	VSd%s - VRI1d%s - VRI2d%s - VRI3d%s - VRA1d%s - VRA2d%s -VRA3d%s -
	VSdd%s- VRI1dd%s -VRI2dd%s -VRI3dd%s -VRA1dd%s-VRA2dd%s-VRA3dd%s -
	VSd_alt%s - VRI1d_alt%s - VRI2d_alt%s - VRI3d_alt%s - VRA1d_alt%s - VRA2d_alt%s - VRA3d_alt%s -
	VSdd_alt%s - VRI1dd_alt%s - VRI2dd_alt%s - VRI3dd_alt%s - VRA1dd_alt%s- VRA2dd_alt%s - VRA3dd_alt%s);



IncidenceAll +=  dN[0] + dN[35] + dN[50] + dN[35+30] + dN[50+30] + dN[1] + dN[36] + dN[51] + dN[36+30] + dN[51+30];
if (!previous_vacc_campaign)
{
	DosesAll  += dN[2] + dN[13] + dN[17] + dN[21] + dN[9] + dN[25] + dN[29] + dN[33] + 2*(dN[3] + dN[14] + dN[18] + dN[22] + dN[10] + dN[26] + dN[30] + dN[34]);
}
CasesAll  +=  dN[0] + dN[35] + dN[50] + dN[35+30] + dN[50+30];
  "

  sirb.rproc <- "
double foi, foi_stoc; // force of infection and its stochastic version
double dw;            // extra-demographic stochasticity on foi
double dB;            // deterministic forward time difference of bacteria in the environment
double k1, k2, k3, k4;  // coefficients of  the Runge-Kutta method
double rate[39];      // vector of all rates in model
double dN[95];        // vector of transitions between classes during integration timestep

double thetaA = thetaI * XthetaA;
double rhoI = rhoA * XrhoI;

int previous_vacc_campaign = TRUE ; /* flag that indicate if we are on the first or second campain */
double r_v_wdn = 0.0;       // rate of vaccination: 0 if out of time window, r_v if not
double p1d = 0;
double mobility = 0;

int scenario =  cases_ext;

double pdd = 0;
double t_eff =  0;
double t_eff_alt = 0;

  "

  for (dp in departements) {
    sirb.rproc = paste0(sirb.rproc, gsub('%s', gsub('-', '_', dp), sirb_file))
  }

  sirb.rproc = pomp::Csnippet(sirb.rproc)

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

  zeronameTemplate = c("C", "W")
  zeronameAll = c('IncidenceAll', 'CasesAll')
  for (dp in departements){
    zeronameAll = append(zeronameAll, lapply(zeronameTemplate, paste0, gsub('-', '_', dp)))

  }
  zeronameAll <- unlist(zeronameAll)

  dt_yrs <- 1 / 365.25 * .2

  pt <- pomp::parameter_trans(
    log = c(
      "mu_B", "thetaI", "lambdaR", "r", "std_W", "k",
      "betaBArtibonite", "foi_addArtibonite", "betaBCentre",
      "foi_addCentre", "betaBGrande_Anse", "foi_addGrande_Anse",
      "betaBNippes", "foi_addNippes", "betaBNord", "foi_addNord",
      "betaBNord_Est", "foi_addNord_Est", "betaBNord_Ouest",
      "foi_addNord_Ouest", "betaBOuest", "foi_addOuest",
      "betaBSud", "foi_addSud", "betaBSud_Est", "foi_addSud_Est"
    ),
    logit = c(
      "XthetaA",
      "epsilon",
      "cas_def"
    )
  )

  all_params["sigma"] <- .25               # Fixed
  all_params["rhoA"] <- 1 / (365 * 8) * 365.25   # Fixed  \rho in table 14
  all_params["XrhoI"] <- 1                 # Fixed, saying loss of infected for those in compartments A and I are the same
  all_params["gammaA"] <- 182.625          # TODO: where does this number come from? I think this is a mistake by the authors. Fixed  divide by 365 to get \gamma in table 14
  all_params["gammaI"] <- 182.625          # Fixed  divide by 365 to get \gamma in table 14
  all_params["Rtot_0"] <- 0.35             # Useless


  all_params['cases_ext'] <- 1
  all_params['mu'] <-  0.01586625546  # divide by 365 to get \mu in table 14
  all_params['alpha'] <- 1.461  # divide by 365 to get \alpha in table 14

  # From calibration
  all_params["mu_B"] <-  133.19716102404308  # This value mu_B / 365 matches calibrated \mu_B in table
  all_params["XthetaA"] <- 0.0436160721505241  # This multiplies thetaI to get \theta_A in table S15
  all_params["thetaI"] <- 3.4476623459780395e-4  # matches \theta_I in table S15
  all_params["lambdaR"] <- 0.2774237712085347  # matches \lambda in table S15
  all_params["r"] <- 0.31360358752214235  # Matches table S15
  all_params["std_W"] <- 0.008172280355938182  # Matches table S15
  all_params["epsilon"] <- 0.9750270707877388  # Matches table S15
  all_params["k"] <- 101.2215999283583  # Matches p in table S15
  all_params["cas_def"] <- 0.10  # TODO: It is clear this wasn't fit using MIF, but the authors say that it was. Matches \epsilon_2 in table S15

  all_params["betaBArtibonite"] =   0.516191
  all_params["betaBSud_Est"] =      1.384372
  all_params["betaBNippes"] =       2.999928
  all_params["betaBNord_Est"] =     3.248645
  all_params["betaBOuest"] =        0.090937
  all_params["betaBCentre"] =       1.977686
  all_params["betaBNord"] =         0.589541
  all_params["betaBSud"] =          1.305966
  all_params["betaBNord_Ouest"] =   1.141691
  all_params["betaBGrande_Anse"] =  2.823539

  all_params["foi_addArtibonite"] =    1.530994e-06
  all_params["foi_addSud_Est"] =     6.105491e-07
  all_params["foi_addNippes"] =       3.056857e-07
  all_params["foi_addNord_Est"] =     8.209611e-07
  all_params["foi_addOuest"] =       1.070717e-06
  all_params["foi_addCentre"] =      0.0000106504579266415
  all_params["foi_addNord"] =         5.319736e-07
  all_params["foi_addSud"] =          1.030357e-06
  all_params["foi_addNord_Ouest"] =   5.855759e-07
  all_params["foi_addGrande_Anse"] =  8.762740e-07

  all_rain <- all_rain %>% dplyr::select(time, dplyr::starts_with('rain_std'))

  sirb_cholera <- pomp::pomp(
    # set data
    data = all_cases %>%
      dplyr::filter(time > t_start & time < (t_end + 0.01)) %>%
      dplyr::select(
        time, casesArtibonite, casesCentre,
        casesGrande_Anse, casesNippes,
        casesNord, casesNord_Est, casesOuest,
        casesSud, casesSud_Est, casesNord_Ouest
      ),
    # time column
    times = "time",
    # initialization time
    t0 = t_start - dt_yrs,
    # paramter vector
    params = all_params,

    # process simulator
    rprocess = pomp::euler(step.fun = sirb.rproc, delta.t = dt_yrs),
    # measurement model simulator
    rmeasure =  rmeas,
    # measurement model density
    dmeasure = dmeas,
    # covariates
    covar = pomp::covariate_table(all_rain, times = 'time'),
    # names of state variables
    statenames = all_state_names,
    # names of accumulator variables to be re-initalized at each observation timestep
    # (C for cases, W for the white noise just for plotting)
    accumvars = zeronameAll,
    # names of paramters
    paramnames = all_param_names,
    # names of covariates
    rinit = initalizeStates,
    partrans = pt,
    # global C definitions
    globals = stringr::str_c(
      sprintf("int n_cases_start = %i;",  nrow(cases_at_t_start)),
      sprintf("int n_cases_other = %i;",  nrow(cases_other_dept)),
      sprintf("double t_start = %f;",  t_start),
      sprintf("double t_end = %f;",  t_end),
      derivativeBacteria.c,
      all_matrix_cases_at_t_start.string,
      all_matrix_cases_other.string,
      eff_v.c,
      sep = " "
    )
  )

  return(sirb_cholera)
}
