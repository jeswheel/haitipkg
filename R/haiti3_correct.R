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
#' @seealso \code{\link{haiti2}} and \code{\link{haiti1}} for other models used
#' to fit the cholera epidemic in Haiti.
#' @return \code{\link[pomp]{pomp}} representation of model 3 described in \href{https://www.sciencedirect.com/science/article/pii/S2214109X20303107}{Lee, Elizabeth et. al.} and it's accompanying \href{https://ars.els-cdn.com/content/image/1-s2.0-S2214109X20303107-mmc3.pdf}{Supplemental Material}.
#'
#' @examples
#' mod3 <- haiti3_correct()
#' @export

haiti3_correct <- function() {

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
    "S", "I", "A", "R1", "R2", "R3",
    "VSd", "VR1d", "VR2d", "VR3d",
    "VSdd", "VR1dd", "VR2dd", "VR3dd",
    # "VAd", "VAdd", "VAd_alt", "VAdd_alt",
    "VSd_alt", "VR1d_alt", "VR2d_alt", "VR3d_alt",
    "VSdd_alt", "VR1dd_alt",
    "VR2dd_alt", "VR3dd_alt",
    "B", "C", "W"
  )

  # All parameters that are common to each departement (disease specific parameters)
  params_common <- c(
    "sigma", "mu_B", "thetaI", "XthetaA", "lambdaR", "r",
    "gamma", "rho", "epsilon", "k",
    "std_W", "cas_def", "Rtot_0", "mu", "alpha", "cases_ext"
  )

  # Parameters that are unique to each department:
  params_diff <- c(
    "foi_add", "betaB", "H", "D", "t_vacc_start",
    "t_vacc_end", "p1d_reg", "r_v_year", "t_vacc_start_alt",
    "t_vacc_end_alt", "p1d_reg_alt", "r_v_year_alt"
  )

  # Load the input parameters
  t_start <- dateToYears(as.Date(MODEL3_INPUT_PARAMETERS$t_start))
  t_end   <- dateToYears(as.Date(MODEL3_INPUT_PARAMETERS$t_end))

  all_state_names <- c('IncidenceAll', 'DosesAll', 'CasesAll')
  all_param_names <- params_common

  all_matrix_cases_at_t_start.string <- ""
  all_matrix_cases_other.string <- ""

  # MODEL3_CASES is imported from the internal data: R/sysdata.rda
  all_cases <- MODEL3_CASES %>%
    dplyr::mutate(
      date = as.Date(date, format = '%Y-%m-%d'),
      time = dateToYears(date)
    )

  # MODEL3_RAIN is imported from the internal data: R/sysdata.rda
  all_rain <- MODEL3_RAIN %>%
    dplyr::mutate(
      date = as.Date(date, format = "%Y-%m-%d"),
      time = dateToYears(date)
    ) %>%
    dplyr::filter(time > t_start - 0.01 & time < (t_end + 0.01))

  # Loop through each of the departements and:
  # (1) Create a dataset with each departement case count
  # (2) Create a dataset with each departement rain history
  for (dp in departements) {

    cases <- MODEL3_CASES %>%
      tidyr::gather(dep, cases, -date) %>%
      dplyr::filter(dep == dp) %>%
      dplyr::mutate(
        date = as.Date(date, format = "%Y-%m-%d"),
        time = dateToYears(date)
      )

    rain <- MODEL3_RAIN %>%
      tidyr::gather(dep, rain,-date) %>%
      dplyr::group_by(dep) %>%
      dplyr::ungroup() %>%
      dplyr::filter(dep == dp) %>%
      dplyr::mutate(
        date = as.Date(date, format = "%Y-%m-%d"),
        time = dateToYears(date)
      ) %>%
      dplyr::filter(time > t_start - 0.01 & time < (t_end + 0.01)) %>%
      dplyr::mutate(max_rain = max(rain), rain_std = rain / max_rain)

    all_rain <- cbind(all_rain, placeholder = rain$max_rain)
    all_rain <- cbind(all_rain, placeholder2 = rain$rain_std)
    names(all_rain)[names(all_rain) == "placeholder"] <- paste0('max_rain', gsub('-','_',dp))
    names(all_rain)[names(all_rain) == "placeholder2"] <- paste0('rain_std', gsub('-','_',dp))

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
    populations <- unlist(purrr::flatten(input_parameters["population"]))
    densities <- unlist(purrr::flatten(input_parameters["density"]))
    all_params[paste0('H', gsub('-', '_', dp))] <- populations[dp]
    all_params[paste0('D', gsub('-', '_', dp))] <- densities[dp]
    p1d_alt_year <- unlist(purrr::flatten(input_parameters['p1d_alt_year']))
    nb_doses_alt_year <- unlist(purrr::flatten(input_parameters['nb_doses_alt_year']))
    t_vacc_start_alt  <- unlist(purrr::flatten(input_parameters["t_vacc_start_alt"]))
    t_vacc_end_alt <- unlist(purrr::flatten(input_parameters["t_vacc_end_alt"]))
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
A%s     = nearbyint((1-sigma)/sigma  * 1/epsilon * cases_at_t_start%s[n_cases_start-1][1]/7 * 365 /(mu + gamma));
I%s     = nearbyint(1/epsilon * cases_at_t_start%s[n_cases_start-1][1]/7 * 365 /(mu+alpha + gamma))  ;  // Steady state, DP says its correct.

B_acc = 0;

for(int i = 0; i < n_cases_start; i++){
R0[0] +=                   cases_at_t_start%s[i][1]/epsilon  * exp((cases_at_t_start%s[i][0] - t_start)  * (rho + mu)); /* because t_i in past so t_ - t_0 negative */
R0[1] += (1-sigma)/sigma * cases_at_t_start%s[i][1]/epsilon  * exp((cases_at_t_start%s[i][0] - t_start)  * (rho + mu));
B_acc += (thetaA * (1-sigma)/sigma * cases_at_t_start%s[i][1]/epsilon + thetaI * cases_at_t_start%s[i][1]/epsilon) *
(1 + lambdaR * pow(0.024, r)) * D%s * exp((cases_at_t_start%s[i][0] - t_start)  * mu_B);

}

B%s = B_acc;
R1%s   = nearbyint((R0[0] + R0[1]) / 3);
R2%s   = nearbyint((R0[0] + R0[1]) / 3);
R3%s   = nearbyint((R0[0] + R0[1]) / 3);

if (A%s + I%s + R1%s + R2%s + R3%s >= H%s)
{
  double R_tot = H%s - A%s - I%s - 100.0;
  if (R_tot <= 0)
  {
  I%s     = nearbyint(H%s - 100);
  A%s     = nearbyint(0);
  R_tot = nearbyint(0);
  }
  R1%s   = nearbyint(R_tot / 3);
  R2%s   = nearbyint(R_tot / 3);
  R3%s   = nearbyint(R_tot / 3);
}
S%s   = nearbyint(H%s - A%s - I%s - R1%s - R2%s - R3%s);
B%s   = (I%s * thetaI/mu_B + A%s * thetaA/mu_B) * D%s * (1 + lambdaR * pow(0.024, r)); // TODO this just overwrites what was done above???
C%s   = 0;
W%s   = 0;
VSd%s = 0;
VAd%s = 0;
VR1d%s = 0;
VR2d%s = 0;
VR3d%s = 0;
VSdd%s = 0;
VAdd%s = 0;
VR1dd%s = 0;
VR2dd%s = 0;
VR3dd%s = 0;

VSd_alt%s = 0;
VR1d_alt%s = 0;
VR2d_alt%s = 0;
VR3d_alt%s = 0;
VSdd_alt%s = 0;
VR1dd_alt%s = 0;
VR2dd_alt%s = 0;
VR3dd_alt%s = 0;

"

  initalizeStatesAll <- "double R0[2] = {0,0};
IncidenceAll = 0;
DosesAll = 0;
CasesAll = 0;
double B_acc = 0;
double thetaA = thetaI * XthetaA;"

  for (dp in departements) {
    initalizeStatesAll = paste0(initalizeStatesAll, gsub('%s', gsub('-', '_', dp), initalizeStatesTemplate))
  }

  initalizeStates <- pomp::Csnippet(initalizeStatesAll)

  ###
  ### rmeas
  ###

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

  ###
  ### dmeas
  ###

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

  ###
  ### rproc
  ###

  sirb_file <- "
previous_vacc_campaign = TRUE;
r_v_wdn = 0.0;
p1d = 0;
pdd = 0;
t_eff =  0;
t_eff_alt = 0;
dw = 0;

mobility =  (IArtibonite + ICentre + IGrande_Anse + INippes + INord +
             INord_Est + INord_Ouest + IOuest + ISud + ISud_Est - I%s +
             AArtibonite + ACentre + AGrande_Anse + ANippes + ANord +
             ANord_Est + ANord_Ouest + AOuest + ASud + ASud_Est - A%s);

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
rate[0] = sigma * foi_stoc;         // infections
rate[1] = (1 - sigma) * foi_stoc;   // asymptomatic infections
rate[2] = p1d * r_v_wdn;            // S -> VSd
rate[3] = pdd * r_v_wdn;            // S -> VSdd
rate[4] = mu;                       // S -> natural death

// I compartment
rate[5] = mu;            // natural deaths
rate[6] = alpha;         // cholera-induced deaths
rate[7] = gamma;         // I -> R

// A compartment
rate[8] = mu;               // natural death
rate[9] = gamma;            // A -> R
rate[10] = p1d * r_v_wdn;   // A -> VAd
rate[11] = pdd * r_v_wdn;   // A -> VAdd

// R1
rate[12] = 3 * rho;        // loss of natural immunity; R1 -> R2
rate[13] = mu;             // natural death R1 -> death
rate[14] = p1d * r_v_wdn;  // R1 -> VR1d
rate[15] = pdd * r_v_wdn;  // R1 -> VR1dd

// R2
rate[16] = 3 * rho;        // loss of natural immunity; R2 -> R3
rate[17] = mu;             // natural death R2 -> death
rate[18] = p1d * r_v_wdn;  // R2 -> VR2d
rate[19] = pdd * r_v_wdn;  // R2 -> VR2dd

// R3
rate[20] = 3 * rho;        // loss of natural immunity; R3 -> S
rate[21] = mu;             // natural death R3 -> death
rate[22] = p1d * r_v_wdn;  // R3 -> VR3d
rate[23] = pdd * r_v_wdn;  // R3 -> VR3dd

// VSd
rate[24] = sigma * (1 - eff_v_1d(t_eff, scenario)) * foi_stoc;       // VSd -> I
rate[25] = (1 - sigma) * (1 - eff_v_1d(t_eff, scenario)) * foi_stoc; // VSd -> A
rate[26] = mu;                                                       // VSd -> death

// VAd
rate[27] = gamma;       // VAd -> VR1d
rate[28] = mu;          // VAd -> natural death

// VR1d
rate[29] = mu;          // natural death:    VR1d -> death
rate[30] = 3 * rho;     // loss of immunity: VR1d -> VR2d

// VR2d
rate[31] = mu;          // natural death:    VR2d -> death
rate[32] = 3 * rho;     // loss of immunity: VR2d -> VR3d

// VR3d
rate[33] = mu;          // natural death:    VR3d -> death
rate[34] = 3 * rho;     // loss of immunity: VR3d -> VSd

// VSdd
rate[35] = sigma * (1 - eff_v_2d(t_eff, scenario)) * foi_stoc;       // VSdd -> I
rate[36] = (1 - sigma) * (1 - eff_v_2d(t_eff, scenario)) * foi_stoc; // VSdd -> A
rate[37] = mu;                                                       // natural death

// VAdd
rate[38] = gamma;       // VAdd -> VR1dd
rate[39] = mu;          // VAdd -> natural death

// VR1dd
rate[40] = mu;          // natural death:    VR1dd -> death
rate[41] = 3 * rho;     // loss of immunity: VR1dd -> VR2dd

// VR2dd
rate[42] = mu;          // natural death:    VR2dd -> death
rate[43] = 3 * rho;     // loss of immunity: VR2dd -> VR3dd

// VR3dd
rate[44] = mu;          // natural death:    VR3dd -> death
rate[45] = 3 * rho;     // loss of immunity: VR3dd -> VSdd


/* For previous vacc campagain */

// VSd_alt
rate[46] = sigma       * (1 - eff_v_1d(t_eff_alt, scenario)) * foi_stoc; // VSd -> I
rate[47] = (1 - sigma) * (1 - eff_v_1d(t_eff_alt, scenario)) * foi_stoc; // VSd -> A
rate[48] = mu;                                                           // natural death

// VSdd_alt
rate[49] = sigma       * (1 - eff_v_2d(t_eff_alt, scenario)) * foi_stoc; // VSdd -> I
rate[50] = (1 - sigma) * (1 - eff_v_2d(t_eff_alt, scenario)) * foi_stoc; // VSdd -> A
rate[51] = mu;                                                           // natural death

// All other alts already defined.

// simulate all transitions
reulermultinom(5, S%s,     &rate[0],  dt, &dN[0]);
reulermultinom(3, I%s,     &rate[5],  dt, &dN[5]);
reulermultinom(4, A%s,     &rate[8],  dt, &dN[8]);
reulermultinom(4, R1%s,    &rate[12], dt, &dN[12]);
reulermultinom(4, R2%s,    &rate[16], dt, &dN[16]);
reulermultinom(4, R3%s,    &rate[20], dt, &dN[20]);

/* Vaccinated 1 dose */
reulermultinom(3,  VSd%s,   &rate[24], dt, &dN[24]);
reulermultinom(2,  VAd%s,   &rate[27], dt, &dN[27]);
reulermultinom(2, VR1d%s,   &rate[29], dt, &dN[29]);
reulermultinom(2, VR2d%s,   &rate[31], dt, &dN[31]);
reulermultinom(2, VR3d%s,   &rate[33], dt, &dN[33]);

/* Vaccinated 2 doses */
reulermultinom(3,  VSdd%s,   &rate[35], dt, &dN[35]);
reulermultinom(2,  VAdd%s,   &rate[38], dt, &dN[38]);
reulermultinom(2, VR1dd%s,   &rate[40], dt, &dN[40]);
reulermultinom(2, VR2dd%s,   &rate[42], dt, &dN[42]);
reulermultinom(2, VR3dd%s,   &rate[44], dt, &dN[44]);

/* For the previous vaccination campain */

/* Alt Vaccinated 1 dose */
reulermultinom(3,  VSd_alt%s,   &rate[46], dt, &dN[46]);
reulermultinom(2, VR1d_alt%s,   &rate[29], dt, &dN[49]);
reulermultinom(2, VR2d_alt%s,   &rate[31], dt, &dN[51]);
reulermultinom(2, VR3d_alt%s,   &rate[33], dt, &dN[53]);

/* Alt Vaccinated 2 doses */
reulermultinom(3,  VSdd_alt%s,   &rate[49], dt, &dN[55]);
reulermultinom(2, VR1dd_alt%s,   &rate[29], dt, &dN[58]);
reulermultinom(2, VR2dd_alt%s,   &rate[31], dt, &dN[60]);
reulermultinom(2, VR3dd_alt%s,   &rate[33], dt, &dN[62]);

// bacteria as continous state variable
// implement Runge-Kutta integration assuming S, I, R, V* stay constant during dt
k1 = dt * fB(I%s, A%s, B%s, mu_B, thetaI, thetaA, lambdaR, rain_std%s, r, D%s);
k2 = dt * fB(I%s, A%s, B%s, mu_B, thetaI, thetaA, lambdaR, rain_std%s, r, D%s);
k3 = dt * fB(I%s, A%s, B%s, mu_B, thetaI, thetaA, lambdaR, rain_std%s, r, D%s);
k4 = dt * fB(I%s, A%s, B%s, mu_B, thetaI, thetaA, lambdaR, rain_std%s, r, D%s);
// bacteria increment
dB = (k1 + 2*k2 + 2*k3 + k4) / 6.0;


I%s += dN[0] + dN[24] + dN[35] + dN[46] + dN[55] - dN[7] - dN[6] - dN[5];      // I += S + VSd + VSdd + VSd_alt + VSdd_alt - R1 - ND - CD
A%s += dN[1] + dN[25] + dN[36] + dN[47] + dN[56] - dN[9] - dN[10] - dN[11] - dN[8];     // A += S + VSd + VSdd + VSd_alt + VSdd_alt - R1 - VR1d - VR1dd - death
R1%s += dN[7] + dN[9] - dN[12] - dN[14] - dN[15] - dN[13];  // I-> R1, A -> R1, R1 -> R2, R1 -> VR1d, R1 -> VR1dd, R1 -> death
R2%s += dN[12] - dN[16] - dN[18] - dN[19] - dN[17]; // R1 -> R2, R2 -> R3, R2 -> VR2d, R2 -> VR2dd, R2 -> death
R3%s += dN[16] - dN[20] - dN[22] - dN[23] - dN[21]; // R2 -> R3, R3 -> S , R3 -> VR3d, R3 -> VR3dd, R3 -> death

//////////////// COMPLETE UNTIL HERE /////////////////

if (previous_vacc_campaign){
	VSd_alt%s    += dN[];
	VR1d_alt%s  += dN[];
	VR2d_alt%s  += dN[];
	VR3d_alt%s  += dN[];
	// VRA1d_alt%s  += dN[] + dN[] ; TODO: WHY two adds?

	VSdd_alt%s    += dN[];
	VR1dd_alt%s  += dN[];
	VR2dd_alt%s  += dN[];
	VR3dd_alt%s  += dN[];

} else {
	VSd%s    += dN[];
	VR1d%s  += dN[];
	VR2d%s  += dN[];
	VR3d%s  += dN[];

	VSdd%s    += dN[];
	VR1dd%s  += dN[];
	VR2dd%s  += dN[];
	VR3dd%s  += dN[];
}

VSd%s    += dN[] + dN[] - dN[] - dN[] - dN[];
VR1d%s  += - dN[] - dN[];
VR2d%s  += dN[] - dN[] - dN[];
VR3d%s  += dN[] - dN[] - dN[];

VSdd%s    +=  dN[] + dN[] - dN[] - dN[] - dN[];
VR1dd%s  += - dN[] - dN[];
VR2dd%s  +=  dN[] - dN[] - dN[] ;
VR3dd%s  +=  dN[] - dN[] - dN[];

/* *previous* vacccination campain */

VSd_alt%s    += dN[] + dN[] - dN[] - dN[] - dN[];
VR1d_alt%s  += - dN[] - dN[];
VR2d_alt%s  += dN[] - dN[] - dN[];
VR3d_alt%s  += dN[] - dN[] - dN[];

VSdd_alt%s    +=  dN[] + dN[] - dN[] - dN[] - dN[];
VR1dd_alt%s  += - dN[] - dN[];
VR2dd_alt%s  +=  dN[] - dN[] - dN[] ;
VR3dd_alt%s  +=  dN[] - dN[] - dN[];

C%s   +=  dN[] + dN[] + dN[] + dN[] + dN[];
W%s   +=  (dw - dt)/std_W;  // standardized i.i.d. white noise
B%s += (((dB) < -B%s) ? (-B%s + 1.0e-3) : (dB)); // condition to ensure B>0

// susceptibles so as to match total population
S%s = nearbyint(H%s - I%s - A%s - RI1%s - RI2%s - RI3%s - RA1%s - RA2%s - RA3%s -
	VSd%s - VRI1d%s - VRI2d%s - VRI3d%s - VRA1d%s - VRA2d%s -VRA3d%s -
	VSdd%s- VRI1dd%s -VRI2dd%s -VRI3dd%s -VRA1dd%s-VRA2dd%s-VRA3dd%s -
	VSd_alt%s - VRI1d_alt%s - VRI2d_alt%s - VRI3d_alt%s - VRA1d_alt%s - VRA2d_alt%s - VRA3d_alt%s -
	VSdd_alt%s - VRI1dd_alt%s - VRI2dd_alt%s - VRI3dd_alt%s - VRA1dd_alt%s- VRA2dd_alt%s - VRA3dd_alt%s);



IncidenceAll +=  dN[] + dN[] + dN[] + dN[] + dN[] + dN[] + dN[] + dN[] + dN[] + dN[];
if (!previous_vacc_campaign)
{
	DosesAll  += dN[] + dN[] + dN[] + dN[] + dN[] + dN[] + dN[] + dN[] + 2*(dN[] + dN[] + dN[] + dN[] + dN[] + dN[] + dN[] + dN[]);
}
CasesAll  +=  dN[] + dN[] + dN[] + dN[] + dN[];
  "

}
