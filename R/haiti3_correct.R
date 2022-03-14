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

haiti3_correct <- function(dt_yrs = 1 / 365.25 * .5) {

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
    "S", "I", "A", "R1", "R2", "R3",
    "VSd", "VR1d", "VR2d", "VR3d",
    "VSdd", "VR1dd", "VR2dd", "VR3dd",
    # "VAd", "VAdd", "VAd_alt", "VAdd_alt",
    "VSd_alt", "VR1d_alt", "VR2d_alt", "VR3d_alt",
    "VSdd_alt", "VR1dd_alt",
    "VR2dd_alt", "VR3dd_alt",
    "B", "C", "W"
    # "MobAI", "MobCases" # TODO: Remove this when done.
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
    "t_vacc_end_alt", "p1d_reg_alt", "r_v_year_alt", "B0"
  )

  # Load the input parameters
  # t_start <- dateToYears(as.Date(MODEL3_INPUT_PARAMETERS$t_start))
  t_start <- dateToYears(as.Date("2010-10-23"))
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
A%s = nearbyint((1-sigma)/sigma * 1/epsilon * cases_at_t_start%s[n_cases_start-1][1]/7 * 365 /(mu + gamma));
I%s = nearbyint(1/epsilon * cases_at_t_start%s[n_cases_start-1][1]/7 * 365 /(mu + alpha + gamma))  ;  // Steady state, DP says its correct.

// Reset R0 vector
R0[0] = 0;
R0[1] = 0;

for(int i = 0; i < n_cases_start; i++){
  R0[0] +=                   cases_at_t_start%s[i][1]/epsilon  * exp((cases_at_t_start%s[i][0] - t_start)  * (rho + mu)); /* because t_i in past so t_ - t_0 negative */
  R0[1] += (1-sigma)/sigma * cases_at_t_start%s[i][1]/epsilon  * exp((cases_at_t_start%s[i][0] - t_start)  * (rho + mu));
}

R1%s = nearbyint((R0[0] + R0[1]) / 3);
R2%s = nearbyint((R0[0] + R0[1]) / 3);
R3%s = nearbyint((R0[0] + R0[1]) / 3);

if (A%s + I%s + R1%s + R2%s + R3%s >= H%s) {
  double R_tot = H%s - A%s - I%s - 100.0;
  if (R_tot <= 0)
  {
  I%s     = nearbyint(H%s - 100);
  A%s     = nearbyint(0);
  R_tot = nearbyint(0);
  }
  R1%s = nearbyint(R_tot / 3);
  R2%s = nearbyint(R_tot / 3);
  R3%s = nearbyint(R_tot / 3);
}

S%s   = nearbyint(H%s - A%s - I%s - R1%s - R2%s - R3%s);
B%s   = (I%s * thetaI/mu_B + A%s * thetaA/mu_B) * D%s * (1 + lambdaR * pow(B0%s, r));
C%s   = 0;
W%s   = 0;

VSd%s = 0;
VR1d%s = 0;
VR2d%s = 0;
VR3d%s = 0;

VSdd%s = 0;
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
double thetaA = thetaI * XthetaA;"

  for (dp in departements) {
    initalizeStatesAll = paste0(initalizeStatesAll, gsub('%s', gsub('-', '_', dp), initalizeStatesTemplate))
  }

  initalizeStates <- pomp::Csnippet(initalizeStatesAll)

  ###
  ### rmeas
  ###

  rmeasTemplate <- "
  if (t > 2018) {
    cases%s = rnbinom_mu(k, epsilon * C%s * cas_def);
  } else {
    cases%s = rnbinom_mu(k, epsilon * C%s);
  }
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
    if (ISNA(cases%s)) {
       lik += (give_log) ? 0 : 1;
    } else {
      if (t > 2018) {
        lik += dnbinom_mu(cases%s, k, epsilon * C%s * cas_def, give_log);
      } else {
        lik += dnbinom_mu(cases%s, k, epsilon * C%s, give_log);
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

  # Every department gets a copy.
  sirb_file <- "
previous_vacc_campaign = TRUE;
r_v_wdn = 0;
dw = 0;

// other_cases = CasesArtibonite + CasesCentre + CasesGrande_Anse + CasesNippes + CasesNord +
//                CasesNord_Est + CasesNord_Ouest + CasesOuest + CasesSud + CasesSud_Est - Cases%s;

// 2019.029
// if (t <= 0) {
//   other_cases = CasesArtibonite + CasesCentre + CasesGrande_Anse + CasesNippes + CasesNord +
//                 CasesNord_Est + CasesNord_Ouest + CasesOuest + CasesSud + CasesSud_Est - Cases%s;
//   if (t >= 2018) {
//     mobility =  (365 * other_cases / (7 * epsilon * cas_def)) * (1 / (mu + alpha + gamma) + (1 - sigma) / (sigma * (mu + gamma)));
//   } else {
//     mobility =  (365 * other_cases / (7 * epsilon)) * (1 / (mu + alpha + gamma) + (1 - sigma) / (sigma * (mu + gamma)));
//   }
// } else {
//   mobility =  (IArtibonite + ICentre + IGrande_Anse + INippes + INord +
//              INord_Est + INord_Ouest + IOuest + ISud + ISud_Est - I%s +
//              AArtibonite + ACentre + AGrande_Anse + ANippes + ANord +
//              ANord_Est + ANord_Ouest + AOuest + ASud + ASud_Est - A%s);
// }

mobility =  (IArtibonite + ICentre + IGrande_Anse + INippes + INord +
             INord_Est + INord_Ouest + IOuest + ISud + ISud_Est - I%s +
             AArtibonite + ACentre + AGrande_Anse + ANippes + ANord +
             ANord_Est + ANord_Ouest + AOuest + ASud + ASud_Est - A%s);


// TODO: Remove when done
// MobAI%s = (IArtibonite + ICentre + IGrande_Anse + INippes + INord +
//              INord_Est + INord_Ouest + IOuest + ISud + ISud_Est - I%s +
//              AArtibonite + ACentre + AGrande_Anse + ANippes + ANord +
//              ANord_Est + ANord_Ouest + AOuest + ASud + ASud_Est - A%s);
// MobCases%s = (365 * other_cases / (7 * epsilon)) * (1 / (mu + alpha + gamma) + (1 - sigma) / (sigma * (mu + gamma)));

// if (t >= 2018) {
//   Mob2%s = (other_cases / (cas_def * epsilon * gamma * 7)) * (1 + 1 / sigma) * 365;
// } else {
//   Mob2%s = (other_cases / (epsilon * gamma * 7)) * (1 + 1 / sigma) * 365;
// }



// force of infection
foi = betaB%s * (B%s / (1 + B%s)) + foi_add%s * mobility;

if(std_W > 0.0) {
  dw = rgammawn(std_W, dt);   // white noise (extra-demographic stochasticity)
  foi_stoc = foi * dw/dt;     // apply stochasticity
} else {
  foi_stoc = foi;
}

if (t <= (t_vacc_end_alt%s + dt)){
	previous_vacc_campaign = TRUE;
	if (t >= t_vacc_start_alt%s && t <= (t_vacc_end_alt%s + dt)) {
    	r_v_wdn = (r_v_year_alt%s / (S%s + A%s + R1%s + R2%s + R3%s));
	}
	p1d = p1d_reg_alt%s;
} else {
	previous_vacc_campaign = FALSE;
	if (t >= t_vacc_start%s && t <= (t_vacc_end%s + dt)) {
    	r_v_wdn = (r_v_year%s / (S%s + A%s + R1%s + R2%s + R3%s));
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
rate[9] = gamma;            // A -> R1
rate[10] = p1d * r_v_wdn;   // A -> VR1d
rate[11] = pdd * r_v_wdn;   // A -> VR1dd

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
rate[33] = sigma * (1 - eff_v_2d(t_eff, scenario)) * foi_stoc;       // VSdd -> I
rate[34] = (1 - sigma) * (1 - eff_v_2d(t_eff, scenario)) * foi_stoc; // VSdd -> A
rate[35] = mu;                                                       // natural death

// VR1dd
rate[36] = mu;          // natural death:    VR1dd -> death
rate[37] = 3 * rho;     // loss of immunity: VR1dd -> VR2dd

// VR2dd
rate[38] = mu;          // natural death:    VR2dd -> death
rate[39] = 3 * rho;     // loss of immunity: VR2dd -> VR3dd

// VR3dd
rate[40] = mu;          // natural death:    VR3dd -> death
rate[41] = 3 * rho;     // loss of immunity: VR3dd -> VSdd


/* For previous vacc campagain */

// VSd_alt
rate[42] = sigma       * (1 - eff_v_1d(t_eff_alt, scenario)) * foi_stoc; // VSd -> I
rate[43] = (1 - sigma) * (1 - eff_v_1d(t_eff_alt, scenario)) * foi_stoc; // VSd -> A
rate[44] = mu;                                                           // natural death

// VSdd_alt
rate[45] = sigma       * (1 - eff_v_2d(t_eff_alt, scenario)) * foi_stoc; // VSdd -> I
rate[46] = (1 - sigma) * (1 - eff_v_2d(t_eff_alt, scenario)) * foi_stoc; // VSdd -> A
rate[47] = mu;                                                           // natural death

// All other alts already defined.

// simulate all transitions
reulermultinom(5, S%s,     &rate[0],  dt, &dN[0]);
reulermultinom(3, I%s,     &rate[5],  dt, &dN[5]);
reulermultinom(4, A%s,     &rate[8],  dt, &dN[8]);
reulermultinom(4, R1%s,    &rate[12], dt, &dN[12]);
reulermultinom(4, R2%s,    &rate[16], dt, &dN[16]);
reulermultinom(4, R3%s,    &rate[20], dt, &dN[20]);

// Vaccinated 1 dose
reulermultinom(3,  VSd%s,   &rate[24], dt, &dN[24]);
reulermultinom(2, VR1d%s,   &rate[27], dt, &dN[27]);
reulermultinom(2, VR2d%s,   &rate[29], dt, &dN[29]);
reulermultinom(2, VR3d%s,   &rate[31], dt, &dN[31]);

// Vaccinated 2 doses
reulermultinom(3,  VSdd%s,   &rate[33], dt, &dN[33]);
reulermultinom(2, VR1dd%s,   &rate[36], dt, &dN[36]);
reulermultinom(2, VR2dd%s,   &rate[38], dt, &dN[38]);
reulermultinom(2, VR3dd%s,   &rate[40], dt, &dN[40]);

/* For the previous vaccination campain */

// Alt Vaccinated 1 dose
reulermultinom(3,  VSd_alt%s,   &rate[42], dt, &dN[42]);
reulermultinom(2, VR1d_alt%s,   &rate[27], dt, &dN[45]);
reulermultinom(2, VR2d_alt%s,   &rate[29], dt, &dN[47]);
reulermultinom(2, VR3d_alt%s,   &rate[31], dt, &dN[49]);

// Alt Vaccinated 2 doses
reulermultinom(3,  VSdd_alt%s,   &rate[45], dt, &dN[51]);
reulermultinom(2, VR1dd_alt%s,   &rate[27], dt, &dN[54]);
reulermultinom(2, VR2dd_alt%s,   &rate[29], dt, &dN[56]);
reulermultinom(2, VR3dd_alt%s,   &rate[31], dt, &dN[58]);

// bacteria as continous state variable
// implement Runge-Kutta integration assuming S, I, R, V* stay constant during dt
k1 = dt * fB(I%s, A%s, B%s, mu_B, thetaI, thetaA, lambdaR, rain_std%s, r, D%s);
k2 = dt * fB(I%s, A%s, B%s, mu_B, thetaI, thetaA, lambdaR, rain_std%s, r, D%s);
k3 = dt * fB(I%s, A%s, B%s, mu_B, thetaI, thetaA, lambdaR, rain_std%s, r, D%s);
k4 = dt * fB(I%s, A%s, B%s, mu_B, thetaI, thetaA, lambdaR, rain_std%s, r, D%s);

// bacteria increment
dB = (k1 + 2*k2 + 2*k3 + k4) / 6.0;

// Update States
I%s  += dN[0] + dN[24] + dN[33] + dN[42] + dN[51] - dN[7] - dN[6] - dN[5];            // S -> I, VSd -> I, VSdd -> I, VSd_alt -> I, VSdd_alt -> I, I -> R, I-> death, I -> death
A%s  += dN[1] + dN[25] + dN[34] + dN[43] + dN[52] - dN[9] - dN[10] - dN[11] - dN[8];  // S -> A, VSd -> A, VSdd -> A, VSd_alt -> A, VSdd_alt -> A, A -> R1, A -> VR1d, A -> VR1dd, natural death.
R1%s += dN[7] + dN[9] - dN[12] - dN[14] - dN[15] - dN[13];  // I-> R1, A -> R1, R1 -> R2, R1 -> VR1d, R1 -> VR1dd, R1 -> death
R2%s += dN[12] - dN[16] - dN[18] - dN[19] - dN[17]; // R1 -> R2, R2 -> R3, R2 -> VR2d, R2 -> VR2dd, R2 -> death
R3%s += dN[16] - dN[20] - dN[22] - dN[23] - dN[21]; // R2 -> R3, R3 -> S , R3 -> VR3d, R3 -> VR3dd, R3 -> death

if (previous_vacc_campaign){
	VSd_alt%s    +=  dN[2];            // S -> VSd_alt
	VR1d_alt%s   +=  dN[14] + dN[10];  // R1 -> VR1d_alt, A -> VR1d_alt
	VR2d_alt%s   +=  dN[18];           // R2 -> VR2d_alt
	VR3d_alt%s   +=  dN[22];           // R3 -> VR3d_alt

	VSdd_alt%s   += dN[3];             // S -> VSdd_alt
	VR1dd_alt%s  += dN[15] + dN[11];   // R1 -> VR1dd_alt, A -> VR1dd_alt
	VR2dd_alt%s  += dN[19];            // R2 -> VR2dd_alt
	VR3dd_alt%s  += dN[23];            // R3 -> VR3dd_alt

} else {
	VSd%s    += dN[2];                 // S -> VSd
	VR1d%s   += dN[14] + dN[10];       // R1 -> VR1d, A -> VR1d
	VR2d%s   += dN[18];                // R2 -> VR2d
	VR3d%s   += dN[22];                // R3 -> VR3d

	VSdd%s   += dN[3];                 // S -> VSdd
	VR1dd%s  += dN[15] + dN[11];       // R1 -> VR1dd, A -> VR1dd
	VR2dd%s  += dN[19];                // R2 -> VR2dd
	VR3dd%s  += dN[23];                // R3 -> VR3dd
}

VSd%s   += dN[32] - dN[24] - dN[25] - dN[26]; // VR3d -> VSd, VSd -> I, VSd -> A, VSd -> death
VR1d%s  += - dN[27] - dN[28];                 // VR1d -> death, VR1d -> VR2d
VR2d%s  += dN[28] - dN[29] - dN[30];          // VR1d -> VR2d, VR2d -> death, VR2d -> VR3d
VR3d%s  += dN[30] - dN[31] - dN[32];          // VR2d -> VR3d, VR3d -> death, VR3d -> VSd

VSdd%s   +=  dN[41] - dN[33] - dN[34] - dN[35]; // VR3dd -> VSdd, VSdd -> I, VSdd -> A, VSdd -> death
VR1dd%s  += - dN[36] - dN[37];                  // VR1dd -> death, VR1dd -> VR2dd
VR2dd%s  +=  dN[37] - dN[38] - dN[39];          // VR1dd -> VR2dd, VR2dd -> death, VR2dd -> VR3dd
VR3dd%s  +=  dN[39] - dN[40] - dN[41];          // VR2dd -> VR3dd, VR3dd -> death, VR3dd -> VSdd

// previous vacccination campain
VSd_alt%s   += dN[50] - dN[42] - dN[43] - dN[44]; // VR3d_alt -> VSd_alt, VSd_alt -> I, VSd_alt -> A, VSd_alt -> death
VR1d_alt%s  += - dN[45] - dN[46];                 // death, VR1d_alt -> VR2d_alt
VR2d_alt%s  += dN[46] - dN[47] - dN[48];          // VR1d_alt -> VR2d_alt, death, VR2d_alt -> VR3d_alt
VR3d_alt%s  += dN[48] - dN[49] - dN[50];          // VR2d_alt -> VR3d_alt, death, VR3d_alt -> VSd_alt

VSdd_alt%s   +=  dN[59] - dN[51] - dN[52] - dN[53]; // VR3dd_alt -> VSdd_alt, VSdd_alt -> I, VSdd_alt -> A, VSdd_alt -> death
VR1dd_alt%s  += - dN[54] - dN[55];
VR2dd_alt%s  +=  dN[55] - dN[56] - dN[57] ;
VR3dd_alt%s  +=  dN[57] - dN[58] - dN[59];

C%s   +=  dN[0] + dN[24] + dN[33] + dN[42] + dN[51]; // S -> I, VSd -> I, VSdd -> I, VSd_alt -> I, VSdd_alt -> I
W%s   +=  (dw - dt) / std_W;  // standardized i.i.d. white noise
B%s   += (((dB) < -B%s) ? (-B%s + 1.0e-3) : (dB)); // condition to ensure B>0

// susceptibles so as to match total population
S%s = nearbyint(H%s - I%s - A%s - R1%s - R2%s - R3%s -
	VSd%s - VR1d%s - VR2d%s - VR3d%s -
	VSdd%s- VR1dd%s -VR2dd%s -VR3dd%s -
	VSd_alt%s - VR1d_alt%s - VR2d_alt%s - VR3d_alt%s -
	VSdd_alt%s - VR1dd_alt%s - VR2dd_alt%s - VR3dd_alt%s);

IncidenceAll +=  dN[0] + dN[24] + dN[33] + dN[42] + dN[51] + dN[1] + dN[25] + dN[34] + dN[43] + dN[52];
if (!previous_vacc_campaign)
{
	DosesAll  += dN[2] + dN[10] + dN[14] + dN[18] + dN[22] + 2 * (dN[3] + dN[11] + dN[15] + dN[19] + dN[23]);
}
CasesAll  +=  dN[0] + dN[24] + dN[33] + dN[42] + dN[51];
  "

rproc_main <- "
double foi, foi_stoc;   // force of infection and its stochastic version
double dw;              // extra-demographic stochasticity on foi
double dB;              // deterministic forward time difference of bacteria in the environment
double k1, k2, k3, k4;  // coefficients of  the Runge-Kutta method
double rate[48];        // vector of all rates in model
double dN[60];          // vector of transitions between classes during integration timestep
double mobility;
double p1d, pdd;
double r_v_wdn = 0.0;   // rate of vaccination: 0 if out of time window, r_v if not
int previous_vacc_campaign; // flag that indicate if we are on the first or second campain
int scenario =  cases_ext;
double t_eff, t_eff_alt;
double other_cases; // TODO: remove this when done

double thetaA = thetaI * XthetaA;

  "

  final_rproc <- ""
  for (dp in departements) {
    final_rproc <- paste0(final_rproc, gsub('%s', gsub('-', '_', dp), sirb_file))
  }

  final_rproc <- paste0(rproc_main, final_rproc)
  final_rproc_c <- pomp::Csnippet(final_rproc)

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

  # dt_yrs <- 1 / 365.25 * .2

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
  all_params["rho"] <- 1 / (365 * 8) * 365.25   # Fixed  \rho in table 14
  all_params["gamma"] <- 182.625          # TODO: where does this number come from? I think this is a mistake by the authors. Fixed  divide by 365 to get \gamma in table 14
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

  all_params["B0Artibonite"] =   0.24
  all_params["B0Sud_Est"] =      0.24
  all_params["B0Nippes"] =       0.24
  all_params["B0Nord_Est"] =     0.24
  all_params["B0Ouest"] =        0.24
  all_params["B0Centre"] =       0.24
  all_params["B0Nord"] =         0.24
  all_params["B0Sud"] =          0.24
  all_params["B0Nord_Ouest"] =   0.24
  all_params["B0Grande_Anse"] =  0.24

  cases_df <- all_cases %>%
    dplyr::filter(time > t_start & time < (t_end + 0.01)) %>%
    dplyr::select(
      time, casesArtibonite, casesCentre,
      casesGrande_Anse, casesNippes,
      casesNord, casesNord_Est, casesOuest,
      casesSud, casesSud_Est, casesNord_Ouest
    )

  tot_cases <- cases_df %>%
    dplyr::rename(
      CasesArtibonite = casesArtibonite,
      CasesCentre = casesCentre,
      CasesGrande_Anse = casesGrande_Anse,
      CasesNippes = casesNippes,
      CasesNord = casesNord,
      CasesNord_Est = casesNord_Est,
      CasesOuest = casesOuest,
      CasesSud = casesSud,
      CasesSud_Est = casesSud_Est,
      CasesNord_Ouest = casesNord_Ouest
    )


  all_rain <- all_rain %>% dplyr::select(time, dplyr::starts_with('rain_std'))
  # covars <- dplyr::full_join(all_rain, tot_cases, by = 'time') %>%
  #   tidyr::fill(CasesArtibonite,
  #               CasesCentre,
  #               CasesGrande_Anse,
  #               CasesNippes,
  #               CasesNord,
  #               CasesNord_Est,
  #               CasesOuest,
  #               CasesSud,
  #               CasesSud_Est,
  #               CasesNord_Ouest, .direction = c('updown'))
  covars <- all_rain %>%
    dplyr::select(time, dplyr::starts_with('rain_std'))

  # %>%
    # tidyr::fill(dplyr::starts_with("Cases"), .direction = 'up')

  sirb_cholera <- pomp::pomp(
    # set data
    data = cases_df,
    # time column
    times = "time",
    # initialization time
    t0 = t_start - dt_yrs,
    # paramter vector
    params = all_params,

    # process simulator
    rprocess = pomp::euler(step.fun = final_rproc_c, delta.t = dt_yrs),
    # measurement model simulator
    rmeasure =  rmeas,
    # measurement model density
    dmeasure = dmeas,
    # covariates
    covar = pomp::covariate_table(covars, times = 'time'),
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
      sprintf("int n_cases_other = %i;",  nrow(cases_other_dept)), # TODO: Check if this can be removed
      sprintf("double t_start = %f;",  t_start),
      sprintf("double t_end = %f;",  t_end),
      derivativeBacteria.c,
      all_matrix_cases_at_t_start.string,
      all_matrix_cases_other.string,  # TODO: Problably don't need these either
      eff_v.c,
      sep = " "
    )
  )

  return(sirb_cholera)
}
