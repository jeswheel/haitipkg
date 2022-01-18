#' Build pomp object for individual departments following model 3
#'
#' @param departement String indicating which departement to fit the
#'   pomp model to.
#' @param start_time Time of which to start the time series. In the
#'    original paper, the authors used a start time of 2014-03-01 for all
#'    departements except for Ouest, which started with a time of 2017-06-10.
#'    This means that there wasn't a national-coupled model fit to the
#'    data until 2017-06-01, which means that there was less than 2 years
#'    of data to fit the model. Still, the authors fit the remaining departemental
#'    models to data from 2014-03-01. This added argument allows us to
#'    fit the Ouest model earlier, or fit the remaining departement models
#'    at 2017-06-01
#' @param delta.t Delta time step used in Euler's approximation
#' @param betaB_trend Whether or not to include a trend parameter for
#'    betaB.
#' @param B0 Boolean indicator. If TRUE, Initial value of the Bacteria
#'    compartment will be estimated; If FALSE, the equilibrium will be used.
#'
#' @return \code{\link[pomp]{pomp}} object.
#'
#' @importFrom magrittr %>%
#' @importFrom foreach %do%
#' @importFrom foreach foreach
#' @return \code{\link[pomp]{pomp}} representation of model 3 described in \href{https://www.sciencedirect.com/science/article/pii/S2214109X20303107}{Lee, Elizabeth et. al.} and it's accompanying \href{https://ars.els-cdn.com/content/image/1-s2.0-S2214109X20303107-mmc3.pdf}{Supplemental Material}.
#' @export

haiti3_dep <- function(departement = 'Artibonite',
                       start_time = "2014-03-01",
                       delta.t = 1/365.25 * .2,
                       betaB_trend = FALSE,
                       B0 = FALSE) {

  # function to convert dates to fractions of years for model
  dateToYears <- function(date, origin = as.Date("2014-01-01"), yr_offset = 2014) {
    julian(date, origin = origin)/365.25 + yr_offset
  }

  yearsToDate <- function(year_frac, origin = as.Date("2014-01-01"), yr_offset = 2014.0) {
    as.Date((year_frac - yr_offset) * 365.25, origin = origin)
  }

  yearsToDateTime <- function(year_frac, origin = as.Date("2014-01-01"), yr_offset = 2014.0) {
    as.POSIXct((year_frac - yr_offset) * 365.25 * 3600 * 24, origin = origin)
  }

  DATA <- haiti3_dep_data(departement = departement, start_time = start_time)

  cases <- DATA$cases
  rain <- DATA$rain
  cases_covar <- DATA$cases_covar
  input_parameters <- MODEL3_INPUT_PARAMETERS

  t_start <- dateToYears(as.Date(start_time))

  # if (departement == 'Ouest' && start_time != '2014-03-01') {
  #   t_start <- dateToYears(as.Date('2017-06-10'))
  # } else {
  #   t_start <- dateToYears(as.Date(input_parameters$t_start))
  # }
  t_end <- dateToYears(as.Date(input_parameters$t_end))

  state_names <- c("S", "I", "A", "RI1", "RI2", "RI3", "RA1", "RA2", "RA3",
                   "VSd", "VRI1d", "VRI2d", "VRI3d", "VRA1d", "VRA2d", "VRA3d",
                   "VSdd", "VRI1dd", "VRI2dd", "VRI3dd", "VRA1dd", "VRA2dd", "VRA3dd",
                   "VSd_alt", "VRI1d_alt", "VRI2d_alt", "VRI3d_alt", "VRA1d_alt", "VRA2d_alt", "VRA3d_alt",
                   "VSdd_alt", "VRI1dd_alt", "VRI2dd_alt", "VRI3dd_alt", "VRA1dd_alt", "VRA2dd_alt", "VRA3dd_alt",
                   "B", "C", "W")

  dmeas <- pomp::Csnippet("
  double mean_cases = epsilon * C;
  // Rprintf(\"%f\\n\", (betaB + betaB_trend * (t - 2016.579)) * (B / (1 + B)));
  // Rprintf(\"%f\\n\", betaB * exp(betaB_trend * (t - 2016.579)) * (B / (1 + B)));
  if (t > 2018)
    mean_cases = mean_cases * cas_def;
  if (ISNA(cases)) {
    lik = (give_log) ? 0 : 1;
    } else {
        if (S < 10000) {
          lik = (give_log) ? -99999 : 1.0e-18;
          } else {
              // Rprintf(\"cases=%f, mean_cases=%f\\n\",cases,mean_cases);
              lik = dnbinom_mu(cases, k, mean_cases, give_log) ;
          }
      }
      ")

  rmeas <- pomp::Csnippet("
  double mean_cases = epsilon * C;
  if (t > 2018)
    mean_cases = mean_cases * cas_def;
  // cases = mean_cases;
  cases = rnbinom_mu(k, mean_cases);
  ")


  if (betaB_trend) {
    sirb.rproc <- pomp::Csnippet("
  double foi, foi_stoc; // force of infection and its stochastic version
double dw;            // extra-demographic stochasticity on foi
double dB;            // deterministic forward time difference of bacteria in the environment
double k1, k2, k3, k4;  // coefficients of  the Runge-Kutta method
double rate[39];      // vector of all rates in model
double dN[94];        // vector of transitions between classes during integration timestep

double thetaA = thetaI * XthetaA;
double rhoI = rhoA * XrhoI;

int previous_vacc_campaign = TRUE ; /* flag that indicate if we are on the first or second campain */
double r_v_wdn = 0.0;       // rate of vaccination: 0 if out of time window, r_v if not
double p1d = 0;
double mobility = 0;

int scenario = 1;

scenario = cases_ext;



/* Compute mobility term */
if (t < t_end) {
      for(int i = 0; i < n_cases_start - 1; i++){
           if (t >= cases_other[i][0] && t <= cases_other[i+1][0])
               mobility = cases_other[i][1];
      }
      if (t > cases_other[n_cases_start-1][0])
          mobility = cases_other[n_cases_start-1][1];
}
else {
    mobility = cases_covar_c;
}


if (t < 2018){
	mobility = mobility / epsilon;
}
else
{
	mobility = mobility / (epsilon*cas_def);
}

// force of infection
// foi = betaB * (B / (1 + B)) + foi_add*mobility;
// foi = (betaB + betaB_trend * (t - 2016.579)) * (B / (1 + B)) + foi_add*mobility;
foi = betaB * exp(betaB_trend * (t - 2016.579)) * (B / (1 + B)) + foi_add*mobility;
if(std_W > 0.0)
{
    dw = rgammawn(std_W, dt);   // white noise (extra-demographic stochasticity)
    foi_stoc = foi * dw/dt;      // apply stochasticity
} else
{
    foi_stoc = foi;
}

if (t <= (t_vacc_end_alt + dt)){
	previous_vacc_campaign = TRUE;
	if (t >= t_vacc_start_alt && t <= (t_vacc_end_alt + dt)) {
    	r_v_wdn = (r_v_alt / (S + A + RI1 + RI2 + RI3 + RA1 + RA2 + RA3));
	}
	p1d = p1d_alt;
} else {
	previous_vacc_campaign = FALSE;
	if (t >= t_vacc_start && t <= (t_vacc_end + dt)) {
    	r_v_wdn = (r_v_year / (S + A + RI1 + RI2 + RI3 + RA1 + RA2 + RA3));
	}
	p1d = p1d_reg;
}
double pdd = 1 - p1d;
//if ((S + A + RI1 + RI2 + RI3 + RA1 + RA2 + RA3) < 1000)
//	r_v_wdn = 0;


// time in the vacc_eff referential. We assume different timing for 1d and 2d
double t_eff =     t - (t_vacc_start + (t_vacc_end - t_vacc_start)/2);
double t_eff_alt = t - (t_vacc_start_alt + (t_vacc_end_alt - t_vacc_start_alt)/2);

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
reulermultinom(4, S,     &rate[0],  dt, &dN[0]);
reulermultinom(3, I,     &rate[4],  dt, &dN[4]);
reulermultinom(4, A,     &rate[7],  dt, &dN[7]);
reulermultinom(4, RI1,   &rate[11], dt, &dN[11]);
reulermultinom(4, RI2,   &rate[11], dt, &dN[15]);
reulermultinom(4, RI3,   &rate[11], dt, &dN[19]);
reulermultinom(4, RA1,   &rate[15], dt, &dN[23]);
reulermultinom(4, RA2,   &rate[15], dt, &dN[27]);
reulermultinom(4, RA3,   &rate[15], dt, &dN[31]);
/* Vaccinated 1 dose */
reulermultinom(3, VSd,   &rate[19], dt, &dN[35]);
reulermultinom(2, VRI1d, &rate[22], dt, &dN[38]);
reulermultinom(2, VRI2d, &rate[22], dt, &dN[40]);
reulermultinom(2, VRI3d, &rate[22], dt, &dN[42]);
reulermultinom(2, VRA1d, &rate[24], dt, &dN[44]);
reulermultinom(2, VRA2d, &rate[24], dt, &dN[46]);
reulermultinom(2, VRA3d, &rate[24], dt, &dN[48]);
/* Vaccinated 2 doses */
reulermultinom(3, VSdd,  &rate[26], dt, &dN[50]);
reulermultinom(2, VRI1dd,&rate[22], dt, &dN[53]);
reulermultinom(2, VRI2dd,&rate[22], dt, &dN[55]);
reulermultinom(2, VRI3dd,&rate[22], dt, &dN[57]);
reulermultinom(2, VRA1dd,&rate[24], dt, &dN[59]);
reulermultinom(2, VRA2dd,&rate[24], dt, &dN[61]);
reulermultinom(2, VRA3dd,&rate[24], dt, &dN[63]);
/* For the previous vaccination campain */
/* Vaccinated 1 dose */
reulermultinom(3, VSd_alt,   &rate[29], dt, &dN[65]);
reulermultinom(2, VRI1d_alt, &rate[32], dt, &dN[68]);
reulermultinom(2, VRI2d_alt, &rate[32], dt, &dN[70]);
reulermultinom(2, VRI3d_alt, &rate[32], dt, &dN[72]);
reulermultinom(2, VRA1d_alt, &rate[34], dt, &dN[74]);
reulermultinom(2, VRA2d_alt, &rate[34], dt, &dN[76]);
reulermultinom(2, VRA3d_alt, &rate[34], dt, &dN[78]);
/* Vaccinated 2 doses */
reulermultinom(3, VSdd_alt,  &rate[36], dt, &dN[80]);
reulermultinom(2, VRI1dd_alt,&rate[32], dt, &dN[83]);
reulermultinom(2, VRI2dd_alt,&rate[32], dt, &dN[85]);
reulermultinom(2, VRI3dd_alt,&rate[32], dt, &dN[87]);
reulermultinom(2, VRA1dd_alt,&rate[34], dt, &dN[89]);
reulermultinom(2, VRA2dd_alt,&rate[34], dt, &dN[91]);
reulermultinom(2, VRA3dd_alt,&rate[34], dt, &dN[93]);

// bacteria as continous state variable
// implement Runge-Kutta integration assuming S, I, R, V* stay constant during dt
k1 = dt * fB(I, A, B, mu_B, thetaI, thetaA, lambdaR, rain, r, D);
k2 = dt * fB(I, A, B, mu_B, thetaI, thetaA, lambdaR, rain, r, D);
k3 = dt * fB(I, A, B, mu_B, thetaI, thetaA, lambdaR, rain, r, D);
k4 = dt * fB(I, A, B, mu_B, thetaI, thetaA, lambdaR, rain, r, D);
// bacteria increment
dB = (k1 + 2*k2 + 2*k3 + k4) / 6.0;


// update state variables
//S = -dN[0] - dN[1] /* FOI */
//- dN[2] - dN[3] /* VACC */
//+ dN[4] + dN[7] + dN[12] + dN[16] + dN[20] + dN[24] + dN[28] + dN[32] /* Mortality of nn vacc: rebirth*/
//+ dN[19] + dN[31]  /* Recovery */
//+ dN[37] + dN[38] + dN[40] + dN[42] + dN[44] + dN[46] + dN[48]        /* Mortality of 1d vacc: rebirth*/
//+ dN[52] + dN[53] + dN[55] + dN[57] + dN[59] + dN[61] + dN[63]		   /* Mortality of dd vacc: rebirth*/
//+ dN[37+30] + dN[38+30] + dN[40+30] + dN[42+30] + dN[44+30] + dN[46+30] + dN[48+30]        /* Prev campain*/
//+ dN[52+30] + dN[53+30] + dN[55+30] + dN[57+30] + dN[59+30] + dN[61+30] + dN[63+30];

I   += dN[0] + dN[35] + dN[50] + dN[35+30] + dN[50+30] - dN[4] - dN[5] - dN[6];
A   += dN[1] + dN[36] + dN[51] + dN[36+30] + dN[51+30] - dN[7] - dN[8] - dN[9] - dN[10];
RI1 += dN[6] -  dN[11] - dN[12] - dN[13] - dN[14] ;
RI2 += dN[11] - dN[15] - dN[16] - dN[17] - dN[18];
RI3 += dN[15] - dN[19] - dN[20] - dN[21] - dN[22];
RA1 += dN[8]  - dN[23] - dN[24] - dN[25] - dN[26];
RA2 += dN[23] - dN[27] - dN[28] - dN[29] - dN[30];
RA3 += dN[27] - dN[31] - dN[32] - dN[33] - dN[34];

if (previous_vacc_campaign){
	VSd_alt    += dN[2];
	VRI1d_alt  += dN[13];
	VRI2d_alt  += dN[17];
	VRI3d_alt  += dN[21];
	VRA1d_alt  += dN[9] + dN[25] ;
	VRA2d_alt  += dN[29];
	VRA3d_alt  += dN[33];

	VSdd_alt    += dN[3];
	VRI1dd_alt  += dN[14];
	VRI2dd_alt  += dN[18];
	VRI3dd_alt  += dN[22];
	VRA1dd_alt  += dN[10] + dN[26];
	VRA2dd_alt  += dN[30];
	VRA3dd_alt  += dN[34];
} else {
	VSd    += dN[2];
	VRI1d  += dN[13];
	VRI2d  += dN[17];
	VRI3d  += dN[21];
	VRA1d  += dN[9] + dN[25] ;
	VRA2d  += dN[29];
	VRA3d  += dN[33];

	VSdd    += dN[3];
	VRI1dd  += dN[14];
	VRI2dd  += dN[18];
	VRI3dd  += dN[22];
	VRA1dd  += dN[10] + dN[26];
	VRA2dd  += dN[30];
	VRA3dd  += dN[34];

}

VSd    += dN[43] + dN[49] - dN[35] - dN[36] - dN[37];
VRI1d  += - dN[38] - dN[39];
VRI2d  += dN[39] - dN[40] - dN[41];
VRI3d  += dN[41] - dN[42] - dN[43];
VRA1d  += - dN[44] - dN[45];
VRA2d  += dN[45] - dN[46] - dN[47];
VRA3d  += dN[47] - dN[48] - dN[49];

VSdd    +=  dN[58] + dN[64] - dN[50] - dN[51] - dN[52];
VRI1dd  += - dN[53] - dN[54];
VRI2dd  +=  dN[54] - dN[55] - dN[56] ;
VRI3dd  +=  dN[56] - dN[57] - dN[58];
VRA1dd  += - dN[59] - dN[60];
VRA2dd  +=  dN[60] - dN[61] - dN[62];
VRA3dd  +=  dN[62] - dN[63] - dN[64];

/* *previous* vacccination campain */

VSd_alt    += dN[43+30] + dN[49+30] - dN[35+30] - dN[36+30] - dN[37+30];
VRI1d_alt  += - dN[38+30] - dN[39+30];
VRI2d_alt  += dN[39+30] - dN[40+30] - dN[41+30];
VRI3d_alt  += dN[41+30] - dN[42+30] - dN[43+30];
VRA1d_alt  += - dN[44+30] - dN[45+30];
VRA2d_alt  += dN[45+30] - dN[46+30] - dN[47+30];
VRA3d_alt  += dN[47+30] - dN[48+30] - dN[49+30];

VSdd_alt    +=  dN[58+30] + dN[64+30] - dN[50+30] - dN[51+30] - dN[52+30];
VRI1dd_alt  += - dN[53+30] - dN[54+30];
VRI2dd_alt  +=  dN[54+30] - dN[55+30] - dN[56+30] ;
VRI3dd_alt  +=  dN[56+30] - dN[57+30] - dN[58+30];
VRA1dd_alt  += - dN[59+30] - dN[60+30];
VRA2dd_alt  +=  dN[60+30] - dN[61+30] - dN[62+30];
VRA3dd_alt  +=  dN[62+30] - dN[63+30] - dN[64+30];


C   +=  dN[0] + dN[35] + dN[50] + dN[35+30] + dN[50+30];
W   +=  (dw - dt)/std_W;  // standardized i.i.d. white noise
B += (((dB) < -B) ? (-B + 1.0e-3) : (dB)); // condition to ensure B>0

// susceptibles so as to match total population
S = nearbyint(H - I - A - RI1 - RI2 - RI3 - RA1 - RA2 - RA3 -
	VSd - VRI1d - VRI2d - VRI3d - VRA1d - VRA2d -VRA3d -
	VSdd- VRI1dd -VRI2dd -VRI3dd -VRA1dd-VRA2dd-VRA3dd -
	VSd_alt - VRI1d_alt - VRI2d_alt - VRI3d_alt - VRA1d_alt - VRA2d_alt - VRA3d_alt -
	VSdd_alt - VRI1dd_alt - VRI2dd_alt - VRI3dd_alt - VRA1dd_alt - VRA2dd_alt - VRA3dd_alt);
                         ")
  } else {
    sirb.rproc <- pomp::Csnippet("
  double foi, foi_stoc; // force of infection and its stochastic version
double dw;            // extra-demographic stochasticity on foi
double dB;            // deterministic forward time difference of bacteria in the environment
double k1, k2, k3, k4;  // coefficients of  the Runge-Kutta method
double rate[39];      // vector of all rates in model
double dN[94];        // vector of transitions between classes during integration timestep

double thetaA = thetaI * XthetaA;
double rhoI = rhoA * XrhoI;

int previous_vacc_campaign = TRUE ; /* flag that indicate if we are on the first or second campain */
double r_v_wdn = 0.0;       // rate of vaccination: 0 if out of time window, r_v if not
double p1d = 0;
double mobility = 0;

int scenario = 1;

scenario = cases_ext;



/* Compute mobility term */
if (t < t_end) {
      for(int i = 0; i < n_cases_start - 1; i++){
           if (t >= cases_other[i][0] && t <= cases_other[i+1][0])
               mobility = cases_other[i][1];
      }
      if (t > cases_other[n_cases_start-1][0])
          mobility = cases_other[n_cases_start-1][1];
}
else {
    mobility = cases_covar_c;
}


if (t < 2018){
	mobility = mobility / epsilon;
}
else
{
	mobility = mobility / (epsilon*cas_def);
}

// force of infection
foi = betaB * (B / (1 + B)) + foi_add*mobility;
// foi = (betaB + betaB_trend * (t - 2016.579)) * (B / (1 + B)) + foi_add*mobility;
if(std_W > 0.0)
{
    dw = rgammawn(std_W, dt);   // white noise (extra-demographic stochasticity)
    foi_stoc = foi * dw/dt;      // apply stochasticity
} else
{
    foi_stoc = foi;
}

if (t <= (t_vacc_end_alt + dt)){
	previous_vacc_campaign = TRUE;
	if (t >= t_vacc_start_alt && t <= (t_vacc_end_alt + dt)) {
    	r_v_wdn = (r_v_alt / (S + A + RI1 + RI2 + RI3 + RA1 + RA2 + RA3));
	}
	p1d = p1d_alt;
} else {
	previous_vacc_campaign = FALSE;
	if (t >= t_vacc_start && t <= (t_vacc_end + dt)) {
    	r_v_wdn = (r_v_year / (S + A + RI1 + RI2 + RI3 + RA1 + RA2 + RA3));
	}
	p1d = p1d_reg;
}
double pdd = 1 - p1d;
//if ((S + A + RI1 + RI2 + RI3 + RA1 + RA2 + RA3) < 1000)
//	r_v_wdn = 0;


// time in the vacc_eff referential. We assume different timing for 1d and 2d
double t_eff =     t - (t_vacc_start + (t_vacc_end - t_vacc_start)/2);
double t_eff_alt = t - (t_vacc_start_alt + (t_vacc_end_alt - t_vacc_start_alt)/2);

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
reulermultinom(4, S,     &rate[0],  dt, &dN[0]);
reulermultinom(3, I,     &rate[4],  dt, &dN[4]);
reulermultinom(4, A,     &rate[7],  dt, &dN[7]);
reulermultinom(4, RI1,   &rate[11], dt, &dN[11]);
reulermultinom(4, RI2,   &rate[11], dt, &dN[15]);
reulermultinom(4, RI3,   &rate[11], dt, &dN[19]);
reulermultinom(4, RA1,   &rate[15], dt, &dN[23]);
reulermultinom(4, RA2,   &rate[15], dt, &dN[27]);
reulermultinom(4, RA3,   &rate[15], dt, &dN[31]);
/* Vaccinated 1 dose */
reulermultinom(3, VSd,   &rate[19], dt, &dN[35]);
reulermultinom(2, VRI1d, &rate[22], dt, &dN[38]);
reulermultinom(2, VRI2d, &rate[22], dt, &dN[40]);
reulermultinom(2, VRI3d, &rate[22], dt, &dN[42]);
reulermultinom(2, VRA1d, &rate[24], dt, &dN[44]);
reulermultinom(2, VRA2d, &rate[24], dt, &dN[46]);
reulermultinom(2, VRA3d, &rate[24], dt, &dN[48]);
/* Vaccinated 2 doses */
reulermultinom(3, VSdd,  &rate[26], dt, &dN[50]);
reulermultinom(2, VRI1dd,&rate[22], dt, &dN[53]);
reulermultinom(2, VRI2dd,&rate[22], dt, &dN[55]);
reulermultinom(2, VRI3dd,&rate[22], dt, &dN[57]);
reulermultinom(2, VRA1dd,&rate[24], dt, &dN[59]);
reulermultinom(2, VRA2dd,&rate[24], dt, &dN[61]);
reulermultinom(2, VRA3dd,&rate[24], dt, &dN[63]);
/* For the previous vaccination campain */
/* Vaccinated 1 dose */
reulermultinom(3, VSd_alt,   &rate[29], dt, &dN[65]);
reulermultinom(2, VRI1d_alt, &rate[32], dt, &dN[68]);
reulermultinom(2, VRI2d_alt, &rate[32], dt, &dN[70]);
reulermultinom(2, VRI3d_alt, &rate[32], dt, &dN[72]);
reulermultinom(2, VRA1d_alt, &rate[34], dt, &dN[74]);
reulermultinom(2, VRA2d_alt, &rate[34], dt, &dN[76]);
reulermultinom(2, VRA3d_alt, &rate[34], dt, &dN[78]);
/* Vaccinated 2 doses */
reulermultinom(3, VSdd_alt,  &rate[36], dt, &dN[80]);
reulermultinom(2, VRI1dd_alt,&rate[32], dt, &dN[83]);
reulermultinom(2, VRI2dd_alt,&rate[32], dt, &dN[85]);
reulermultinom(2, VRI3dd_alt,&rate[32], dt, &dN[87]);
reulermultinom(2, VRA1dd_alt,&rate[34], dt, &dN[89]);
reulermultinom(2, VRA2dd_alt,&rate[34], dt, &dN[91]);
reulermultinom(2, VRA3dd_alt,&rate[34], dt, &dN[93]);

// bacteria as continous state variable
// implement Runge-Kutta integration assuming S, I, R, V* stay constant during dt
k1 = dt * fB(I, A, B, mu_B, thetaI, thetaA, lambdaR, rain, r, D);
k2 = dt * fB(I, A, B, mu_B, thetaI, thetaA, lambdaR, rain, r, D);
k3 = dt * fB(I, A, B, mu_B, thetaI, thetaA, lambdaR, rain, r, D);
k4 = dt * fB(I, A, B, mu_B, thetaI, thetaA, lambdaR, rain, r, D);
// bacteria increment
dB = (k1 + 2*k2 + 2*k3 + k4) / 6.0;


// update state variables
//S = -dN[0] - dN[1] /* FOI */
//- dN[2] - dN[3] /* VACC */
//+ dN[4] + dN[7] + dN[12] + dN[16] + dN[20] + dN[24] + dN[28] + dN[32] /* Mortality of nn vacc: rebirth*/
//+ dN[19] + dN[31]  /* Recovery */
//+ dN[37] + dN[38] + dN[40] + dN[42] + dN[44] + dN[46] + dN[48]        /* Mortality of 1d vacc: rebirth*/
//+ dN[52] + dN[53] + dN[55] + dN[57] + dN[59] + dN[61] + dN[63]		   /* Mortality of dd vacc: rebirth*/
//+ dN[37+30] + dN[38+30] + dN[40+30] + dN[42+30] + dN[44+30] + dN[46+30] + dN[48+30]        /* Prev campain*/
//+ dN[52+30] + dN[53+30] + dN[55+30] + dN[57+30] + dN[59+30] + dN[61+30] + dN[63+30];

I   += dN[0] + dN[35] + dN[50] + dN[35+30] + dN[50+30] - dN[4] - dN[5] - dN[6];
A   += dN[1] + dN[36] + dN[51] + dN[36+30] + dN[51+30] - dN[7] - dN[8] - dN[9] - dN[10];
RI1 += dN[6] -  dN[11] - dN[12] - dN[13] - dN[14] ;
RI2 += dN[11] - dN[15] - dN[16] - dN[17] - dN[18];
RI3 += dN[15] - dN[19] - dN[20] - dN[21] - dN[22];
RA1 += dN[8]  - dN[23] - dN[24] - dN[25] - dN[26];
RA2 += dN[23] - dN[27] - dN[28] - dN[29] - dN[30];
RA3 += dN[27] - dN[31] - dN[32] - dN[33] - dN[34];

if (previous_vacc_campaign){
	VSd_alt    += dN[2];
	VRI1d_alt  += dN[13];
	VRI2d_alt  += dN[17];
	VRI3d_alt  += dN[21];
	VRA1d_alt  += dN[9] + dN[25] ;
	VRA2d_alt  += dN[29];
	VRA3d_alt  += dN[33];

	VSdd_alt    += dN[3];
	VRI1dd_alt  += dN[14];
	VRI2dd_alt  += dN[18];
	VRI3dd_alt  += dN[22];
	VRA1dd_alt  += dN[10] + dN[26];
	VRA2dd_alt  += dN[30];
	VRA3dd_alt  += dN[34];
} else {
	VSd    += dN[2];
	VRI1d  += dN[13];
	VRI2d  += dN[17];
	VRI3d  += dN[21];
	VRA1d  += dN[9] + dN[25] ;
	VRA2d  += dN[29];
	VRA3d  += dN[33];

	VSdd    += dN[3];
	VRI1dd  += dN[14];
	VRI2dd  += dN[18];
	VRI3dd  += dN[22];
	VRA1dd  += dN[10] + dN[26];
	VRA2dd  += dN[30];
	VRA3dd  += dN[34];

}

VSd    += dN[43] + dN[49] - dN[35] - dN[36] - dN[37];
VRI1d  += - dN[38] - dN[39];
VRI2d  += dN[39] - dN[40] - dN[41];
VRI3d  += dN[41] - dN[42] - dN[43];
VRA1d  += - dN[44] - dN[45];
VRA2d  += dN[45] - dN[46] - dN[47];
VRA3d  += dN[47] - dN[48] - dN[49];

VSdd    +=  dN[58] + dN[64] - dN[50] - dN[51] - dN[52];
VRI1dd  += - dN[53] - dN[54];
VRI2dd  +=  dN[54] - dN[55] - dN[56] ;
VRI3dd  +=  dN[56] - dN[57] - dN[58];
VRA1dd  += - dN[59] - dN[60];
VRA2dd  +=  dN[60] - dN[61] - dN[62];
VRA3dd  +=  dN[62] - dN[63] - dN[64];

/* *previous* vacccination campain */

VSd_alt    += dN[43+30] + dN[49+30] - dN[35+30] - dN[36+30] - dN[37+30];
VRI1d_alt  += - dN[38+30] - dN[39+30];
VRI2d_alt  += dN[39+30] - dN[40+30] - dN[41+30];
VRI3d_alt  += dN[41+30] - dN[42+30] - dN[43+30];
VRA1d_alt  += - dN[44+30] - dN[45+30];
VRA2d_alt  += dN[45+30] - dN[46+30] - dN[47+30];
VRA3d_alt  += dN[47+30] - dN[48+30] - dN[49+30];

VSdd_alt    +=  dN[58+30] + dN[64+30] - dN[50+30] - dN[51+30] - dN[52+30];
VRI1dd_alt  += - dN[53+30] - dN[54+30];
VRI2dd_alt  +=  dN[54+30] - dN[55+30] - dN[56+30] ;
VRI3dd_alt  +=  dN[56+30] - dN[57+30] - dN[58+30];
VRA1dd_alt  += - dN[59+30] - dN[60+30];
VRA2dd_alt  +=  dN[60+30] - dN[61+30] - dN[62+30];
VRA3dd_alt  +=  dN[62+30] - dN[63+30] - dN[64+30];


C   +=  dN[0] + dN[35] + dN[50] + dN[35+30] + dN[50+30];
W   +=  (dw - dt)/std_W;  // standardized i.i.d. white noise
B += (((dB) < -B) ? (-B + 1.0e-3) : (dB)); // condition to ensure B>0

// susceptibles so as to match total population
S = nearbyint(H - I - A - RI1 - RI2 - RI3 - RA1 - RA2 - RA3 -
	VSd - VRI1d - VRI2d - VRI3d - VRA1d - VRA2d -VRA3d -
	VSdd- VRI1dd -VRI2dd -VRI3dd -VRA1dd-VRA2dd-VRA3dd -
	VSd_alt - VRI1d_alt - VRI2d_alt - VRI3d_alt - VRA1d_alt - VRA2d_alt - VRA3d_alt -
	VSdd_alt - VRI1dd_alt - VRI2dd_alt - VRI3dd_alt - VRA1dd_alt - VRA2dd_alt - VRA3dd_alt);
                         ")
  }



  derivativeBacteria.c <- " double fB(int I, int A, double B,
    double mu_B, double thetaI, double XthetaA, double lambdaR, double rain, double r, double D) {
  double thetaA = thetaI * XthetaA;
  double dB;
  dB = -mu_B * B +  (1 + lambdaR * pow(rain, r)) * D * (thetaI * (double) I + thetaA * (double) A);
  return(dB);
};
"


  cases_at_t_start <- cases %>% dplyr::filter(dateToYears(date) <= t_start)
  cases_at_t_start.string <- foreach(
    r = iterators::iter(cases_at_t_start, by = "row"),
    .combine = c
  ) %do% {
    sprintf(" {%f, %f} ", r$time, r$cases)
  } %>%
    stringr::str_c(collapse = ", \n")

  matrix_cases_at_t_start.string <- stringr::str_c(
    sprintf("double cases_at_t_start[%i][%i] = {\n", nrow(cases_at_t_start), 2),
    cases_at_t_start.string,
    " \n };"
  )

  cases_other.string <- foreach(
    r = iterators::iter(cases_covar, by = "row"),
    .combine = c
  ) %do% {
    sprintf(" {%f, %f} ", r$time, r$cases)
  } %>%
    stringr::str_c(collapse = ", \n")

  matrix_cases_other.string <- stringr::str_c(
    sprintf("double cases_other[%i][%i] = {\n", nrow(cases_covar), 2),
    cases_other.string,
    " \n };"
  )

  if (B0) {
    initalizeStates <- pomp::Csnippet("
  A     = nearbyint((1-sigma)/sigma  * 1/epsilon * cases_at_t_start[n_cases_start-1][1]/7 * 365 /(mu+gammaA));
  I     = nearbyint(1/epsilon * cases_at_t_start[n_cases_start-1][1]/7 * 365 /(mu+alpha+gammaI))  ;  // Steady state, DP says its correct.
  double R0[2] = {0,0};
  //FILE *fp;
  //fp = fopen('Ouuuutput.txt', 'w');

  double B_acc = 0;
  double rhoI = rhoA * XrhoI;
  double thetaA = thetaI * XthetaA;
  for(int i = 0; i < n_cases_start; i++){
    R0[0] +=                   cases_at_t_start[i][1]/epsilon  * exp((cases_at_t_start[i][0] - t_start)  * (rhoI+mu)); /* because t_i in past so t_ - t_0 negative */
    R0[1] += (1-sigma)/sigma * cases_at_t_start[i][1]/epsilon  * exp((cases_at_t_start[i][0] - t_start)  * (rhoA+mu));
    B_acc += (thetaA * (1-sigma)/sigma * cases_at_t_start[i][1]/epsilon + thetaI * cases_at_t_start[i][1]/epsilon) *
              (1 + lambdaR * pow(0.024, r)) * D * exp((cases_at_t_start[i][0] - t_start)  * mu_B);

    //fprintf(fp, '%f %f %f', cases_at_t_start[i][0] - t_start, cases_at_t_start[i][0], t_start);
  }
 //fclose(fp);

  B = B_acc;
  RI1   = nearbyint(R0[0]/3);
  RI2   = nearbyint(R0[0]/3);
  RI3   = nearbyint(R0[0]/3);
  RA1   = nearbyint(R0[1]/3);
  RA2   = nearbyint(R0[1]/3);
  RA3   = nearbyint(R0[1]/3);
  if (A + I + RI1 + RI2 + RI3 + RA1 + RA2 + RA3 >= H)
  {
    double R_tot = H - A - I - 100.0;
    if (R_tot <= 0)
    {
      I     = nearbyint(H - 100);
      A     = nearbyint(0);
      R_tot = nearbyint(0);
    }
    RI1   = nearbyint(sigma * R_tot/3.0);
    RI2   = nearbyint(sigma * R_tot/3.0);
    RI3   = nearbyint(sigma * R_tot/3.0);
    RA1   = nearbyint((1-sigma) * R_tot/3.0);
    RA2   = nearbyint((1-sigma) * R_tot/3.0);
    RA3   = nearbyint((1-sigma) * R_tot/3.0);
  }
  S   = nearbyint(H - A - I - RI1 - RI2 - RI3 - RA1 - RA2 - RA3);
  // B   = (I * thetaI/mu_B + A * thetaA/mu_B) * D * (1 + lambdaR * pow(B0, r)); // TODO custom initial conditions equivalent to the 'forcing' in the continous model
  B = B0;
  C   = 0;
  W   = 0;
  VSd = 0;
  VRI1d = 0;
  VRI2d = 0;
  VRI3d = 0;
  VRA1d = 0;
  VRA2d = 0;
  VRA3d = 0;
  VSdd = 0;
  VRI1dd = 0;
  VRI2dd = 0;
  VRI3dd = 0;
  VRA1dd = 0;
  VRA2dd = 0;
  VRA3dd = 0;
  VSd_alt = 0;
  VRI1d_alt = 0;
  VRI2d_alt = 0;
  VRI3d_alt = 0;
  VRA1d_alt = 0;
  VRA2d_alt = 0;
  VRA3d_alt = 0;
  VSdd_alt = 0;
  VRI1dd_alt = 0;
  VRI2dd_alt = 0;
  VRI3dd_alt = 0;
  VRA1dd_alt = 0;
  VRA2dd_alt = 0;
  VRA3dd_alt = 0;
   ")
  } else {
    initalizeStates <- pomp::Csnippet("
  A     = nearbyint((1-sigma)/sigma  * 1/epsilon * cases_at_t_start[n_cases_start-1][1]/7 * 365 /(mu+gammaA));
  I     = nearbyint(1/epsilon * cases_at_t_start[n_cases_start-1][1]/7 * 365 /(mu+alpha+gammaI))  ;  // Steady state, DP says its correct.
  double R0[2] = {0,0};
  //FILE *fp;
  //fp = fopen('Ouuuutput.txt', 'w');

  double B_acc = 0;
  double rhoI = rhoA * XrhoI;
  double thetaA = thetaI * XthetaA;
  for(int i = 0; i < n_cases_start; i++){
    R0[0] +=                   cases_at_t_start[i][1]/epsilon  * exp((cases_at_t_start[i][0] - t_start)  * (rhoI+mu)); /* because t_i in past so t_ - t_0 negative */
    R0[1] += (1-sigma)/sigma * cases_at_t_start[i][1]/epsilon  * exp((cases_at_t_start[i][0] - t_start)  * (rhoA+mu));
    B_acc += (thetaA * (1-sigma)/sigma * cases_at_t_start[i][1]/epsilon + thetaI * cases_at_t_start[i][1]/epsilon) *
              (1 + lambdaR * pow(0.024, r)) * D * exp((cases_at_t_start[i][0] - t_start)  * mu_B);

    //fprintf(fp, '%f %f %f', cases_at_t_start[i][0] - t_start, cases_at_t_start[i][0], t_start);
  }
 //fclose(fp);

  B = B_acc;
  RI1   = nearbyint(R0[0]/3);
  RI2   = nearbyint(R0[0]/3);
  RI3   = nearbyint(R0[0]/3);
  RA1   = nearbyint(R0[1]/3);
  RA2   = nearbyint(R0[1]/3);
  RA3   = nearbyint(R0[1]/3);
  if (A + I + RI1 + RI2 + RI3 + RA1 + RA2 + RA3 >= H)
  {
    double R_tot = H - A - I - 100.0;
    if (R_tot <= 0)
    {
      I     = nearbyint(H - 100);
      A     = nearbyint(0);
      R_tot = nearbyint(0);
    }
    RI1   = nearbyint(sigma * R_tot/3.0);
    RI2   = nearbyint(sigma * R_tot/3.0);
    RI3   = nearbyint(sigma * R_tot/3.0);
    RA1   = nearbyint((1-sigma) * R_tot/3.0);
    RA2   = nearbyint((1-sigma) * R_tot/3.0);
    RA3   = nearbyint((1-sigma) * R_tot/3.0);
  }
  S   = nearbyint(H - A - I - RI1 - RI2 - RI3 - RA1 - RA2 - RA3);
  B   = (I * thetaI/mu_B + A * thetaA/mu_B) * D * (1 + lambdaR * pow(0.024, r)); // TODO custom initial conditions equivalent to the 'forcing' in the continous model
  C   = 0;
  W   = 0;
  VSd = 0;
  VRI1d = 0;
  VRI2d = 0;
  VRI3d = 0;
  VRA1d = 0;
  VRA2d = 0;
  VRA3d = 0;
  VSdd = 0;
  VRI1dd = 0;
  VRI2dd = 0;
  VRI3dd = 0;
  VRA1dd = 0;
  VRA2dd = 0;
  VRA3dd = 0;
  VSd_alt = 0;
  VRI1d_alt = 0;
  VRI2d_alt = 0;
  VRI3d_alt = 0;
  VRA1d_alt = 0;
  VRA2d_alt = 0;
  VRA3d_alt = 0;
  VSdd_alt = 0;
  VRI1dd_alt = 0;
  VRI2dd_alt = 0;
  VRI3dd_alt = 0;
  VRA1dd_alt = 0;
  VRA2dd_alt = 0;
  VRA3dd_alt = 0;
   ")
  }

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
}  double eff_v_1d(double t_since_vacc, int scenario) {
  if (t_since_vacc < 1)
    return eff_v_2d(t_since_vacc, scenario);
  else
   return 0;
};"

  populations  <- unlist(purrr::flatten(MODEL3_INPUT_PARAMETERS["population"]))
  densities <- unlist(purrr::flatten(input_parameters["density"]))
  # param_proc_fixed['H'] <- populations[departement]
  # param_proc_fixed['D'] <- densities[departement]

  p1d_alt_year  <- unlist(purrr::flatten(input_parameters["p1d_alt_year"]))
  nb_doses_alt_year <- unlist(purrr::flatten(input_parameters["nb_doses_alt_year"]))
  t_vacc_start_alt  <- unlist(purrr::flatten(input_parameters["t_vacc_start_alt"]))
  t_vacc_end_alt <- unlist(purrr::flatten(input_parameters["t_vacc_end_alt"]))

  t_vacc_start_alt = dateToYears(as.Date(t_vacc_start_alt[departement]))
  t_vacc_end_alt   = dateToYears(as.Date(t_vacc_end_alt[departement]))
  r_v_alt_year = nb_doses_alt_year[departement]/(t_vacc_end_alt - t_vacc_start_alt)
  p1d_alt = p1d_alt_year[departement]

  cases_ext_mean <- cases_covar %>% dplyr::filter(time > t_start)
  cases_ext_mean <- mean(cases_ext_mean$cases)

  # rate of simulation in fractions of years
  dt_yrs <- delta.t

  # TODO: Put actual parameter values here.
  params <- DATA$params
  params['H'] <- populations[departement]
  params['D'] <- densities[departement]

  if (betaB_trend) {
    params['betaB_trend'] <- 0
  }

  if (B0) {
    params['B0'] <- 0.2
  }

  rain <- rain %>%
    dplyr::filter(time > (t_start - 0.01) & time < (t_end + 0.01)) %>%
    dplyr::select(time, rain_std) %>%
    dplyr::rename(rain = rain_std)

  cases_covar <- cases_covar %>%
    dplyr::filter(time > (t_start - 0.01) & time < (t_end + 0.01)) %>%
    dplyr::select(time, cases_other) %>%
    dplyr::rename(cases_covar_c = cases_other)

  covar <- dplyr::full_join(rain, cases_covar)

  if (B0) {
    pt <- pomp::parameter_trans(
      log = c(
        "betaB", "mu_B", "thetaI", "rhoA", "lambdaR", "r",
        "std_W", "k", "foi_add", "gammaA", "gammaI"
      ),
      logit = c(
        "sigma",
        "XthetaA",
        "XrhoI",
        "epsilon",
        "cas_def",
        "Rtot_0",
        "B0"
      )
    )
  } else {
    pt <- pomp::parameter_trans(
      log = c(
        "betaB", "mu_B", "thetaI", "rhoA", "lambdaR", "r",
        "std_W", "k", "foi_add", "gammaA", "gammaI"
      ),
      logit = c(
        "sigma",
        "XthetaA",
        "XrhoI",
        "epsilon",
        "cas_def",
        "Rtot_0"
      )
    )
  }


  sirb_cholera <- pomp::pomp(
    # set data
    data = cases %>%
      dplyr::filter(time > t_start & time < (t_end + 0.01)) %>%
      dplyr::select(time, cases) ,
    # time column
    times = "time",
    # initialization time
    t0 = t_start - dt_yrs,
    # paramter vector
    params = params,
    # process simulator
    rprocess = pomp::euler(step.fun = sirb.rproc, delta.t = dt_yrs),
    # measurement model simulator
    rmeasure =  rmeas,
    # measurement model density
    dmeasure = dmeas,
    # covariates
    covar = pomp::covariate_table(covar, times = 'time'),
    # names of state variables
    statenames = state_names,
    # names of accumulator variables to be re-initalized at each observation timestep
    # (C for cases, W for the white noise just for plotting)
    accumvars = c("C", "W"),
    # names of paramters
    paramnames = names(params),
    rinit = initalizeStates,
    partrans = pt,
    # global C definitions
    globals = stringr::str_c(
      sprintf("double t_vacc_start_alt = %f; double t_vacc_end_alt = %f;", t_vacc_start_alt, t_vacc_end_alt),
      sprintf("double r_v_alt = %f;", r_v_alt_year),
      sprintf("double p1d_alt = %f;", p1d_alt),
      sprintf("int n_cases_start = %i;",  nrow(cases_at_t_start)),
      sprintf("int n_cases_other = %i;",  nrow(cases_covar)),
      sprintf("double t_start = %f;",  t_start),
      sprintf("double t_end = %f;",  t_end),
      derivativeBacteria.c,
      matrix_cases_at_t_start.string,
      matrix_cases_other.string,
      eff_v.c,
      sep = " ")
  )

  sirb_cholera
}
