#' Build single department POMP object following model 3
#'
#' This code is a slight modification of the \code{\link{haiti3_dep}} function,
#' correcting minor flaws in the original code.
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

haiti3_dep_correct <- function(departement = 'Artibonite',
                               start_time = "2010-10-23",
                               delta.t = 1/365.25) {

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
  # cases_other_dept <- DATA$cases_other_dept

  t_start <- dateToYears(as.Date(start_time))

  # if (departement == 'Ouest' && start_time != '2014-03-01') {
  #   t_start <- dateToYears(as.Date('2017-06-10'))
  # } else {
  #   t_start <- dateToYears(as.Date(input_parameters$t_start))
  # }
  t_end <- dateToYears(as.Date(input_parameters$t_end))

  state_names <- c(
    "S", "I", "A", "R1", "R2", "R3",
    "VSd", "VR1d", "VR2d", "VR3d",
    "VSdd", "VR1dd", "VR2dd", "VR3dd",
    "VSd_alt", "VR1d_alt", "VR2d_alt", "VR3d_alt",  # For previous vaccination campagins
    "VSdd_alt", "VR1dd_alt", "VR2dd_alt", "VR3dd_alt",
    "B", "C"
  )

  dmeas <- pomp::Csnippet("
  if (ISNA(cases)) {
    lik = (give_log) ? 0 : 1;
    } else {
      if (t > 2018) {
        lik = dnbinom_mu(cases, k, epsilon * C * cas_def, give_log);
      } else {
        lik = dnbinom_mu(cases, k, epsilon * C, give_log) ;
      }
    }
    ")

  rmeas <- pomp::Csnippet("
  if (t > 2018) {
    cases = rnbinom_mu(k, epsilon * C * cas_def);
  } else {
    cases = rnbinom_mu(k, epsilon * C);
  }
  ")


  sirb.rproc <- pomp::Csnippet("

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
int scenario = cases_ext;
double t_eff, t_eff_alt;

double thetaA = thetaI * XthetaA;

if (t < 2018){
  /* We want estimate of instantaneous I and A. To do this, we divide the total
     number of observed cases by reporting rate to get total number of newly
     infected individuals for the week. Then, we convert the weekly cases into
     daily cases by dividing by 7 (assuming each day is approximately the same
     number of cases). Finally, individuals remain infected for an average of
     (1/gamma) days
  */
	mobility = (365 * cases_other / (7 * epsilon)) * (1 / (mu + alpha + gamma) + (1 - sigma) / (sigma * (mu + gamma)));
}
else
{
	mobility = (365 * cases_other / (7 * epsilon * cas_def)) * (1 / (mu + alpha + gamma) + (1 - sigma) / (sigma * (mu + gamma)));
}

// force of infection
foi = betaB * (B / (1 + B)) + foi_add * mobility;

if(std_W > 0.0) {
  dw = rgammawn(std_W, dt);    // white noise (extra-demographic stochasticity)
  foi_stoc = foi * dw/dt;      // apply stochasticity
} else {
  foi_stoc = foi;
}

if (t <= (t_vacc_end_alt + dt)){
	previous_vacc_campaign = TRUE;
	if (t >= t_vacc_start_alt && t <= (t_vacc_end_alt + dt)) {
    	r_v_wdn = (r_v_alt / (S + A + R1 + R2 + R3));
	}
	p1d = p1d_alt;
} else {
	previous_vacc_campaign = FALSE;
	if (t >= t_vacc_start && t <= (t_vacc_end + dt)) {
    	r_v_wdn = (r_v_year / (S + A + R1 + R2 + R3));
	}
	p1d = p1d_reg;
}

pdd = 1 - p1d;

// time in the vacc_eff referential. We assume different timing for 1d and 2d
t_eff =     t - (t_vacc_start + (t_vacc_end - t_vacc_start)/2);
t_eff_alt = t - (t_vacc_start_alt + (t_vacc_end_alt - t_vacc_start_alt)/2);

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
reulermultinom(5, S,     &rate[0],  dt, &dN[0]);
reulermultinom(3, I,     &rate[5],  dt, &dN[5]);
reulermultinom(4, A,     &rate[8],  dt, &dN[8]);
reulermultinom(4, R1,    &rate[12], dt, &dN[12]);
reulermultinom(4, R2,    &rate[16], dt, &dN[16]);
reulermultinom(4, R3,    &rate[20], dt, &dN[20]);

// Vaccinated 1 dose
reulermultinom(3,  VSd,   &rate[24], dt, &dN[24]);
reulermultinom(2, VR1d,   &rate[27], dt, &dN[27]);
reulermultinom(2, VR2d,   &rate[29], dt, &dN[29]);
reulermultinom(2, VR3d,   &rate[31], dt, &dN[31]);

// Vaccinated 2 doses
reulermultinom(3,  VSdd,   &rate[33], dt, &dN[33]);
reulermultinom(2, VR1dd,   &rate[36], dt, &dN[36]);
reulermultinom(2, VR2dd,   &rate[38], dt, &dN[38]);
reulermultinom(2, VR3dd,   &rate[40], dt, &dN[40]);

/* For the previous vaccination campain */

// Alt Vaccinated 1 dose
reulermultinom(3,  VSd_alt,   &rate[42], dt, &dN[42]);
reulermultinom(2, VR1d_alt,   &rate[27], dt, &dN[45]);
reulermultinom(2, VR2d_alt,   &rate[29], dt, &dN[47]);
reulermultinom(2, VR3d_alt,   &rate[31], dt, &dN[49]);

// Alt Vaccinated 2 doses
reulermultinom(3,  VSdd_alt,   &rate[45], dt, &dN[51]);
reulermultinom(2, VR1dd_alt,   &rate[27], dt, &dN[54]);
reulermultinom(2, VR2dd_alt,   &rate[29], dt, &dN[56]);
reulermultinom(2, VR3dd_alt,   &rate[31], dt, &dN[58]);

// bacteria as continous state variable
// implement Runge-Kutta integration assuming S, I, R, V* stay constant during dt
k1 = dt * fB(I, A, B, mu_B, thetaI, thetaA, lambdaR, rain, r, D);
k2 = dt * fB(I, A, B, mu_B, thetaI, thetaA, lambdaR, rain, r, D);
k3 = dt * fB(I, A, B, mu_B, thetaI, thetaA, lambdaR, rain, r, D);
k4 = dt * fB(I, A, B, mu_B, thetaI, thetaA, lambdaR, rain, r, D);

// bacteria increment
dB = (k1 + 2*k2 + 2*k3 + k4) / 6.0;

// Update States
I  += dN[0] + dN[24] + dN[33] + dN[42] + dN[51] - dN[7] - dN[6] - dN[5];            // S -> I, VSd -> I, VSdd -> I, VSd_alt -> I, VSdd_alt -> I, I -> R, I-> death, I -> death
A  += dN[1] + dN[25] + dN[34] + dN[43] + dN[52] - dN[9] - dN[10] - dN[11] - dN[8];  // S -> A, VSd -> A, VSdd -> A, VSd_alt -> A, VSdd_alt -> A, A -> R1, A -> VR1d, A -> VR1dd, natural death.
R1 += dN[7] + dN[9] - dN[12] - dN[14] - dN[15] - dN[13];  // I -> R1, A -> R1, R1 -> R2, R1 -> VR1d, R1 -> VR1dd, R1 -> death
R2 += dN[12] - dN[16] - dN[18] - dN[19] - dN[17]; // R1 -> R2, R2 -> R3, R2 -> VR2d, R2 -> VR2dd, R2 -> death
R3 += dN[16] - dN[20] - dN[22] - dN[23] - dN[21]; // R2 -> R3, R3 -> S , R3 -> VR3d, R3 -> VR3dd, R3 -> death

if (previous_vacc_campaign){
	VSd_alt    +=  dN[2];            // S -> VSd_alt
	VR1d_alt   +=  dN[14] + dN[10];  // R1 -> VR1d_alt, A -> VR1d_alt
	VR2d_alt   +=  dN[18];           // R2 -> VR2d_alt
	VR3d_alt   +=  dN[22];           // R3 -> VR3d_alt

	VSdd_alt   += dN[3];             // S -> VSdd_alt
	VR1dd_alt  += dN[15] + dN[11];   // R1 -> VR1dd_alt, A -> VR1dd_alt
	VR2dd_alt  += dN[19];            // R2 -> VR2dd_alt
	VR3dd_alt  += dN[23];            // R3 -> VR3dd_alt

} else {
	VSd    += dN[2];                 // S -> VSd
	VR1d   += dN[14] + dN[10];       // R1 -> VR1d, A -> VR1d
	VR2d   += dN[18];                // R2 -> VR2d
	VR3d   += dN[22];                // R3 -> VR3d

	VSdd   += dN[3];                 // S -> VSdd
	VR1dd  += dN[15] + dN[11];       // R1 -> VR1dd, A -> VR1dd
	VR2dd  += dN[19];                // R2 -> VR2dd
	VR3dd  += dN[23];                // R3 -> VR3dd
}

VSd   += dN[32] - dN[24] - dN[25] - dN[26]; // VR3d -> VSd, VSd -> I, VSd -> A, VSd -> death
VR1d  += - dN[27] - dN[28];                 // VR1d -> death, VR1d -> VR2d
VR2d  += dN[28] - dN[29] - dN[30];          // VR1d -> VR2d, VR2d -> death, VR2d -> VR3d
VR3d  += dN[30] - dN[31] - dN[32];          // VR2d -> VR3d, VR3d -> death, VR3d -> VSd

VSdd   +=  dN[41] - dN[33] - dN[34] - dN[35]; // VR3dd -> VSdd, VSdd -> I, VSdd -> A, VSdd -> death
VR1dd  += - dN[36] - dN[37];                  // VR1dd -> death, VR1dd -> VR2dd
VR2dd  +=  dN[37] - dN[38] - dN[39];          // VR1dd -> VR2dd, VR2dd -> death, VR2dd -> VR3dd
VR3dd  +=  dN[39] - dN[40] - dN[41];          // VR2dd -> VR3dd, VR3dd -> death, VR3dd -> VSdd

// previous vacccination campain
VSd_alt   += dN[50] - dN[42] - dN[43] - dN[44]; // VR3d_alt -> VSd_alt, VSd_alt -> I, VSd_alt -> A, VSd_alt -> death
VR1d_alt  += - dN[45] - dN[46];                 // death, VR1d_alt -> VR2d_alt
VR2d_alt  += dN[46] - dN[47] - dN[48];          // VR1d_alt -> VR2d_alt, death, VR2d_alt -> VR3d_alt
VR3d_alt  += dN[48] - dN[49] - dN[50];          // VR2d_alt -> VR3d_alt, death, VR3d_alt -> VSd_alt

VSdd_alt   +=  dN[59] - dN[51] - dN[52] - dN[53]; // VR3dd_alt -> VSdd_alt, VSdd_alt -> I, VSdd_alt -> A, VSdd_alt -> death
VR1dd_alt  += - dN[54] - dN[55];
VR2dd_alt  +=  dN[55] - dN[56] - dN[57] ;
VR3dd_alt  +=  dN[57] - dN[58] - dN[59];

// Rprintf(\"I: %f\\n\", I);

C   +=  dN[0] + dN[24] + dN[33] + dN[42] + dN[51]; // S -> I, VSd -> I, VSdd -> I, VSd_alt -> I, VSdd_alt -> I
B   += (((dB) < -B) ? (-B + 1.0e-3) : (dB)); // condition to ensure B>0

// susceptibles so as to match total population
S = nearbyint(H - I - A - R1 - R2 - R3 -
	VSd - VR1d - VR2d - VR3d -
	VSdd- VR1dd -VR2dd -VR3dd -
	VSd_alt - VR1d_alt - VR2d_alt - VR3d_alt -
	VSdd_alt - VR1dd_alt - VR2dd_alt - VR3dd_alt);
if (S < 0) {
  S = 0;
}
                         ")



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

  # cases_other.string <- foreach(
  #   r = iterators::iter(cases_covar, by = "row"),
  #   .combine = c
  # ) %do% {
  #   sprintf(" {%f, %f} ", r$time, r$cases_other)
  # } %>%
  #   stringr::str_c(collapse = ", \n")
  #
  # matrix_cases_other.string <- stringr::str_c(
  #   sprintf("double cases_other[%i][%i] = {\n", nrow(cases_covar), 2),
  #   cases_other.string,
  #   " \n };"
  # )

  initalizeStates <- pomp::Csnippet("
  A = nearbyint((1-sigma)/sigma * 1/epsilon * cases_at_t_start[n_cases_start-1][1]/7 * 365 /(mu+gamma));
  I = nearbyint(1/epsilon * cases_at_t_start[n_cases_start-1][1]/7 * 365 /(mu+alpha+gamma))  ;  // Steady state
  double R0[2] = {0,0};

  double thetaA = thetaI * XthetaA;
  for(int i = 0; i < n_cases_start; i++){
    R0[0] +=                   cases_at_t_start[i][1]/epsilon  * exp((cases_at_t_start[i][0] - t_start)  * (rho+mu)); /* because t_i in past so t_ - t_0 negative */
    R0[1] += (1-sigma)/sigma * cases_at_t_start[i][1]/epsilon  * exp((cases_at_t_start[i][0] - t_start)  * (rho+mu));
  }

  R1 = nearbyint((R0[0] + R0[1]) / 3);
  R2 = nearbyint((R0[0] + R0[1]) / 3);
  R3 = nearbyint((R0[0] + R0[1]) / 3);

  if (A + I + R1 + R2 + R3 >= H) {
    double R_tot = H - A - I - 100.0;
    if (R_tot <= 0)
    {
      I     = nearbyint(H - 100);
      A     = nearbyint(0);
      R_tot = nearbyint(0);
    }
    R1 = nearbyint(R_tot / 3);
    R2 = nearbyint(R_tot / 3);
    R3 = nearbyint(R_tot / 3);
  }

  S   = nearbyint(H - A - I - R1 - R2 - R3);
  B   = (I * thetaI/mu_B + A * thetaA/mu_B) * D * (1 + lambdaR * pow(B0, r));
  C   = 0;

  VSd  = 0;
  VR1d = 0;
  VR2d = 0;
  VR3d = 0;

  VSdd = 0;
  VR1dd = 0;
  VR2dd = 0;
  VR3dd = 0;

  VSd_alt = 0;
  VR1d_alt = 0;
  VR2d_alt = 0;
  VR3d_alt = 0;

  VSdd_alt = 0;
  VR1dd_alt = 0;
  VR2dd_alt = 0;
  VR3dd_alt = 0;
   ")

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

  populations  <- unlist(purrr::flatten(input_parameters["population"]))
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
  gamma <- params['gammaI']
  rho <- params['rhoA']
  params <- params[!names(params) %in% c('rhoA', 'XrhoI', 'gammaI', 'gammaA')]
  params['gamma'] <- gamma
  params['rho'] <- rho

  params['H'] <- populations[departement]
  params['D'] <- densities[departement]
  params['B0'] <- 0.024

  rain <- rain %>%
    dplyr::filter(time > (t_start - 0.01) & time < (t_end + 0.01)) %>%
    dplyr::select(time, rain_std) %>%
    dplyr::rename(rain = rain_std)

  cases_covar <- cases_covar %>%
    dplyr::filter(time > (t_start - 0.01) & time < (t_end + 0.01))

  covar <- dplyr::full_join(rain, cases_covar, by = 'time') %>%
    tidyr::fill(cases_other, .direction = c('updown'))

  pt <- pomp::parameter_trans(
    log = c(
      "betaB", "mu_B", "thetaI", "rho", "lambdaR", "r",
      "std_W", "k", "foi_add", "gamma"
    ),
    logit = c(
      "sigma",
      "XthetaA",
      "epsilon",
      "cas_def",
      "Rtot_0",
      "B0"
    )
  )

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
    # (C for cases)
    accumvars = c("C"),
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
      # sprintf("int n_cases_other = %i;",  nrow(cases_other_dept)),
      sprintf("double t_start = %f;",  t_start),
      sprintf("double t_end = %f;",  t_end),
      derivativeBacteria.c,
      matrix_cases_at_t_start.string,
      # matrix_cases_other.string,
      eff_v.c,
      sep = " ")
  )

  sirb_cholera
}
