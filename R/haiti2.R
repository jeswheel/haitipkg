#' Build SpatPomp object for Model 2
#'
#' Generate a \sQuote{spatPomp} object for fitting to Haiti cholera data. Based on
#' the model developed by Laura Matrajt et al. at the University of Florida.
#'
#' @importFrom pomp Csnippet
#' @param cutoff a numeric. Used to define the cutoff time for the model
#' @param region a character string. Specifies whether the time region is before
#' or after the cutoff time
#' @return An object of class \sQuote{spatPomp}.
#' @examples
#' c_epidemic <- haiti2(cutoff=2014.25, region="before")
#' @export


haiti2 <- function(cutoff=2020, region=c("before","after")){
  cholera_unit_statenames <- c(
    "S","E","I","A","R","RA","C",
    "VOD","EOD","IOD","AOD","ROD","RAOD",
    "VTD","ETD","ITD","ATD","RTD","RATD","W")
  cholera_IVPnames <- paste0(c("InitInfected"),1:10)
  cholera_RPnames <- c(
    "sigma","gammaE","k","gamma","Beta","Omega1",
    "Omega2","VE1","VE2","Delta","Mu","AlphaS","ps",
    "Sat","BetaW","sigmaSE","VR","WR","Psi","Rho")
  cholera_paramnames <- c(cholera_RPnames,cholera_IVPnames)

  # make components of spatPomp object
  cholera_globals <- Csnippet("
    const double Wmat[10][10] = {
    {0,0.158,0,0,0.0711,0,0.0124,0,0,0},
    {0.0294,0,0,0,0.00278,0.00776,0,0.00236,0,0},
    {0,0,0,0,0,0,0,0,0.0246,0},
    {0,0,0,0,0,0,0,0.00157,0.0796,0.00792},
    {0.0218,0.175,0,0,0,0.0285,0.00926,0,0,0},
    {0,0.0733,0,0,0.00860,0,0,0,0,0},
    {0.00105,0,0,0,0.00860,0,0,0,0,0},
    {0,0.0183,0,0.00264,0,0,0,0,0,0.0717},
    {0,0,0.0276,0.0402,0,0,0,0,0,0.0305},
    {0,0,0,0.00118,0,0,0,0.0167,0.0156,0}
    };

    const int pop[10] = {1727524,746236,468301,342525,1067177,393967,
    728807,4029705,774976,632601};

    const double pi = 3.14159265358979323846;
    const double redb = 0.001;
    const double redmu = 0.000001;

    const double Tmat[10][10] = {
    {0,81990000,4103000,9178000,147000000,20610000,151900000,367900000,11900000,20640000},
    {81990000,0,2148000,5622000,46620000,10990000,14720000,320400000,6660000,13840000},
    {4103000,2148000,0,4436000,1631000,574000,1220000,20090000,30720000,2761000},
    {9178000,5622000,4436000,0,2768000,950500,2175000,102000000,40010000,11470000},
    {147000000,46620000,1631000,2768000,0,75920000,37350000,70460000,4179000,5901000},
    {20610000,10990000,574000,950500,75920000,0,14240000,23410000,1454000,2055000},
    {151900000,14720000,1220000,2175000,37350000,14240000,0,59380000,3199000,4647000},
    {367900000,320400000,20090000,102000000,70460000,23410000,59380000,0,79810000,248600000},
    {11900000,6660000,30720000,40010000,4179000,1454000,3199000,79810000,0,10230000},
    {20640000,13840000,2761000,11470000,5901000,2055000,4647000,248600000,10230000,0}
    };

    double inflow(const double mat[10][10], int j, const double *A){
      int u;
      double ret = 0;
      for(u=0; u<U; u++){
        ret += mat[u][j] * A[u];
      }
      return ret;
    }

    double outflow(const double mat[10][10], int j, const double *A){
      int u;
      double ret = 0;
      for(u=0; u<U; u++){
        ret += mat[j][u] * A[j];
      }
      return ret;
    }
  ")

  cholera_rinit_before <- Csnippet("
    double *S = &S1;
    double *E = &E1;
    double *I = &I1;
    double *A = &A1;
    double *R = &R1;
    double *RA = &RA1;
    double *C = &C1;
    double *VOD = &VOD1;
    double *EOD = &EOD1;
    double *IOD = &IOD1;
    double *AOD = &AOD1;
    double *ROD = &ROD1;
    double *RAOD = &RAOD1;
    double *VTD = &VTD1;
    double *ETD = &ETD1;
    double *ITD = &ITD1;
    double *ATD = &ATD1;
    double *RTD = &RTD1;
    double *RATD = &RATD1;
    double *W = &W1;
    const double *InitInfected = &InitInfected1;

    double m;
    int u;

    for (u = 0; u < U; u++) {
      S[u] = pop[u] - InitInfected[u];
      I[u] = InitInfected[u];
      E[u] = 0;
      A[u] = 0;
      R[u] = 0;
      RA[u] = 0;

      VOD[u] = 0;
      IOD[u] = 0;
      EOD[u] = 0;
      AOD[u] = 0;
      ROD[u] = 0;
      RAOD[u] = 0;

      VTD[u] = 0;
      ITD[u] = 0;
      ETD[u] = 0;
      ATD[u] = 0;
      RTD[u] = 0;
      RATD[u] = 0;

      C[u] = 0;

      W[u]=0;
    }
  ")

  cholera_rinit_after <- Csnippet("
    double *S = &S1;
    double *E = &E1;
    double *I = &I1;
    double *A = &A1;
    double *R = &R1;
    double *RA = &RA1;
    double *C = &C1;
    double *VOD = &VOD1;
    double *EOD = &EOD1;
    double *IOD = &IOD1;
    double *AOD = &AOD1;
    double *ROD = &ROD1;
    double *RAOD = &RAOD1;
    double *VTD = &VTD1;
    double *ETD = &ETD1;
    double *ITD = &ITD1;
    double *ATD = &ATD1;
    double *RTD = &RTD1;
    double *RATD = &RATD1;
    double *W = &W1;
    const double *InitInfected = &InitInfected1;

    double m;
    int u;

    for (u = 0; u < U; u++) {
      I[u] = InitInfected[u];
      A[u] = nearbyint((1-k)/k * I[u]);
      R[u] = nearbyint(k*0.25*(pop[u]-InitInfected[u]-(I[u]+A[u])));
      RA[u] = nearbyint((1-k)*0.25*(pop[u]-InitInfected[u]-(I[u]+A[u])));
      E[u] = 0;
      S[u] = nearbyint(0.75*(pop[u]-InitInfected[u]-(I[u]+A[u])));

      VOD[u] = 0;
      IOD[u] = 0;
      EOD[u] = 0;
      AOD[u] = 0;
      ROD[u] = 0;
      RAOD[u] = 0;

      VTD[u] = 0;
      ITD[u] = 0;
      ETD[u] = 0;
      ATD[u] = 0;
      RTD[u] = 0;
      RATD[u] = 0;

      C[u] = 0;

      W[u]=0;
    }
  ")

  cholera_rprocess <- Csnippet("
    double *S = &S1;
    double *E = &E1;
    double *I = &I1;
    double *A = &A1;
    double *R = &R1;
    double *RA = &RA1;
    double *C = &C1;
    double *VOD = &VOD1;
    double *EOD = &EOD1;
    double *IOD = &IOD1;
    double *AOD = &AOD1;
    double *ROD = &ROD1;
    double *RAOD = &RAOD1;
    double *VTD = &VTD1;
    double *ETD = &ETD1;
    double *ITD = &ITD1;
    double *ATD = &ATD1;
    double *RTD = &RTD1;
    double *RATD = &RATD1;
    double *W = &W1;

    double rate[11], trans[29];

    int u;
    double foi, dw;

    for (u = 0 ; u < U ; u++) {

      // expected force of infection
      foi = Beta*(I[u]+IOD[u]+ITD[u]) + redb*Beta*(A[u]+AOD[u]+ATD[u]);

      // water effect
      foi += 0.5*(1 + AlphaS*cos(2*pi*t/ps))*(BetaW*W[u])/(Sat+W[u]);

      // white noise (extrademographic stochasticity)
      dw = rgammawn(sigmaSE,dt);

      // rates
      rate[0] = foi*dw/dt;  // infection rate - no dose
      rate[1] = gammaE; // latent period
      rate[2] = gamma; // infectious period
      rate[3] = sigma; // rate of natural immunity
      rate[4] = Omega1; // rate of lose of vaccine immunity - one dose
      rate[5] = Omega2; // rate of lost of vaccine immunity - two dose
      rate[6] = (1-VE1)*foi*dw/dt;  // infection rate - one dose
      rate[7] = (1-VE2)*foi*dw/dt;  // infection rate - two dose
      rate[8] = Delta; // decay of cholera in water
      rate[9] = Mu; // symptomatic shedding rate
      rate[10] = Mu*redmu; // asymptomatic shedding rate

      // transitions between classes
      reulermultinom(1,S[u],&rate[0],dt,&trans[0]); // S to E
      reulermultinom(1,E[u],&rate[1],dt,&trans[1]); // E to (I,A)
      reulermultinom(1,I[u],&rate[2],dt,&trans[2]); // I to R
      reulermultinom(1,A[u],&rate[2],dt,&trans[3]); // A to RA
      reulermultinom(1,R[u],&rate[3],dt,&trans[4]); // R to S
      reulermultinom(1,RA[u],&rate[3],dt,&trans[5]); // RA to S

      reulermultinom(1,VOD[u],&rate[4],dt,&trans[6]); // V1 to S
      reulermultinom(1,VTD[u],&rate[5],dt,&trans[7]); // V2 to S

      reulermultinom(1,VOD[u],&rate[6],dt,&trans[8]); // V1 to E1
      reulermultinom(1,EOD[u],&rate[1],dt,&trans[9]); // E1 to (I1,A1)
      reulermultinom(1,IOD[u],&rate[2],dt,&trans[10]); // I1 to R1
      reulermultinom(1,AOD[u],&rate[2],dt,&trans[11]); // A1 to RA1
      reulermultinom(1,ROD[u],&rate[3],dt,&trans[12]); // R1 to S1
      reulermultinom(1,RAOD[u],&rate[3],dt,&trans[13]); // RA1 to S1

      reulermultinom(1,VTD[u],&rate[7],dt,&trans[14]); // V2 to E2
      reulermultinom(1,ETD[u],&rate[1],dt,&trans[15]); // E2 to (I2,A2)
      reulermultinom(1,ITD[u],&rate[2],dt,&trans[16]); // I2 to R2
      reulermultinom(1,ATD[u],&rate[2],dt,&trans[17]); // A2 to RA2
      reulermultinom(1,RTD[u],&rate[3],dt,&trans[18]); // R2 to S2
      reulermultinom(1,RATD[u],&rate[3],dt,&trans[19]); // RA2 to S2

      trans[20] = (Mu*(I[u]+IOD[u]+ITD[u]) + redmu*Mu*(A[u]+AOD[u]+ATD[u]))*dt; // (I,A) to W
      trans[21] = Delta*W[u]*dt; // water decay

      trans[22] = nearbyint(dt*VR*(inflow(Tmat,u,S) - outflow(Tmat,u,S))); // travel for S
      trans[23] = nearbyint(dt*VR*(inflow(Tmat,u,E) - outflow(Tmat,u,E))); // travel for E
      trans[24] = nearbyint(dt*VR*(inflow(Tmat,u,I) - outflow(Tmat,u,I))); // travel for I
      trans[25] = nearbyint(dt*VR*(inflow(Tmat,u,A) - outflow(Tmat,u,A))); // travel for A
      trans[26] = nearbyint(dt*VR*(inflow(Tmat,u,R) - outflow(Tmat,u,R))); // travel for R
      trans[27] = nearbyint(dt*VR*(inflow(Tmat,u,RA) - outflow(Tmat,u,RA))); // travel for RA

      trans[28] = nearbyint(dt*WR*(inflow(Wmat,u,W) - outflow(Wmat,u,W))); // water movement

      // make transitions
      S[u] += trans[4] + trans[5] + trans[6] + trans[7] + trans[22] - trans[0];
      E[u] += trans[0] + trans[23] - trans[1];
      I[u] += nearbyint(k*trans[1]) + trans[24] - trans[2];
      A[u] += nearbyint((1-k)*trans[1]) + trans[25] - trans[3];
      R[u] += trans[2] + trans[26] - trans[4];
      RA[u] += trans[3] + trans[27] - trans[5];

      VOD[u] += trans[12] + trans[13] - trans[6] - trans[8];
      EOD[u] += trans[8] - trans[9];
      IOD[u] += nearbyint(k*trans[9]) - trans[10];
      AOD[u] += nearbyint((1-k)*trans[9]) - trans[11];
      ROD[u] += trans[10] - trans[12];
      RAOD[u] += trans[11] - trans[13];

      VTD[u] += trans[18] + trans[19] - trans[7] - trans[14];
      ETD[u] += trans[14] - trans[15];
      ITD[u] += nearbyint(k*trans[15]) - trans[16];
      ATD[u] += nearbyint((1-k)*trans[15]) - trans[17];
      RTD[u] += trans[16] - trans[18];
      RATD[u] += trans[17] - trans[19];

      C[u] += trans[2] + trans[10] + trans[16];

      W[u] += trans[20] + trans[28] - trans[21];

      // printf('%f ',trans[33]); // for debugging
    }
  ")

  cholera_rmeasure <- Csnippet("
    const double *C = &C1;
    double *cases = &cases1;
    double m,v;
    double tol = 1.0e-300;
    int u;
    for (u = 0; u < U; u++) {
      m = (C[u]+tol)*Rho;
      v = m*(1.0-Rho + Psi*Psi*m);
      cases[u] = rnorm(m,sqrt(v)+tol);
      if (cases[u] > 0.0) {
        cases[u] = nearbyint(cases[u]);
      } else {
        cases[u] = 0.0;
      }
    }
  ")

  cholera_dmeasure <- Csnippet("
    const double *C = &C1;
    const double *cases = &cases1;
    double m,v;
    double tol = 1e-300;
    double mytol = 1e-5;
    int u;
    lik = 0;
    for (u = 0; u < U; u++) {
      m = Rho*(C[u]+mytol);
      v = m*(1.0-Rho+Psi*Psi*m);
      // C < 0 can happen in bootstrap methods such as bootgirf
      if (C < 0) {lik += log(tol);} else {
        if (cases[u] > tol) {
          lik += log(pnorm(cases[u]+0.5,m,sqrt(v)+tol,1,0)-
            pnorm(cases[u]-0.5,m,sqrt(v)+tol,1,0)+tol);
        } else {
            lik += log(pnorm(cases[u]+0.5,m,sqrt(v)+tol,1,0)+tol);
        }
      }
    }
    if(!give_log) lik = (lik > log(tol)) ? exp(lik) : tol;
  ")

  cholera_dunit_measure <- Csnippet("
    double mytol = 1e-5;
    double m = Rho*(C+mytol);
    double v = m*(1.0-Rho+Psi*Psi*m);
    double tol = 1e-300;
    // C < 0 can happen in bootstrap methods such as bootgirf
    if (C < 0) {lik = 0;} else {
      if (cases > tol) {
        lik = pnorm(cases+0.5,m,sqrt(v)+tol,1,0)-
          pnorm(cases-0.5,m,sqrt(v)+tol,1,0)+tol;
      } else {
        lik = pnorm(cases+0.5,m,sqrt(v)+tol,1,0)+tol;
      }
    }
    if(give_log) lik = log(lik);
  ")

  cholera_eunit_measure <- Csnippet("
    ey = Rho*C;
  ")

  cholera_vunit_measure <- Csnippet("
    double m = Rho*C;
    vc = m*(1.0-Rho + Psi*Psi*m);
  ")

  cholera_munit_measure <- Csnippet("
    double binomial_var;
    double m;
    double mytol = 1e-5;
    m = Rho*(C+mytol);
    binomial_var = Rho*(1-Rho)*C;
    if(vc > binomial_var) {
      M_Psi = sqrt(vc - binomial_var)/m;
    }
  ")

  cholera_skel <- Csnippet('
    const double *S = &S1;
    const double *E = &E1;
    const double *I = &I1;
    const double *A = &A1;
    const double *R = &R1;
    const double *RA = &RA1;
    const double *C = &C1;
    const double *VOD = &VOD1;
    const double *EOD = &EOD1;
    const double *IOD = &IOD1;
    const double *AOD = &AOD1;
    const double *ROD = &ROD1;
    const double *RAOD = &RAOD1;
    const double *VTD = &VTD1;
    const double *ETD = &ETD1;
    const double *ITD = &ITD1;
    const double *ATD = &ATD1;
    const double *RTD = &RTD1;
    const double *RATD = &RATD1;
    const double *W = &W1;

    double *DS = &DS1;
    double *DE = &DE1;
    double *DI = &DI1;
    double *DA = &DA1;
    double *DR = &DR1;
    double *DRA = &DRA1;
    double *DC = &DC1;
    double *DVOD = &DVOD1;
    double *DEOD = &DEOD1;
    double *DIOD = &DIOD1;
    double *DAOD = &DAOD1;
    double *DROD = &DROD1;
    double *DRAOD = &DRAOD1;
    double *DVTD = &DVTD1;
    double *DETD = &DETD1;
    double *DITD = &DITD1;
    double *DATD = &DATD1;
    double *DRTD = &DRTD1;
    double *DRATD = &DRATD1;
    double *DW = &DW1;

    double trans[29];

    int u;
    double foi;

    for (u = 0 ; u < U ; u++) {

      // expected force of infection
      foi = Beta*(I[u]+IOD[u]+ITD[u]) + redb*Beta*(A[u]+AOD[u]+ATD[u]);

      // water effect
      foi += 0.5*(1 + AlphaS*cos(2*pi*t/ps))*(BetaW*W[u])/(Sat+W[u]);

      // transitions between classes
      trans[0] = foi*S[u]; // S to E
      trans[1] = gammaE*E[u]; // E to (I,A)
      trans[2] = gamma*I[u]; // I to R
      trans[3] = gamma*A[u]; // A to RA
      trans[4] = sigma*R[u]; // R to S
      trans[5] = sigma*RA[u]; // RA to S

      trans[6] = Omega1*VOD[u]; // V1 to S
      trans[7] = Omega2*VTD[u]; // V2 to S

      trans[8] = (1-VE1)*foi*VOD[u]; // V1 to E1
      trans[9] = gammaE*EOD[u]; // E1 to (I1,A1)
      trans[10] = gamma*IOD[u]; // I1 to R1
      trans[11] = gamma*AOD[u]; // A1 to RA1
      trans[12] = sigma*ROD[u]; // R1 to S1
      trans[13] = sigma*RAOD[u]; // RA1 to S1

      trans[14] = (1-VE2)*foi*VTD[u]; // V2 to E2
      trans[15] = gammaE*ETD[u]; // E2 to (I2,A2)
      trans[16] = gamma*ITD[u]; // I2 to R2
      trans[17] = gamma*ATD[u]; // A2 to RA2
      trans[18] = sigma*RTD[u]; // R2 to S2
      trans[19] = sigma*RATD[u]; // RA2 to S2

      trans[20] = Mu*(I[u]+IOD[u]+ITD[u]) + redmu*Mu*(A[u]+AOD[u]+ATD[u]); // (I,A) to W
      trans[21] = Delta*W[u]; // water decay

      trans[22] = VR*(inflow(Tmat,u,S) - outflow(Tmat,u,S)); // travel for S
      trans[23] = VR*(inflow(Tmat,u,E) - outflow(Tmat,u,E)); // travel for E
      trans[24] = VR*(inflow(Tmat,u,I) - outflow(Tmat,u,I)); // travel for I
      trans[25] = VR*(inflow(Tmat,u,A) - outflow(Tmat,u,A)); // travel for A
      trans[26] = VR*(inflow(Tmat,u,R) - outflow(Tmat,u,R)); // travel for R
      trans[27] = VR*(inflow(Tmat,u,RA) - outflow(Tmat,u,RA)); // travel for RA

      trans[28] = WR*(inflow(Wmat,u,W) - outflow(Wmat,u,W)); // water movement

      // make transitions
      DS[u] += trans[4] + trans[5] + trans[6] + trans[7] + trans[22] - trans[0];
      DE[u] += trans[0] + trans[23] - trans[1];
      DI[u] += k*trans[1] + trans[24] - trans[2];
      DA[u] += (1-k)*trans[1] + trans[25] - trans[3];
      DR[u] += trans[2] + trans[26] - trans[4];
      DRA[u] += trans[3] + trans[27] - trans[5];

      DVOD[u] += trans[12] + trans[13] - trans[6] - trans[8];
      DEOD[u] += trans[8] - trans[9];
      DIOD[u] += k*trans[9] - trans[10];
      DAOD[u] += (1-k)*trans[9] - trans[11];
      DROD[u] += trans[10] - trans[12];
      DRAOD[u] += trans[11] - trans[13];

      DVTD[u] += trans[18] + trans[19] - trans[7] - trans[14];
      DETD[u] += trans[14] - trans[15];
      DITD[u] += k*trans[15] - trans[16];
      DATD[u] += (1-k)*trans[15] - trans[17];
      DRTD[u] += trans[16] - trans[18];
      DRATD[u] += trans[17] - trans[19];

      DC[u] += trans[2] + trans[10] + trans[16];

      DW[u] += trans[20] + trans[28] - trans[21];
    }
  ')

  cholera_partrans <- pomp::parameter_trans(
    log=c("Mu","Beta","BetaW")
  )

  haiti <- haiti2_data()
  if (region == "before"){
    haiti <- haiti[haiti$year <= cutoff,]
    cholera_rinit <- cholera_rinit_before
  } else {
    haiti <- haiti[haiti$year > cutoff,]
    cholera_rinit <- cholera_rinit_after
  }

  ret <- spatPomp::spatPomp(
    data = haiti,
    units = "department",
    times = "year",
    t0 = min(haiti$year)-1/52,
    unit_statenames = cholera_unit_statenames,
    rprocess=pomp::euler(cholera_rprocess, delta.t=2/365),
    skeleton=pomp::vectorfield(cholera_skel),
    unit_accumvars = c("C"),
    paramnames=cholera_paramnames,
    partrans=cholera_partrans,
    globals=cholera_globals,
    rinit=cholera_rinit,
    dmeasure=cholera_dmeasure,
    eunit_measure=cholera_eunit_measure,
    munit_measure=cholera_munit_measure,
    vunit_measure=cholera_vunit_measure,
    rmeasure=cholera_rmeasure,
    dunit_measure=cholera_dunit_measure
  )

  return(ret)
}
