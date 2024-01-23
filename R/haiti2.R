#' Build SpatPomp object for Model 2
#'
#' Generate a \sQuote{spatPomp} object for fitting to Haiti cholera data. Based on
#' the model developed by Laura Matrajt et al. at the University of Florida.
#'
#' Note that VR is \eqn{v_{rate}} in the manuscript. The units were calculated
#' using the fact that VR multiplies the travel matrix Tmat, which has units
#' in \eqn{km^{-2}}.
#'
#' @importFrom pomp Csnippet
#' @param cutoff a numeric. Used to define the cutoff time for the model
#' @param region a character string. Specifies whether the time region is before
#' or after the cutoff time.
#' @param measure a character string. Specifies the measurement model used. Valid
#' entries are "linear" and "log".
#' @return An object of class \sQuote{spatPomp}.
#' @examples
#' c_epidemic <- haiti2(cutoff=2014.25, region="before")
#' @export


haiti2 <- function(cutoff=2014.161, region="before", measure="linear"){
  cholera_unit_statenames <- c("S","E","I","A","R","RA","C",
                               "VOD","EOD","IOD","AOD","ROD","RAOD",
                               "VTD","ETD","ITD","ATD","RTD","RATD",
                               "VODu5","EODu5","IODu5","AODu5","RODu5","RAODu5",
                               "VTDu5","ETDu5","ITDu5","ATDu5","RTDu5","RATDu5",
                               "W")
  cholera_IVPnames <- paste0(c("InitInfected"),1:10)
  cholera_RPnames <- c("sigma","gammaE","k","gamma","Beta","Omega1",
                       "Omega2","VE1","VE2","Delta","Mu","AlphaS","ps",
                       "Sat","BetaW","sigmaSE","VR","WR","Psi","Rho",
                       "scenario","f","v","phase", "redb", "redmu")
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
    double p_u5 = 0.118;

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
    double *VODu5 = &VODu51;
    double *EODu5 = &EODu51;
    double *IODu5 = &IODu51;
    double *AODu5 = &AODu51;
    double *RODu5 = &RODu51;
    double *RAODu5 = &RAODu51;
    double *VTDu5 = &VTDu51;
    double *ETDu5 = &ETDu51;
    double *ITDu5 = &ITDu51;
    double *ATDu5 = &ATDu51;
    double *RTDu5 = &RTDu51;
    double *RATDu5 = &RATDu51;
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

      VODu5[u] = 0;
      IODu5[u] = 0;
      EODu5[u] = 0;
      AODu5[u] = 0;
      RODu5[u] = 0;
      RAODu5[u] = 0;

      VTDu5[u] = 0;
      ITDu5[u] = 0;
      ETDu5[u] = 0;
      ATDu5[u] = 0;
      RTDu5[u] = 0;
      RATDu5[u] = 0;

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
    double *VODu5 = &VODu51;
    double *EODu5 = &EODu51;
    double *IODu5 = &IODu51;
    double *AODu5 = &AODu51;
    double *RODu5 = &RODu51;
    double *RAODu5 = &RAODu51;
    double *VTDu5 = &VTDu51;
    double *ETDu5 = &ETDu51;
    double *ITDu5 = &ITDu51;
    double *ATDu5 = &ATDu51;
    double *RTDu5 = &RTDu51;
    double *RATDu5 = &RATDu51;
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

      VODu5[u] = 0;
      IODu5[u] = 0;
      EODu5[u] = 0;
      AODu5[u] = 0;
      RODu5[u] = 0;
      RAODu5[u] = 0;

      VTDu5[u] = 0;
      ITDu5[u] = 0;
      ETDu5[u] = 0;
      ATDu5[u] = 0;
      RTDu5[u] = 0;
      RATDu5[u] = 0;

      C[u] = 0;

      W[u]=0;
    }
  ")

  cholera_dmeasure <- Csnippet("
    double *C = &C1;
    double *cases = &cases1;
    double m;
    int u;
    lik = 0;
    for (u = 0; u < U; u++) {
      if(ISNA(cases[u])){
        lik += 0;
      } else {
        m = Rho*C[u];
        lik += dnorm(cases[u],m,v,1);
      }
    }
    if(!give_log) lik = exp(lik);
  ")

  cholera_log_dmeasure <- Csnippet("
    double *C = &C1;
    double *cases = &cases1;
    double m;
    int u;
    lik = 0;
    for (u = 0; u < U; u++) {
      if(ISNA(cases[u])){
        lik += 0;
      } else {
        m = log(Rho*C[u] + 1);
        lik += dnorm(log(cases[u] + 1),m,v,1) - log(cases[u] + 1);
      }
    }
    if(!give_log) lik = exp(lik);
  ")

  cholera_skel <- Csnippet('
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
    double *VODu5 = &VODu51;
    double *EODu5 = &EODu51;
    double *IODu5 = &IODu51;
    double *AODu5 = &AODu51;
    double *RODu5 = &RODu51;
    double *RAODu5 = &RAODu51;
    double *VTDu5 = &VTDu51;
    double *ETDu5 = &ETDu51;
    double *ITDu5 = &ITDu51;
    double *ATDu5 = &ATDu51;
    double *RTDu5 = &RTDu51;
    double *RATDu5 = &RATDu51;
    double *W = &W1;

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
    double *DVODu5 = &DVODu51;
    double *DEODu5 = &DEODu51;
    double *DIODu5 = &DIODu51;
    double *DAODu5 = &DAODu51;
    double *DRODu5 = &DRODu51;
    double *DRAODu5 = &DRAODu51;
    double *DVTDu5 = &DVTDu51;
    double *DETDu5 = &DETDu51;
    double *DITDu5 = &DITDu51;
    double *DATDu5 = &DATDu51;
    double *DRTDu5 = &DRTDu51;
    double *DRATDu5 = &DRATDu51;
    double *DW = &DW1;

    double trans[67];
    double vacc[2]; //vacc[0] is 2 dose rate, vacc[1] is 1 dose rate

    int u;
    double foi;

    for (u = 0 ; u < U ; u++) {

      if (scenario==1 && t>2019.03 && t<2021.03 && (u==0 || u==1)){
        // 2-dep
        vacc[0] = 0.7*pop[u]/2;
        vacc[1] = 0.1*pop[u]/2;
      }
      else if (scenario==2 && t>2019.03 && t<2021.03 && (u==0||u==1||u==7)){
      // 3-dep
        vacc[0] = 0.7*pop[u]/2;
        vacc[1] = 0.1*pop[u]/2;
      }
      else if (scenario==3 && t>2019.03 && t<2024.03){
        // slow national
        vacc[0] = 0.7*pop[u]/5;
        vacc[1] = 0.1*pop[u]/5;
      }
      else if (scenario==4 && t>2019.03 && t<2021.03){
        // fast national
        vacc[0] = 0.7*pop[u]/2;
        vacc[1] = 0.1*pop[u]/2;
      }
      else if (scenario==5 && t>2019.03 && t<2021.03){
        // fast high coverage national
        vacc[0] = 0.95*pop[u]/2;
        vacc[1] = 0.0167*pop[u]/2;
      }
      else {vacc[0]=0; vacc[1]=0;}

      // expected force of infection
      foi = Beta*(I[u]+IOD[u]+ITD[u]) + redb*Beta*(A[u]+AOD[u]+ATD[u]);

      // water effect
      foi += 0.5*(1 + AlphaS*cos(2*pi*t/ps + phase))*(BetaW*W[u])/(Sat+W[u]);

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

      trans[20] = Mu*(I[u]+IOD[u]+ITD[u]+IODu5[u]+ITDu5[u]) +
                  redmu*Mu*(A[u]+AOD[u]+ATD[u]+AODu5[u]+ATDu5[u]); // (I,A) to W
      trans[21] = Delta*W[u]; // water decay

      trans[22] = VR*(inflow(Tmat,u,S) - outflow(Tmat,u,S)); // travel for S
      trans[23] = VR*(inflow(Tmat,u,E) - outflow(Tmat,u,E)); // travel for E
      trans[24] = VR*(inflow(Tmat,u,I) - outflow(Tmat,u,I)); // travel for I
      trans[25] = VR*(inflow(Tmat,u,A) - outflow(Tmat,u,A)); // travel for A
      trans[26] = VR*(inflow(Tmat,u,R) - outflow(Tmat,u,R)); // travel for R
      trans[27] = VR*(inflow(Tmat,u,RA) - outflow(Tmat,u,RA)); // travel for RA

      trans[28] = WR*(inflow(Wmat,u,W) - outflow(Wmat,u,W)); // water movement

      trans[29] = Omega1*VODu5[u]; // V1 to S (u5)
      trans[30] = Omega2*VTDu5[u]; // V2 to S (u5)

      trans[31] = (1-VE1)*foi*VODu5[u]; // V1 to E1 (u5)
      trans[32] = gammaE*EODu5[u]; // E1 to (I1,A1) (u5)
      trans[33] = gamma*IODu5[u]; // I1 to R1 (u5)
      trans[34] = gamma*AODu5[u]; // A1 to RA1 (u5)
      trans[35] = sigma*RODu5[u]; // R1 to S1 (u5)
      trans[36] = sigma*RAODu5[u]; // RA1 to S1 (u5)

      trans[37] = (1-VE2)*foi*VTDu5[u]; // V2 to E2 (u5)
      trans[38] = gammaE*ETDu5[u]; // E2 to (I2,A2) (u5)
      trans[39] = gamma*ITDu5[u]; // I2 to R2 (u5)
      trans[40] = gamma*ATDu5[u]; // A2 to RA2 (u5)
      trans[41] = sigma*RTDu5[u]; // R2 to S2 (u5)
      trans[42] = sigma*RATDu5[u]; // RA2 to S2 (u5)

      trans[43] = VR*(inflow(Tmat,u,VOD) - outflow(Tmat,u,VOD)); // travel for VOD
      trans[44] = VR*(inflow(Tmat,u,EOD) - outflow(Tmat,u,EOD)); // travel for EOD
      trans[45] = VR*(inflow(Tmat,u,IOD) - outflow(Tmat,u,IOD)); // travel for IOD
      trans[46] = VR*(inflow(Tmat,u,AOD) - outflow(Tmat,u,AOD)); // travel for AOD
      trans[47] = VR*(inflow(Tmat,u,ROD) - outflow(Tmat,u,ROD)); // travel for ROD
      trans[48] = VR*(inflow(Tmat,u,RAOD) - outflow(Tmat,u,RAOD)); // travel for RAOD

      trans[49] = VR*(inflow(Tmat,u,VTD) - outflow(Tmat,u,VTD)); // travel for VTD
      trans[50] = VR*(inflow(Tmat,u,ETD) - outflow(Tmat,u,ETD)); // travel for ETD
      trans[51] = VR*(inflow(Tmat,u,ITD) - outflow(Tmat,u,ITD)); // travel for ITD
      trans[52] = VR*(inflow(Tmat,u,ATD) - outflow(Tmat,u,ATD)); // travel for ATD
      trans[53] = VR*(inflow(Tmat,u,RTD) - outflow(Tmat,u,RTD)); // travel for RTD
      trans[54] = VR*(inflow(Tmat,u,RATD) - outflow(Tmat,u,RATD)); // travel for RATD

      trans[55] = VR*(inflow(Tmat,u,VODu5) - outflow(Tmat,u,VODu5)); // travel for VODu5
      trans[56] = VR*(inflow(Tmat,u,EODu5) - outflow(Tmat,u,EODu5)); // travel for EODu5
      trans[57] = VR*(inflow(Tmat,u,IODu5) - outflow(Tmat,u,IODu5)); // travel for IODu5
      trans[58] = VR*(inflow(Tmat,u,AODu5) - outflow(Tmat,u,AODu5)); // travel for AODu5
      trans[59] = VR*(inflow(Tmat,u,RODu5) - outflow(Tmat,u,RODu5)); // travel for RODu5
      trans[60] = VR*(inflow(Tmat,u,RAODu5) - outflow(Tmat,u,RAODu5)); // travel for RAODu5

      trans[61] = VR*(inflow(Tmat,u,VTDu5) - outflow(Tmat,u,VTDu5)); // travel for VTDu5
      trans[62] = VR*(inflow(Tmat,u,ETDu5) - outflow(Tmat,u,ETDu5)); // travel for ETDu5
      trans[63] = VR*(inflow(Tmat,u,ITDu5) - outflow(Tmat,u,ITDu5)); // travel for ITDu5
      trans[64] = VR*(inflow(Tmat,u,ATDu5) - outflow(Tmat,u,ATDu5)); // travel for ATDu5
      trans[65] = VR*(inflow(Tmat,u,RTDu5) - outflow(Tmat,u,RTDu5)); // travel for RTDu5
      trans[66] = VR*(inflow(Tmat,u,RATDu5) - outflow(Tmat,u,RATDu5)); // travel for RATDu5


      // make transitions
      DS[u] = trans[4] + trans[5] + trans[6] + trans[7] + trans[22] +
              trans[29] + trans[30] - trans[0] - vacc[0] - vacc[1];
      DE[u] = trans[0] + trans[23] - trans[1];
      DI[u] = k*trans[1] + trans[24] - trans[2];
      DA[u] = (1-k)*trans[1] + trans[25] - trans[3];
      DR[u] = trans[2] + trans[26] - trans[4];
      DRA[u] = trans[3] + trans[27] - trans[5];

      DVOD[u] = trans[12] + trans[13] + trans[43] - trans[6] - trans[8] + vacc[1]*(1-p_u5);
      DEOD[u] = trans[8] + trans[44] - trans[9];
      DIOD[u] = k*trans[9] + trans[45] - trans[10];
      DAOD[u] = (1-k)*trans[9] + trans[46] - trans[11];
      DROD[u] = trans[10] + trans[47] - trans[12];
      DRAOD[u] = trans[11] + trans[48] - trans[13];

      DVTD[u] = trans[18] + trans[19] + trans[49] - trans[7] - trans[14] + vacc[0]*(1-p_u5);
      DETD[u] = trans[14] + trans[50] - trans[15];
      DITD[u] = k*trans[15] + trans[51] - trans[16];
      DATD[u] = (1-k)*trans[15] + trans[52] - trans[17];
      DRTD[u] = trans[16] + trans[53] - trans[18];
      DRATD[u] = trans[17] + trans[54] - trans[19];

      DVODu5[u] = trans[35] + trans[36] + trans[55] - trans[29] - trans[31] + vacc[1]*p_u5;
      DEODu5[u] = trans[31] + trans[56] - trans[32];
      DIODu5[u] = k*trans[32] + trans[57] - trans[33];
      DAODu5[u] = (1-k)*trans[32] + trans[58] - trans[34];
      DRODu5[u] = trans[33] + trans[59] - trans[35];
      DRAODu5[u] = trans[34] + trans[60] - trans[36];

      DVTDu5[u] = trans[41] + trans[42] + trans[61] - trans[30] - trans[37] + vacc[0]*p_u5;
      DETDu5[u] = trans[37] + trans[62] - trans[38];
      DITDu5[u] = k*trans[38] + trans[63] - trans[39];
      DATDu5[u] = (1-k)*trans[38] + trans[64] - trans[40];
      DRTDu5[u] = trans[39] + trans[65] - trans[41];
      DRATDu5[u] = trans[40] + trans[66] - trans[42];

      DC[u] = trans[2] + trans[10] + trans[16] + trans[33] + trans[39];

      DW[u] = trans[20] + trans[28] - trans[21];
    }
  ')

  cholera_partrans <- pomp::parameter_trans(
    log=c("sigma","gammaE","gamma","Omega1","Omega2","AlphaS",
          "Delta","ps","Sat","sigmaSE","VR","WR","Psi",
          "Mu","Beta","BetaW","v", paste0("InitInfected", 1:10)),
    logit=c("k","VE1","VE2","Rho","f", "redb", "redmu")
  )

  haiti <- haiti2_data()
  if (region == "before"){
    haiti <- haiti[haiti$year <= cutoff,]
    cholera_rinit <- cholera_rinit_before
  } else {
    haiti <- haiti[haiti$year > cutoff,]
    cholera_rinit <- cholera_rinit_after
  }

  if (measure == "log"){
    cholera_dmeasure <- cholera_log_dmeasure
  }

  # Initialize unit values based on reported cases
  c_params_IVPS <- haiti[haiti$year == min(haiti$year),]$cases / 0.2
  names(c_params_IVPS) <- paste0(c("InitInfected"),1:10)

  # MLE for epidemic model
  start_params <- c(
    sigma = 1/5, gammaE = 52.18/0.18, k = 0.2, gamma = 52.18,
    Beta = 52.18 * 9.92853814e-07, Omega1 = 1, Omega2 = 1 / 5,
    VE1 = 0.43, VE2 = 0.52, Delta = 52 / 3, Mu = 52.18 * 3.98 * 10^3,
    AlphaS = 0.4, ps = 1.00, Sat = 1e5, BetaW = 52.18 * 4.03310744e-02,
    sigmaSE = 0.01, VR = 52.18 * 10^(-12), WR = 52.18 * 0.5, Psi = 0.5,
    Rho = 0.2, scenario = 0,f = 0.75, v = 5.371e2, phase = 0,
    redb = 0.001, redmu = 0.0000001
  )
  # MLE for endemic model
  if (region=="after") {
    start_params[c("VR","Mu","v")] <- c(52.18 * 10^(-12), 52.18 * 2.58208148, 8.201e1)
  }

  par <- c(start_params,c_params_IVPS)

  ret <- spatPomp::spatPomp(
    data = haiti,
    units = "department",
    times = "year",
    t0 = min(haiti$year)-1/52,
    unit_statenames = cholera_unit_statenames,
    skeleton = pomp::vectorfield(cholera_skel),
    unit_accumvars = c("C"),
    paramnames = cholera_paramnames,
    params = par,
    partrans = cholera_partrans,
    globals = cholera_globals,
    rinit = cholera_rinit,
    dmeasure = cholera_dmeasure
  )
  ret@params <- par
  return(ret)
}
