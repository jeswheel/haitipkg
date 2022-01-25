# This file creates the internal data file: R/sysdata.rda

#### MODEL 1 INTERNAL DATA #####################################################
model1_params_yaml <- "
# supplement model parameters from original model 1 team
supp_pars_epi:
  rho: 0.34
  tau: 4.0
  beta1: 5.8
  beta2: 3.5
  beta3: 4.4
  beta4: 2.8
  beta5: 6.2
  beta6: 2.4
  gamma: 3.5
  sigma: 5.0
  theta0: 0.0
  alpha: 0.0023973
  mu: 0.0004287
  delta: 0.0001433
  nu: 0.96
  S_0: 0.9999154
  E_0: 3.7115718e-05
  I_0: 4.7471462e-05
  A_0: 0.0e+00
  R_0: 0.0e+00
  pop_0: 1.0911819e+07

supp_pars_end:
  rho: 0.34
  tau: 4.1
  beta1: 5.4
  beta2: 3.3
  beta3: 4.3
  beta4: 3.5
  beta5: 5.0
  beta6: 2.7
  gamma: 3.5
  sigma: 5.0
  theta0: 0.0
  alpha: 0.0023973
  mu: 0.0004287
  delta: 0.0001433
  nu: 0.96
  S_0: 0.9641389
  E_0: 2.1271687e-05
  I_0: 2.7790429e-05
  A_0: 0.0e+00
  R_0: 1.2158502e-01
  pop_0: 1.165869e+07

# adjusted model parameters
adj_pars_epi:
  rho: 0.3148019
  tau: 376.7801647
  beta1: 5.3323503
  beta2: 2.6565848
  beta3: 3.8324673
  beta4: 2.7666027
  beta5: 5.0973542
  beta6: 1.8033846
  nu: 0.9840841
  gamma: 3.5
  sigma: 5.0
  theta0: 0.0
  alpha: 0.0023973
  mu: 0.0004287
  delta: 0.0001433
  sig_sq: 0.1016082
  S_0: 0.9985646
  E_0: 2.1477506e-06
  I_0: 1.43324e-03
  A_0: 0.0e+00
  R_0: 0.0e+00
  pop_0: 1.0911819e+07

adj_pars_end:
  rho: 9.517093e-01
  tau: 8.5445968e+01
  beta1: 2.4295772e+00
  beta2: 4.1215816e+00
  beta3: 2.0810796e+00
  beta4: 3.7738013e+00
  beta5: 2.4401692e+00
  beta6: 3.603727e+00
  nu: 9.8693735e-01
  gamma: 3.5e+00
  sigma: 5.0e+00
  theta0: 0.0e+00
  alpha: 2.3972603e-03
  mu: 4.2871488e-04
  delta: 1.4331704e-04
  sig_sq: 1.1219656e-01
  S_0: 9.5181736e-01
  E_0: 1.4667053e-05
  I_0: 2.0413793e-05
  A_0: 0.0e+00
  R_0: 4.814756e-02
  pop_0: 1.1658784e+07


## joint model parameters
joint_pars:
  rho_epi: 0.4765437
  rho_end: 0.4496893
  tau_epi: 688.7796
  tau_end: 105.3583
  sig_sq_epi: 0.1105648
  sig_sq_end: 0.1677307
  beta1: 4.014758
  beta2: 2.7089
  beta3: 2.742331
  beta4: 3.058927
  beta5: 3.57466
  beta6: 2.230872
  nu: 0.9976078
  gamma: 3.5
  sigma: 5.0
  theta0: 0.0
  alpha: 0.0023973
  mu: 0.0004287
  delta: 0.0001433
  kappa: 0.0
  S_0: 0.9990317
  E_0: 4.604823e-06
  I_0: 9.63733e-04
  A_0: 0.0e+00
  R_0: 0.0e+00
  pop_0: 1.0911819e+07

dep_params:
  Artibonite:
    rho_epi: 0.383749
    tau_epi: 145.1022583
    beta1: 7.5787035
    beta2: 4.3469957
    beta3: 6.5306949
    beta4: 3.3930507
    beta5: 7.9882798
    beta6: 3.9106257
    nu: 0.9304411
    gamma: 3.5
    sigma: 5.0
    theta0: 0.0
    alpha: 0.0023973
    mu: 0.0004287
    delta: 0.0001433
    sig_sq_epi: 0.1185599
    kappa: 0.0
    S_0: 0.9874133
    E_0: 4.0808052e-05
    I_0: 1.2545842e-02
    A_0: 0.0e+00
    R_0: 0.0e+00
    pop_0: 7.46236e+05
    rho_end: 1.3076759e-01
    tau_end: 1.8921246e+01
    sig_sq_end: 2.1124445e-01
    mob_c: 0.0e+00
  Centre:
    rho_epi: 0.3340463
    tau_epi: 78.6370229
    beta1: 3.6488766
    beta2: 3.2993265
    beta3: 3.1930774
    beta4: 3.5528296
    beta5: 3.7327194
    beta6: 2.9777486
    nu: 0.9722747
    gamma: 3.5
    sigma: 5.0
    theta0: 0.0
    alpha: 0.0023973
    mu: 0.0004287
    delta: 0.0001433
    sig_sq_epi: 0.1308516
    kappa: 0.0
    S_0: 0.9991466
    E_0: 1.5191083e-05
    I_0: 8.3816666e-04
    A_0: 0.0e+00
    R_0: 0.0e+00
    pop_0: 1.727524e+06
    rho_end: 2.2137364e-01
    tau_end: 2.1720161e+01
    sig_sq_end: 1.2110847e-01
    mob_c: 0.0e+00
  Grand_Anse:
    rho_epi: 0.0096261
    tau_epi: 8971.9304261
    beta1: 10.9681899
    beta2: 5.8077791
    beta3: 7.0497384
    beta4: 4.6017879
    beta5: 9.4181223
    beta6: 2.4087516
    nu: 0.9730264
    gamma: 3.5
    sigma: 5.0
    theta0: 0.0
    alpha: 0.0023973
    mu: 0.0004287
    delta: 0.0001433
    sig_sq_epi: 0.2526313
    kappa: 0.0
    S_0: 0.9999989
    E_0: 5.957156e-07
    I_0: 5.2597349e-07
    A_0: 0.0e+00
    R_0: 0.0e+00
    pop_0: 4.029705e+06
    rho_end: 7.259955e-03
    tau_end: 9.0939744e+00
    sig_sq_end: 2.9558153e-01
    mob_c: 0.0e+00

  Nippes:
    rho_epi: 0.0220382
    tau_epi: 74.1086007
    beta1: 12.2065332
    beta2: 5.1698711
    beta3: 8.2117033
    beta4: 2.5695345
    beta5: 12.1397699
    beta6: 1.1453644
    nu: 0.9203798
    gamma: 3.5
    sigma: 5.0
    theta0: 0.0
    alpha: 0.0023973
    mu: 0.0004287
    delta: 0.0001433
    sig_sq_epi: 0.284506
    kappa: 0.0
    S_0: 0.9999948
    E_0: 3.3703587e-06
    I_0: 1.8125202e-06
    A_0: 0.0e+00
    R_0: 0.0e+00
    pop_0: 7.28807e+05
    rho_end: 6.5109067e-02
    tau_end: 4.8582933e+00
    sig_sq_end: 6.8280776e-01
    mob_c: 0.0e+00
  Nord:
    rho_epi: 0.15454
    tau_epi: 3012.64988
    beta1: 4.7313415
    beta2: 6.6430733
    beta3: 2.8213788
    beta4: 5.8651783
    beta5: 4.4868882
    beta6: 4.0658614
    nu: 0.9587818
    gamma: 3.5
    sigma: 5.0
    theta0: 0.0
    alpha: 0.0023973
    mu: 0.0004287
    delta: 0.0001433
    sig_sq_epi: 0.153497
    kappa: 0.0
    S_0: 0.999932
    E_0: 3.8942537e-06
    I_0: 6.4063859e-05
    A_0: 0.0e+00
    R_0: 0.0e+00
    pop_0: 1.067177e+06
    rho_end: 6.2013574e-02
    tau_end: 4.7345686e+02
    sig_sq_end: 2.2718365e-01
    mob_c: 0.0e+00
  Nord_Est:
    rho_epi: 0.0999359
    tau_epi: 645.7778126
    beta1: 6.2914301
    beta2: 5.446388
    beta3: 4.4770558
    beta4: 5.9093195
    beta5: 4.3522892
    beta6: 5.3328708
    nu: 0.9462916
    gamma: 3.5
    sigma: 5.0
    theta0: 0.0
    alpha: 0.0023973
    mu: 0.0004287
    delta: 0.0001433
    sig_sq_epi: 0.2362206
    kappa: 0.0
    S_0: 0.9999953
    E_0: 1.8550553e-06
    I_0: 2.812548e-06
    A_0: 0.0e+00
    R_0: 0.0e+00
    pop_0: 7.74976e+05
    rho_end: 1.0607023e-02
    tau_end: 9.3793768e+00
    sig_sq_end: 3.9945573e-01
    mob_c: 0.0e+00
  Nord_Ouest:
    rho_epi: 0.1325456
    tau_epi: 23.1387759
    beta1: 1.6937492
    beta2: 14.1274657
    beta3: 1.1028379
    beta4: 12.2788406
    beta5: 0.7625294
    beta6: 12.2616334
    nu: 0.9523155
    gamma: 3.5
    sigma: 5.0
    theta0: 0.0
    alpha: 0.0023973
    mu: 0.0004287
    delta: 0.0001433
    sig_sq_epi: 0.2488908
    kappa: 0.0
    S_0: 0.9998923
    E_0: 2.8839475e-05
    I_0: 7.882876e-05
    A_0: 0.0e+00
    R_0: 0.0e+00
    pop_0: 3.42525e+05
    rho_end: 4.5470814e-02
    tau_end: 1.6969595e+01
    sig_sq_end: 3.1347886e-01
    mob_c: 0.0e+00
  Ouest:
    rho_epi: 0.9826225
    tau_epi: 90.5567277
    beta1: 9.6423596
    beta2: 4.7487868
    beta3: 6.9562077
    beta4: 6.0106314
    beta5: 8.5357427
    beta6: 3.5816965
    nu: 0.9467228
    gamma: 3.5
    sigma: 5.0
    theta0: 0.0
    alpha: 0.0023973
    mu: 0.0004287
    delta: 0.0001433
    sig_sq_epi: 0.1935664
    kappa: 0.0
    S_0: 0.99944
    E_0: 0.0002373
    I_0: 0.0003227
    A_0: 0.0
    R_0: 0.0
    pop_0: 393967.0
    rho_end: 0.3309908
    tau_end: 11.2499227
    sig_sq_end: 0.2430906
    mob_c: 0.0
  Sud:
    rho_epi: 0.0949872
    tau_epi: 338.1798452
    beta1: 4.5883302
    beta2: 10.9440715
    beta3: 0.2777675
    beta4: 9.5491781
    beta5: 3.1053673
    beta6: 8.4701333
    nu: 0.9384123
    gamma: 3.5
    sigma: 5.0
    theta0: 0.0
    alpha: 0.0023973
    mu: 0.0004287
    delta: 0.0001433
    sig_sq_epi: 0.2374588
    kappa: 0.0
    S_0: 0.9999979
    E_0: 1.6680161e-06
    I_0: 4.3539489e-07
    A_0: 0.0e+00
    R_0: 0.0e+00
    pop_0: 6.32601e+05
    rho_end: 4.1989146e-02
    tau_end: 4.1129264e+00
    sig_sq_end: 3.4645195e-01
    mob_c: 0.0e+00
  Sud_Est:
    rho_epi: 0.0635203
    tau_epi: 10.0311336
    beta1: 4.59917
    beta2: 8.0904975
    beta3: 1.5316885
    beta4: 7.5833201
    beta5: 4.2558375
    beta6: 5.4098747
    nu: 0.9400739
    gamma: 3.5
    sigma: 5.0
    theta0: 0.0
    alpha: 0.0023973
    mu: 0.0004287
    delta: 0.0001433
    sig_sq_epi: 0.2554795
    kappa: 0.0
    S_0: 0.9999883
    E_0: 8.2360658e-06
    I_0: 3.5107046e-06
    A_0: 0.0e+00
    R_0: 0.0e+00
    pop_0: 4.68301e+05
    rho_end: 3.5157508e-02
    tau_end: 4.7404628e+00
    sig_sq_end: 3.0963919e-01
    mob_c: 0.0e+00

## panel model parameters
panel_shared_params:
  mu: 0.0004287
  delta: 0.0001433
  gamma: 3.5
  sigma: 5.0
  theta0: 0.0
  alpha: 0.0023973
  kappa: 0
  beta1: 4.014758
  beta2: 2.7089
  beta3: 2.742331
  beta4: 3.058927
  beta5: 3.57466
  beta6: 2.230872
  nu: 0.9976078
  tau_epi: 688.7796
  tau_end: 105.3583
  sig_sq_epi: 0.1105648
  sig_sq_end: 0.1677307
"


#### MODEL 3 INTERNAL DATA #####################################################
model3_params_yaml <- "
# INPUTA PARAMTERS FOR CHOLERA MODEL

# Note: all parameters are given in per day (/day) units except if mentioned otherwise

# total population

population:
  - Artibonite: 1727524
  - Centre: 746236
  - Grande_Anse: 468301
  - Nippes: 342525
  - Nord: 1067177
  - Nord-Est: 393967
  - Nord-Ouest: 728807
  - Ouest: 4029705
  - Sud: 774976
  - Sud-Est: 632601

density:
  - Artibonite: 353
  - Centre: 214
  - Grande_Anse: 245
  - Nippes: 270
  - Nord: 505
  - Nord-Est: 243
  - Nord-Ouest: 347
  - Ouest: 809
  - Sud: 292
  - Sud-Est: 311

# Summary previous vaccination campagain:
# Grande Anse:
#Target pop :  37616 + 32635 + 8900 + 27342 + 40041 + 138802 + 24155 + 32197 + 46149 = 387837
# p2d =       (21255 + 17270 + 5688 + 18979 + 26480 + 90245 + 17882 + 19732 + 26897) /387837 = 0.63
# Sud:
# Target Pop: 107695 + 46547 + 26083 + 22014 + 30888 + 156762 + 19554 + 19736 + 20437 = 449716
#p2d = (70633 +  39029 +  15989 +16487 + 19300 + 110318 + 11748 + 16017 + 11366) / 449716 = 0.69
# Ouest:
# target pop 4053
#p2d = 2829/ 4053 = 0.70
# Centre
# Target pop 98563
# p2d  69905/98563 = 0.71
# Artibonite
# Target pop 62241
# p2d  59537/62241 = 0.95.

t_vacc_start_alt:
  - Artibonite: 2018-04-17
  - Centre: 2017-11-15
  - Grande_Anse: 2016-11-08
  - Nippes: 2010-01-01 # No previous campagin
  - Nord: 2010-01-01 # No previous campagin
  - Nord-Est: 2010-01-01 # No previous campagin
  - Nord-Ouest: 2010-01-01 # No previous campagin
  - Ouest: 2017-07-25
  - Sud: 2016-11-08
  - Sud-Est: 2010-01-01 # No previous campagin

t_vacc_end_alt:
  - Artibonite: 2018-05-15
  - Centre: 2017-12-16
  - Grande_Anse: 2017-06-02
  - Nippes: 2000-01-01 # No previous campagin
  - Nord: 2000-01-01 # No previous campagin
  - Nord-Est: 2000-01-01 # No previous campagin
  - Nord-Ouest: 2000-01-01 # No previous campagin
  - Ouest: 2017-08-25
  - Sud: 2017-05-31
  - Sud-Est: 2000-01-01 # No previous campagin

p1d_alt_year:
  - Artibonite: 0.05
  - Centre: 0.29
  - Grande_Anse: 0.37
  - Nippes: 0 # No previous campagin
  - Nord: 0 # No previous campagin
  - Nord-Est: 0 # No previous campagin
  - Nord-Ouest: 0 # No previous campagin
  - Ouest: 0.3
  - Sud: 0.31
  - Sud-Est: 0 # No previous campagin

nb_doses_alt_year:
  - Artibonite: 62241
  - Centre: 98563
  - Grande_Anse: 387837
  - Nippes: 0 # No previous campagin
  - Nord: 0 # No previous campagin
  - Nord-Est: 0 # No previous campagin
  - Nord-Ouest: 0 # No previous campagin
  - Ouest: 4053
  - Sud: 449716
  - Sud-Est: 0 # No previous campagin

# rate of loss of infection
#gammaA: 0.2
#gammaI: 0.2
# natural mortality rate in year-1
mu: 0.01586625546    # population natality and mortality rate (life exp is 63.07 year) (day^-1) 1/(63.07*365) *365.25
# rate of cholera-induced mortality in year^-1
alpha: 1.461            # mortality rate due to cholera (day^-1) (2% case fatality rate) * 365.25

# start and end times of epidemic
t_start: 2014-03-01
#t_start: 2017-06-10   # For West
t_end: 2019-01-12

#t_end: 2017-01-01 # that is before peak
"

model3_vacc_scenarios_txt <- "
ID,Roll-out,VE,Coverage,Priority
1,1,1,1,1
2,2,1,1,1
3,3,1,1,1
4,4,1,1,1
5,1,2,1,0
6,2,2,1,0
7,3,2,1,0
8,4,2,1,0
9,1,3,1,0
10,2,3,1,0
11,3,3,1,0
12,4,3,1,0
13,1,1,2,0
14,2,1,2,0
15,3,1,2,0
16,4,1,2,0
17,1,2,2,0
18,2,2,2,0
19,3,2,2,0
20,4,2,2,0
21,1,3,2,0
22,2,3,2,0
23,3,3,2,0
24,4,3,2,0
25,1,1,3,1
26,2,1,3,0
27,3,1,3,0
28,4,1,3,0
29,1,2,3,0
30,2,2,3,0
31,3,2,3,0
32,4,2,3,0
33,1,3,3,0
34,2,3,3,0
35,3,3,3,0
36,4,3,3,0
"

model3_params_text <- '
betaB,mu_B,XthetaA,thetaI,lambdaR,std_W,epsilon,k,cas_def,foi_add,sigma,r,gammaI,gammaA,rhoA,XrhoI,Rtot_0,H,D,mu,alpha,cases_ext,t_vacc_start,t_vacc_end,p1d_reg,r_v_year
0.354169973579405,111.359282045485,0.0577953719949065,0.000432765266676422,0.499169039613328,0.01400411832052,0.993716987777945,516.87552896729,0.934120456253182,1.18214116603088e-06,0.25,0.6359411,182.625,182.625,0.125085616438356,1,0.35,1727524,353,0.01586625546,1.461,1,0,0,0,0
'


MODEL1_INPUT_PARAMETERS <- yaml::read_yaml(text = model1_params_yaml)
MODEL3_INPUT_PARAMETERS <- yaml::read_yaml(text = model3_params_yaml)
MODEL3_VACC_SCENARIOS   <- read.csv(text = model3_vacc_scenarios_txt)
MODEL3_PARAMS           <- read.csv(text = model3_params_text)


usethis::use_data(
  MODEL1_INPUT_PARAMETERS,
  MODEL3_INPUT_PARAMETERS,
  MODEL3_VACC_SCENARIOS,
  MODEL3_PARAMS,
  internal = TRUE,
  overwrite = TRUE
)
