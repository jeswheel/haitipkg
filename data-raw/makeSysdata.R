# This file creates the internal data file: R/sysdata.rda

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

MODEL3_INPUT_PARAMETERS <- yaml::read_yaml(text = model3_params_yaml)

usethis::use_data(
  MODEL3_INPUT_PARAMETERS,
  internal = TRUE,
  overwrite = TRUE
)
