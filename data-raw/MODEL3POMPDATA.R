## code to prepare `MODEL3POMPDATA` dataset goes here

library(magrittr)

# First Define some helper functions:
dateToYears <- function(date, origin = as.Date("2014-01-01"), yr_offset = 2014) {
  julian(date, origin = origin) / 365.25 + yr_offset
}

yearsToDate <- function(year_frac, origin = as.Date("2014-01-01"), yr_offset = 2014.0) {
  as.Date((year_frac - yr_offset) * 365.25, origin = origin)
}

yearsToDateTime <- function(year_frac, origin = as.Date("2014-01-01"), yr_offset = 2014.0) {
  as.POSIXct((year_frac - yr_offset) * 365.25 * 3600 * 24, origin = origin)
}


# List of all departements in Haiti
departements <- c(
  'Artibonite', 'Centre', 'Grande_Anse', 'Nippes', 'Nord',
  'Nord-Est','Nord-Ouest', 'Ouest', 'Sud', 'Sud-Est'
)

# input parameters to the model
input_parameters <- yaml::read_yaml("https://raw.githubusercontent.com/jcblemai/haiti-mass-ocv-campaign/master/haiti-data/input_parameters.yaml")

# Rain data for all of the departements:
ALL_RAIN <- readr::read_csv("https://raw.githubusercontent.com/jcblemai/haiti-mass-ocv-campaign/master/haiti-data/fromAzman/rainfall.csv")
ALL_CASES <- readr::read_csv("https://raw.githubusercontent.com/jcblemai/haiti-mass-ocv-campaign/master/haiti-data/fromAzman/cases_corrected.csv")

usethis::use_data(
  departements, input_parameters, ALL_CASES, ALL_RAIN,
  internal = TRUE, overwrite = TRUE
)



# usethis::use_data(MODEL3POMPDATA, overwrite = TRUE)
