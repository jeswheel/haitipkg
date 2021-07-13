#' Function used to get department data for haiti3_dep model
#'
#' @importFrom magrittr %>%
#' @importFrom foreach %do%
#' @importFrom foreach foreach

haiti3_dep_data <- function(departement = 'Artibonite') {

  dateToYears <- function(date, origin = as.Date("2014-01-01"), yr_offset = 2014) {
    julian(date, origin = origin)/365.25 + yr_offset
  }

  yearsToDate <- function(year_frac, origin = as.Date("2014-01-01"), yr_offset = 2014.0) {
    as.Date((year_frac - yr_offset) * 365.25, origin = origin)
  }

  yearsToDateTime <- function(year_frac, origin = as.Date("2014-01-01"), yr_offset = 2014.0) {
    as.POSIXct((year_frac - yr_offset) * 365.25 * 3600 * 24, origin = origin)
  }

  input_parameters <- yaml::read_yaml("https://raw.githubusercontent.com/ionides/haiti/main/model3/input/haiti-data/input_parameters.yaml?token=AQ356PMPCSG54WKPJDZZ6I3A64K2C")

  # Start and end dates of epidemic
  if (departement == 'Ouest') {
    t_start <- dateToYears(as.Date('2017-06-10'))
  } else {
    t_start <- dateToYears(as.Date(input_parameters$t_start))
  }
  t_end <- dateToYears(as.Date(input_parameters$t_end))


  CASES_TEMP <- readr::read_csv("https://raw.githubusercontent.com/ionides/haiti/main/model3/input/haiti-data/fromAzman/cases_corrected.csv?token=AQ356PNOJAOF4MGB2BIWCK3A64K7A")

  cases <- CASES_TEMP %>%
    tidyr::gather(dep, cases, -date) %>%
    dplyr::filter(dep == departement) %>%
    dplyr::mutate(date = as.Date(date, format = "%Y-%m-%d"),
                  time = dateToYears(date))

  case_dates <- with(
    cases %>%
      dplyr::filter(time > t_start - 0.01 & time < (t_end + 0.01)),
    seq.Date(min(date), max(date), by = "1 week")
  )

  missing_dates <- setdiff(case_dates, cases$date) %>% as.Date(origin = as.Date("1970-01-01"))

  rain <- readr::read_csv("https://raw.githubusercontent.com/ionides/haiti/main/model3/input/haiti-data/fromAzman/rainfall.csv?token=AQ356PM3DCWJXDUWFXGA5ETA64LLM")  %>%
    tidyr::gather(dep, rain, -date) %>%
    dplyr::group_by(dep) %>%
    dplyr::ungroup() %>%
    dplyr::filter(dep == departement) %>%
    dplyr::mutate(date = as.Date(date, format = "%Y-%m-%d"),
           time = dateToYears(date)) %>%
    dplyr::filter(time > t_start - 0.01 & time < (t_end + 0.01)) %>%
    dplyr::mutate(max_rain = max(rain), rain_std = rain/max_rain)


  cases_other_dept <- CASES_TEMP %>%
    tidyr::gather(dep, cases, -date) %>%
    dplyr::filter(dep != departement) %>%
    dplyr::mutate(date = as.Date(date, format = "%Y-%m-%d"),
           time = dateToYears(date))

  cases_other_dept <- aggregate(cases_other_dept$cases, by=list(Category=cases_other_dept$time), FUN=sum, na.rm=TRUE, na.action=NULL) %>%
    dplyr::mutate(time = Category) %>%
    dplyr::mutate(cases = x)

  cases_covar <- cases_other_dept

  params <- read.csv('https://raw.githubusercontent.com/ionides/haiti/main/model3/output/best_artibonite_params.csv?token=AQ356PMHAND3J5L5OSGYLR3A64QJA')

  ret <- list()
  ret$cases_covar <- cases_covar
  ret$rain <- rain
  ret$cases <- cases
  ret$input_parameters <- input_parameters
  ret$cases_other_dept <- cases_other_dept
  ret$params <- params

  ret
}
