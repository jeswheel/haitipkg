#' Function used to get department data for haiti3_dep model
#'
#' @param departement Name of department to fit POMP model to
#' @param start_time Time of which to start the time series. In the
#'    original paper, the authors used a start time of 2014-03-01 for all
#'    departements except for Ouest, which started with a time of 2017-06-10.
#'    This means that there wasn't a national-coupled model fit to the
#'    data until 2017-06-01, which means that there was less than 2 years
#'    of data to fit the model. Still, the authors fit the remaining departemental
#'    models to data from 2014-03-01. This added argument allows us to
#'    fit the Ouest model earlier, or fit the remaining departement models
#'    at 2017-06-01
#'
#' @importFrom magrittr %>%
#' @importFrom foreach %do%
#' @importFrom foreach foreach

haiti3_dep_data <- function(departement = 'Artibonite', start_time = '2014-03-01') {

  dateToYears <- function(date, origin = as.Date("2014-01-01"), yr_offset = 2014) {
    julian(date, origin = origin)/365.25 + yr_offset
  }

  yearsToDate <- function(year_frac, origin = as.Date("2014-01-01"), yr_offset = 2014.0) {
    as.Date((year_frac - yr_offset) * 365.25, origin = origin)
  }

  yearsToDateTime <- function(year_frac, origin = as.Date("2014-01-01"), yr_offset = 2014.0) {
    as.POSIXct((year_frac - yr_offset) * 365.25 * 3600 * 24, origin = origin)
  }

  # # Start and end dates of epidemic
  # if (departement == 'Ouest') {
  #   t_start <- dateToYears(as.Date('2017-06-10'))
  # } else {
  #   t_start <- dateToYears(as.Date(MODEL3_INPUT_PARAMETERS$t_start))
  # }

  t_start <- dateToYears(as.Date(start_time))
  t_end <- dateToYears(as.Date(MODEL3_INPUT_PARAMETERS$t_end))


  # haitiCholera is a saved data.frame in the package
  MODEL3_CASES <- haitiCholera %>%
    dplyr::rename(
      date = date_saturday, Grande_Anse = Grand.Anse,
      `Nord-Est` = Nord.Est, `Nord-Ouest` = Nord.Ouest,
      `Sud-Est` = Sud.Est
    ) %>%
    dplyr::mutate(date = as.Date(date)) %>%
    dplyr::select(-report)

  cases <- MODEL3_CASES %>%
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

  rain <- haitiRainfall  %>%
    tidyr::gather(dep, rain, -date) %>%
    dplyr::group_by(dep) %>%
    dplyr::ungroup() %>%
    dplyr::filter(dep == departement) %>%
    dplyr::mutate(date = as.Date(date, format = "%Y-%m-%d"),
           time = dateToYears(date)) %>%
    dplyr::filter(time > t_start - 0.01 & time < (t_end + 0.01)) %>%
    dplyr::mutate(max_rain = max(rain), rain_std = rain/max_rain)


  cases_other_dept <- MODEL3_CASES %>%
    tidyr::gather(dep, cases, -date) %>%
    dplyr::filter(dep != departement) %>%
    dplyr::mutate(date = as.Date(date, format = "%Y-%m-%d"),
           time = dateToYears(date))

  cases_other_dept <- aggregate(cases_other_dept$cases, by=list(Category=cases_other_dept$time), FUN=sum, na.rm=TRUE, na.action=NULL) %>%
    dplyr::mutate(time = Category) %>%
    dplyr::mutate(cases = x)

  cases_covar <- cases_other_dept %>%
    dplyr::select(time, cases) %>%
    dplyr::rename(cases_other = cases)

  ret <- list()
  ret$cases_covar <- cases_covar
  ret$rain <- rain
  ret$cases <- cases
  # ret$input_parameters <- MODEL3_INPUT_PARAMETERS
  # ret$cases_other_dept <- cases_other_dept
  ret$params <- MODEL3_PARAMS

  ret
}
