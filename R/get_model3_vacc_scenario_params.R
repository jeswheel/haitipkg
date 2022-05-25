#' Model3 Vaccination Scenario Parameters
#'
#' This function returns the parameters that are necessary to run the vaccination
#' scenarios (in spatPomp form).
#'
#' @param scenario Vaccination scenario that is to be ran.
#'
#' @export

get_model3_vacc_scenario_params <- function(scenario = c('noVacc', '2dep', '3dep',
                                                         'slowNation', 'fastNation',
                                                         'highFastNation')) {

  # Auxillary Helper functions
  dateToYears <- function(date, origin = as.Date("2014-01-01"), yr_offset = 2014) {
    # Converts Date object to decimal years.
    julian(date, origin = origin) / 365.25 + yr_offset
  }

  yearsToDate <- function(year_frac, origin = as.Date("2014-01-01"), yr_offset = 2014.0) {
    # Converts decimal years to Date object
    as.Date((year_frac - yr_offset) * 365.25, origin = origin)
  }

  yearsToDateTime <- function(year_frac, origin = as.Date("2014-01-01"), yr_offset = 2014.0) {
    # Converts decimal years to DateTime object.
    as.POSIXct((year_frac - yr_offset) * 365.25 * 3600 * 24, origin = origin)
  }


  if (scenario == '2dep') {
    params_list <- VaccinationScenario3(
      course_year = 2,
      percent_uv = 20,
      percent_1d = 10,
      percent_2d = 70,
      ve = 1,
      not_dep = c('Ouest', 'Nord_Ouest', 'Nord', 'Sud', 'Nippes', 'Nord_Est', 'Sud_Est', 'Grande_Anse')
    )
  } else if (scenario == '3dep') {
    params_list <- VaccinationScenario3(
      course_year = 2,
      percent_uv = 20,
      percent_1d = 10,
      percent_2d = 70,
      ve = 1,
      not_dep = c('Nord_Ouest', 'Nord', 'Sud', 'Nippes', 'Nord_Est', 'Sud_Est', 'Grande_Anse')
    )
  } else if (scenario == 'slowNation') {
    params_list <- VaccinationScenario3(
      course_year = 5,
      percent_uv = 20,
      percent_1d = 10,
      percent_2d = 70,
      ve = 1,
      not_dep = c()
    )
  } else if (scenario == 'fastNation') {
    params_list <- VaccinationScenario3(
      course_year = 2,
      percent_uv = 20,
      percent_1d = 10,
      percent_2d = 70,
      ve = 1,
      not_dep = c()
    )
  } else if (scenario == 'highFastNation') {
    params_list <- VaccinationScenario3(
      course_year = 2,
      percent_uv = 3.33,
      percent_1d = 1.67,
      percent_2d = 95,
      ve = 1,
      not_dep = c()
    )
  } else if (scenario == 'noVacc') {
    params_list <- VaccinationScenario3(
      50,
      99.9999999,
      0.00000001,
      0.00,
      ve = 1
    )
  } else {
    stop("Not a valid Vaccination Scenario")
  }

  departements <-
    c(
      'Artibonite',
      'Centre',
      'Grande_Anse',
      'Nippes',
      'Nord',
      'Nord_Est',
      'Nord_Ouest',
      'Ouest',
      'Sud',
      'Sud_Est'
    )



  t_vacc_start <- params_list$t_vacc_start
  t_vacc_end <- params_list$t_vacc_end
  p1d_reg <- params_list$p1d_reg
  r_v_year <- params_list$r_v_year
  # cases_ext <- params_list$ve

  params <- rep(0, 10 * 4)
  names(params) <- c(
    paste0('t_vacc_start', 1:10),
    paste0('t_vacc_end', 1:10),
    paste0('p1d_reg', 1:10),
    paste0('r_v_year', 1:10)
  )

  for (i in 1:length(departements)) {
    dp <- departements[i]
    params[paste0('t_vacc_start', i)] <- as.numeric(dateToYears(as.Date(t_vacc_start[[dp]])))
    params[paste0('t_vacc_end', i)] <- as.numeric(dateToYears(as.Date(t_vacc_end[[dp]])))
    params[paste0('p1d_reg', i)] <- p1d_reg[[dp]]
    params[paste0('r_v_year', i)] <- r_v_year[[dp]]
  }

  params
}
