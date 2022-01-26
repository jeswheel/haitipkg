#' Get vaccination scenario parameters
#'
#' This function returns the parameters that are needed to run the various
#' vaccination scenarios for model 3.
#'
#' @param scenario_str character vector containing the name of each vaccination
#' scenario that will be simulated. The options are:
#' TODO: update this list
#' * S0: Default vaccination scenario
#' * S1
#' * S2
#' * S3
#' * S4
#' * S25
#'
#' @importFrom magrittr %>%
#' @importFrom foreach %do%
#' @importFrom foreach foreach
#' @export

get_vacc_scenario_params3 <- function(scenario_str = "S1") {

  # Create 6 different vaccination scenarios:
  #   S0, S1, S2, S3, S4, S25.
  # These scenarios are "encoded" in the file MODEL3_VACC_SCENARIOS,
  # and the lines below simply "decode" the scenarios.


  # Auxillary Helper functions
  dateToYears <- function(date, origin = as.Date("2014-01-01"), yr_offset = 2014) {
    # Converts Date object to decimal years.
    julian(date, origin = origin) / 365.25 + yr_offset
  }

  # Save convenient list of all of the department names
  DEPARTMENTS <-
    c(
      'Artibonite',
      'Centre',
      'Grande_Anse',
      'Nippes',
      'Nord',
      'Nord-Est',
      'Nord-Ouest',
      'Ouest',
      'Sud',
      'Sud-Est'
    )

  # Create list to store vaccination scenarios
  scenarios <- list()

  # Loop through all encoded vaccination scenarios and only keep the ones that
  # we are actually interested. They are marked with a priority == 1.
  for (i in 1:nrow(MODEL3_VACC_SCENARIOS)) {
    not_dep <- c()  # By default, include all departements.
    course_year <- 2  # By default, the campaign is 2 years.

    if (MODEL3_VACC_SCENARIOS[i, 'Roll.out'] == 2) {
      not_dep <- c('Ouest', 'Nord-Ouest', 'Nord', 'Sud', 'Nippes', 'Nord-Est', 'Sud-Est', 'Grande_Anse')
    } else if (MODEL3_VACC_SCENARIOS[i, 'Roll.out'] == 3) {
      course_year <- 5
    } else if (MODEL3_VACC_SCENARIOS[i, 'Roll.out'] == 4) {
      not_dep <- c('Nord-Ouest', 'Nord', 'Sud', 'Nippes', 'Nord-Est', 'Sud-Est', 'Grande_Anse')
    }

    # Once again setting defaults.
    percent_completely_unvaccinated <- 0
    percent_onedose <- 0
    percent_twodoses <- 0

    if (MODEL3_VACC_SCENARIOS[i, 'Coverage'] == 1) {
      percent_completely_unvaccinated <- 20
      percent_onedose <- 10
      percent_twodoses <- 70
    } else if (MODEL3_VACC_SCENARIOS[i, 'Coverage'] == 2) {
      percent_completely_unvaccinated <- 40
      percent_onedose <- 20
      percent_twodoses <- 40
    } else if (MODEL3_VACC_SCENARIOS[i, 'Coverage'] == 3) {
      percent_completely_unvaccinated <- 3.33
      percent_onedose <- 1.67
      percent_twodoses <- 95
    }

    # 1 for all campaigns we are interested in.
    ve <- MODEL3_VACC_SCENARIOS[i, 'VE']

    # Priority == 1 corresponds to the S1, S2, S3, S4, and S25 campaigns.
    if (MODEL3_VACC_SCENARIOS[i, 'Priority'] == 1) {
      sid <- paste0('S', MODEL3_VACC_SCENARIOS[i, 'ID'])

      # Only save and run input scenarios.
      if (sid == scenario_str) {
        scenarios[[sid]] <- VaccinationScenario3(
          course_year = course_year,
          percent_uv = percent_completely_unvaccinated,
          percent_1d = percent_onedose, percent_2d = percent_twodoses,
          ve = ve, not_dep = not_dep
        )
        break
      }
    }
  }


  if (scenario_str == "S0") {
    # Create S0 campagin and name it.
    S0 <- VaccinationScenario3(50,
                               99.9999999,
                               0.00000001,
                               0.00,
                               ve = 1)

    scenarios[['S0']] <- S0
  }

  run_scenario <- scenarios[[scenario_str]]

  t_vacc_start <- run_scenario$t_vacc_start
  t_vacc_end <- run_scenario$t_vacc_end
  p1d_reg <- run_scenario$p1d_reg
  r_v_year <- run_scenario$r_v_year
  cases_ext <- run_scenario$ve

  params_out <- c()
  for (dp in DEPARTMENTS) {
    tv_start <- dateToYears(as.Date(t_vacc_start[gsub('-','_', dp)][[1]][1]))
    names(tv_start) <- paste0("t_vacc_start", gsub('-','_', dp))

    tv_end <- dateToYears(as.Date(t_vacc_end[gsub('-','_', dp)][[1]][1]))
    names(tv_end) <- paste0("t_vacc_end", gsub('-','_', dp))

    p1d <- as.numeric(p1d_reg[gsub('-','_', dp)][[1]][1])
    names(p1d) <- paste0("p1d_reg", gsub('-','_', dp))

    rv <- as.numeric(r_v_year[gsub('-','_', dp)][[1]][1])
    names(rv) <- paste0("r_v_year", gsub('-','_', dp))

    params_out <- c(params_out, tv_start, tv_end, p1d, rv)
  }

  params_out
}
