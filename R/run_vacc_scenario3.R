#' Simulate Vaccination Campaign.
#'
#' @param scenario_strs character vector containing the name of each vaccination
#' scenario that will be simulated. The options are:
#' TODO: update this list
#' * S0: Default vaccination scenario
#' * S1
#' * S2
#' * S3
#' * S4
#' * S25
#' @param nsim integer number of simulations for each vaccination scenario.
#'
#' @importFrom magrittr %>%
#' @importFrom foreach %do%
#' @importFrom foreach foreach
#' @export

run_vacc_scenario3 <- function(scenario_strs = c("S0", "S1", "S2", "S3",
                                                 "S4", "S25"),
                               nsim = 20) {

  # Create 6 different vaccination scenarios:
  #   S0, S1, S2, S3, S4, S25.
  # These scenarios are "encoded" in the file MODEL3_VACC_SCENARIOS,
  # and the lines below simply "decode" the scenarios.


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

  # Load the input parameters
  input_parameters <- MODEL3_INPUT_PARAMETERS

  # Start and end dates of epidemic
  t_start <- dateToYears(as.Date(input_parameters$t_start))
  t_end <- dateToYears(as.Date(input_parameters$t_end))
  t_forecast <- dateToYears(as.Date("2029-12-21"))  # TODO: Make as a function parameter.

  # Get the times that we would like to forecast
  time_forecast <- dateToYears(seq.Date(yearsToDate(t_start), yearsToDate(t_forecast), by = "1 week"))

  std_rain <- function(x) {
    # This function simply standardizes the rain for us.
    x / max(x)
  }

  all_rain <- haitiRainfall %>%
    summarize(
      date = date, across(Artibonite:`Sud-Est`, std_rain)
    ) %>%
    mutate(
      time = dateToYears(date)
    ) %>%
    dplyr::filter(time > t_start - 0.01 & time < (t_end + 0.01))

  colnames(all_rain) <- c(
    "date",
    paste0(
      'rain_std', c(
        'Artibonite', 'Centre', 'Grande_Anse',
        'Nippes', 'Nord', 'Nord_Est', 'Nord_Ouest',
        'Ouest', 'Sud', 'Sud_Est'
      )
    ),
    'time'
  )

  all_rain_forecast <- project_rain() %>%
    dplyr::filter(time > t_start - 0.01 & time < (t_forecast + 0.01))

  # Create list to store vaccination scenarios
  scenarios <- list()

  # Loop through all encoded vaccination scenarios and only keep the ones that
  # we are actually interested. They are marked with a priority == 1.
  for (i in 1:nrow(MODEL3_VACC_SCENARIOS)) {
    not_dep <- c()  # By default, include all departements.
    course_year <- 2  # By default, the campaign is 2 years.

    if (MODEL3_VACC_SCENARIOS[i, 'Roll-out'] == 2) {
      not_dep <- c('Ouest', 'Nord-Ouest', 'Nord', 'Sud', 'Nippes', 'Nord-Est', 'Sud-Est', 'Grande_Anse')
    } else if (MODEL3_VACC_SCENARIOS[i, 'Roll-out'] == 3) {
      course_year <- 5
    } else if (MODEL3_VACC_SCENARIOS[i, 'Roll-out'] == 4) {
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
      if (sid %in% scenario_strs) {
        scenarios[[sid]] <- VaccinationScenario3(
          course_year = course_year,
          percent_uv = percent_completely_unvaccinated,
          percent_1d = percent_onedose, percent_2d = percent_twodoses,
          ve = ve, not_dep = not_dep
        )
      }
    }
  }


  if (any(scenario_strs == "S0")) {
    # Create S0 campagin and name it.
    S0 <- VaccinationScenario3(50,
                               99.9999999,
                               0.00000001,
                               0.00,
                               ve = 1)

    scenarios[['S0']] <- S0
  }

  mod <- haiti3()

  first_flag <- TRUE

  # Loop through each of the scenarios and simulate the model.
  all_sims <- foreach(scenario_str = scenario_strs,
                      .combine = dplyr::bind_rows) %do% {

    if (first_flag) {
      first_flag <- FALSE
      include_data <- TRUE
    } else {
      include_data <- FALSE
    }

    cat('Working on Scenario', scenario_str, '\n')
    run_scenario <- scenarios[[scenario_str]]

    t_vacc_start <- run_scenario$t_vacc_start
    t_vacc_end <- run_scenario$t_vacc_end
    p1d_reg <- run_scenario$p1d_reg
    r_v_year <- run_scenario$r_v_year
    cases_ext <- run_scenario$ve

    for (dp in DEPARTMENTS) {
      pomp::coef(mod)[paste0("t_vacc_start", gsub('-','_', dp))] <- dateToYears(as.Date(t_vacc_start[gsub('-','_', dp)][[1]][1]))
      pomp::coef(mod)[paste0("t_vacc_end", gsub('-','_', dp))] <- dateToYears(as.Date(t_vacc_end[gsub('-','_', dp)][[1]][1]))
      pomp::coef(mod)[paste0("p1d_reg", gsub('-','_', dp))] <-  as.numeric(p1d_reg[gsub('-','_', dp)][[1]][1])
      pomp::coef(mod)[paste0("r_v_year", gsub('-','_', dp))] <- as.numeric(r_v_year[gsub('-','_', dp)][[1]][1])
    }
    pomp::coef(mod)["cases_ext"] <- as.numeric(cases_ext)

    sirb_cholera <- pomp::pomp(
      mod,
      covar = pomp::covariate_table(all_rain_forecast, times = 'time')
    )

    sims <- pomp::simulate(sirb_cholera, nsim = nsim,
                           format = 'data.frame', include.data = FALSE,  # TODO: maybe make this so that it can change
                           times = time_forecast)
    # %>%
    #   select(.id, time, date, IncidenceAll, DosesAll, CasesAll, starts_with('cases'))

    sims$scenario <- scenario_str

    sims
  }


  all_sims
  # res_df <- data.frame('time' = temp[temp$isdata == 'simulation' & temp$variable == paste0("S", "Sud"), 'time'])
  #
  # important_compartments <- COMPARTMENTS[1:13]
  # comps <- paste0('(', paste(COMPARTMENTS, collapse = '|'), ')')
  #
  # temp %>%
  #   filter(!grepl('(^max)|(^rain)', variable)) %>%
  #   filter(isdata == 'simulation') %>%
  #   filter(grepl('cases', variable)) %>%
  #   mutate(
  #     dp = gsub('cases', '', variable)
  #   ) %>%
  #   select(-variable) %>%
  #   filter(time <= 2019.01) %>%
  #   ggplot(aes(x = time)) +
  #   geom_line(aes(y = q50, group = dp)) +
  #   facet_wrap(~dp) +
  #   geom_ribbon(aes(ymin = q05, ymax = q95),
  #               fill = 'red', alpha = 0.3)

  # for (dp in DEPARTMENTS) {
  #   for (comp in important_compartments) {
  #     res_df[[paste0(dp, '_q05_', comp)]] <- temp[temp$isdata == 'simulation' & temp$variable == paste0(comp, dp), ]$q05
  #     res_df[[paste0(dp, '_q50_', comp)]] <- temp[temp$isdata == 'simulation' & temp$variable == paste0(comp, dp), ]$q50
  #     res_df[[paste0(dp, '_q95_', comp)]] <- temp[temp$isdata == 'simulation' & temp$variable == paste0(comp, dp), ]$q95
  #   }
  # }
  #
  # deps <- paste0('(', paste(DEPARTMENTS, collapse = '|'), ')')
  #
  # res_df %>%
  #   tidyr::pivot_longer(
  #     data = .,
  #     cols = -c(time)
  #   ) %>%
  #   mutate(
  #     dp = gsub('_[[:alnum:]]{3}_[[:alnum:]]+$', '', name),
  #     comp = gsub(paste0(deps, '{1}_[[:alnum:]]+_'), '', name),
  #     level = gsub(paste0(deps, '_{1}|_[[:alnum:]]+$'), '', name)
  #   ) %>%
  #   select(-name) %>%
  #   ggplot(aes(x = time, y = value)) +
  #   geom_line()
}
