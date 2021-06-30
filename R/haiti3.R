#' Build pomp object for model 3.
#'
#'
#'
#'
#'
#'

haiti3 <- function() {

  library(magrittr)
  library(foreach)

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

  # List all state-names in pomp object:
  state_names <- c(
    "S", "I", "A", "RI1", "RI2", "RI3", "RA1", "RA2", "RA3",
    "VSd", "VRI1d", "VRI2d", "VRI3d", "VRA1d", "VRA2d", "VRA3d",
    "VSdd", "VRI1dd", "VRI2dd", "VRI3dd", "VRA1dd", "VRA2dd",
    "VRA3dd", "VSd_alt", "VRI1d_alt", "VRI2d_alt", "VRI3d_alt",
    "VRA1d_alt", "VRA2d_alt", "VRA3d_alt", "VSdd_alt", "VRI1dd_alt",
    "VRI2dd_alt", "VRI3dd_alt", "VRA1dd_alt", "VRA2dd_alt",
    "VRA3dd_alt", "B", "C", "W"
  )

  # All parameters that are common to each departement (disease specific parameters)
  params_common <- c(
    "sigma", "mu_B", "thetaI", "XthetaA", "lambdaR", "r",
    "gammaI", "gammaA", "rhoA", "XrhoI", "epsilon", "k",
    "std_W", "cas_def", "Rtot_0", "mu", "alpha", "cases_ext"
  )

  # Parameters that are unique to each department:
  params_diff <- c(
    "foi_add", "betaB", "H", "D","t_vacc_start",
    "t_vacc_end", "p1d_reg", "r_v_year", "t_vacc_start_alt",
    "t_vacc_end_alt", "p1d_reg_alt", "r_v_year_alt"
  )

  # Loads the input parameters
  load('R/sysdata.rda')
  t_start <- dateToYears(as.Date(input_parameters$t_start))
  t_end   <- dateToYears(as.Date(input_parameters$t_end))

  all_state_names <- c('IncidenceAll', 'DosesAll', 'CasesAll')
  all_param_names <- params_common

  all_matrix_cases_at_t_start.string <- ""
  all_matrix_cases_other.string <- ""

  all_cases <- ALL_CASES %>%
    dplyr::mutate(
      date = as.Date(date, format = '%Y-%m-%d'),
      time = dateToYears(date)
    )

  all_rain <- ALL_RAIN %>%
    dplyr::mutate(
      date = as.Date(date, format = "%Y-%m-%d"),
      time = dateToYears(date)
    ) %>%
    dplyr::filter(time > t_start - 0.01 & time < (t_end + 0.01))

  for (dp in departements) {

    cases <- ALL_CASES %>%
      tidyr::gather(dep, cases, -date) %>%
      dplyr::filter(dep == dp) %>%
      dplyr::mutate(
        date = as.Date(date, format = "%Y-%m-%d"),
        time = dateToYears(date)
      )

    rain <- ALL_RAIN %>%
      tidyr::gather(dep, rain,-date) %>%
      dplyr::group_by(dep) %>%
      dplyr::ungroup() %>%
      dplyr::filter(dep == dp) %>%
      dplyr::mutate(
        date = as.Date(date, format = "%Y-%m-%d"),
        time = dateToYears(date)
      ) %>%
      dplyr::filter(time > t_start - 0.01 & time < (t_end + 0.01)) %>%
      dplyr::mutate(max_rain = max(rain), rain_std = rain / max_rain)

    all_rain <- cbind(all_rain, placeholder = rain$max_rain)
    all_rain <- cbind(all_rain, placeholder2 = rain$rain_std)
    names(all_rain)[names(all_rain) == "placeholder"] <- paste0('max_rain', gsub('-','_',dp))
    names(all_rain)[names(all_rain) == "placeholder2"] <- paste0('rain_std', gsub('-','_',dp))

    all_cases <- cbind(all_cases, placeholder = cases$cases)
    names(all_cases)[names(all_cases) == "placeholder"] <- paste0('cases', gsub('-','_',dp))

    cases_other_dept <- ALL_CASES  %>%
      tidyr::gather(dep, cases,-date) %>%
      dplyr::filter(dep != dp) %>%
      dplyr::mutate(
        date = as.Date(date, format = "%Y-%m-%d"),
        time = dateToYears(date)
      )

    cases_other_dept <- aggregate(
        cases_other_dept$cases,
        by = list(Category = cases_other_dept$time),
        FUN = sum,
        na.rm = TRUE,
        na.action = NULL
      ) %>%
      dplyr::mutate(time = Category) %>%
      dplyr::mutate(cases = x)


    cases_at_t_start <- cases %>% dplyr::filter(dateToYears(date) <= t_start)
    cases_at_t_start.string <- foreach::foreach(r = iterators::iter(cases_at_t_start, by = "row"),
                                                .combine = c) %do% {
                                                  sprintf(" {%f, %f} ", r$time, r$cases)
                                                } %>%
      stringr::str_c(collapse = ", \n")

    matrix_cases_at_t_start.string <- stringr::str_c(
      sprintf(
        "double cases_at_t_start%s[%i][%i] = {\n",
        gsub('-', '_', dp),
        nrow(cases_at_t_start),
        2
      ),
      cases_at_t_start.string,
      " \n };"
    )

    # Cases from other departement as a mobility rational
    cases_other.string <-
      foreach::foreach(r = iterators::iter(cases_other_dept, by = "row"),
                       .combine = c) %do% {
                         sprintf(" {%f, %f} ", r$time, r$cases)
                       } %>%
      stringr::str_c(collapse = ", \n")

    matrix_cases_other.string <-
      stringr::str_c(sprintf("double cases_other%s[%i][%i] = {\n", gsub('-', '_', dp), nrow(cases_other_dept), 2),
            cases_other.string,
            " \n };")

    all_matrix_cases_at_t_start.string = stringr::str_c(all_matrix_cases_at_t_start.string, matrix_cases_at_t_start.string)
    all_matrix_cases_other.string = stringr::str_c(all_matrix_cases_other.string, matrix_cases_other.string)

    all_state_names = append(all_state_names, lapply(state_names, paste0, gsub('-', '_',dp)))
    all_param_names = append(all_param_names, lapply(params_diff, paste0, gsub('-', '_',dp)))

    all_state_names = unlist(all_state_names)
    all_param_names = unlist(all_param_names)

    all_params <- purrr::set_names(seq_along(all_param_names) * 0, all_param_names)

  }


  # TODO: Continue starting with line 286 in pomp_all_dept.R

}
