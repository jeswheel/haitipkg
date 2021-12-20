#' Vaccination Scenarios
#'
#' \code{VaccinationScenario3} returns lists containing various vaccination scenarios.
#'
#' @param course_year Number of years that the vaccination campaign is conducted.
#' @param percent_uv Percent of completly unvaccinated individuals.
#' @param percent_1d Percent of individuals that received a single dose of vaccination.
#' @param percent_2d Percent of individuals that received two doses of the vaccination.
#' @param ve \code{c(1, 2, 3)}. Vaccination Efficacy. Only 1 is used by original authors.
#' @param not_dep Departments that are not part of the vaccination campaign.
#'
#' @return List containing parameters of the vaccination campaign.
#' @export

VaccinationScenario3 <- function(
  course_year, percent_uv, percent_1d, percent_2d, ve = 1, not_dep = c()
) {

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

  DEPARTMENTS <- c(
    'Artibonite', 'Centre', 'Grande_Anse', 'Nippes', 'Nord',
    'Nord_Est', 'Nord_Ouest', 'Ouest', 'Sud', 'Sud_Est'
  )

  # Population of each of the departments. Matches the values you get when you
  # just google the department populations, except Nippes, which says 342325.
  pop <- c(
    1727524, 746236, 468301, 342525, 1067177,
    393967, 728807, 4029705, 774976, 632601
  )

  names(pop) <- DEPARTMENTS

  # Order that the departments receive vaccinations.
  ocv_order <- c(
    'Centre', 'Artibonite', 'Ouest', 'Nord_Ouest', 'Nord', 'Sud', 'Nippes', 'Nord_Est', 'Sud_Est',
    'Grande_Anse'
  )

  # Create a list that has a value for each of the departements.
  t_vacc_start <- list()  # When vaccination starts
  t_vacc_end <- list()  # When vaccination ends
  p1d_reg <- list()  # TODO: I think this is percent of the vaccinations that are only one dose. Verify.
  r_v_year <- list()  # TODO: I think this is annual vaccination rate. Verify.

  # 20% completely unvaccinated, 10% one-dose only, 70% two doses
  t_init <- as.Date(
    ISOdate(
      year = 2019,
      month = 1,
      day = 12
    )
  )

  days_per_department <- round((course_year * 365) / length(ocv_order))

  for (i in 1:length(ocv_order)) {
    dp <- ocv_order[i]
    if (!( dp %in% not_dep )) {
      t_vacc_start[[dp]] <- as.character(t_init + (i * days_per_department))
      t_vacc_end[[dp]] <- as.character(t_init + (i + 1) * days_per_department)
      p1d_reg[[dp]] <- percent_1d / (percent_1d + percent_2d)
      r_v_year[[dp]] <- pop[dp] * (100 - percent_uv) / 100 / days_per_department * 365.25
    } else {
      t_vacc_start[[dp]] <- as.character(t_init + (i * days_per_department))
      t_vacc_end[[dp]] <- as.character(t_init + (i + 1) * days_per_department)
      p1d_reg[[dp]] <- 0
      r_v_year[[dp]] <- 0
    }
  }

  ret <- list()
  ret$t_vacc_start <- t_vacc_start
  ret$t_vacc_end <- t_vacc_end
  ret$p1d_reg <- p1d_reg
  ret$r_v_year <- r_v_year
  ret$ve <- ve
  ret$not_dep <- not_dep

  ret
}


