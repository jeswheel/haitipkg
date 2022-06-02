#' Aggregate Model 3 Simulations
#'
#' Because Model 3 has spatial structure, simulations from model 3 result in
#' simulations of the process for each departement. For comparitive purposes,
#' it's useful to aggregate these results into national counts, and this
#' function does precisly that.
#'
#' @param sims data.frame containing simulations from model 3.
#'
#' @return A data.frame with columns: \itemize{
#'   \item{time: }{Numerical year}
#'   \item{q05: }{2.5 percentile of reported cholera cases}
#'   \item{mean: }{Mean of reported cholera cases}
#'   \item{q50: }{Median of reported cholera cases}
#'   \item{q95: }{97.5 percentile of reported cholera cases}
#' }
#'
#' @importFrom magrittr %>%
#' @export
agg_mod3_sims <- function(sims) {

  # First Define some helper functions:
  dateToYears <- function(date, origin = as.Date("2014-01-01"), yr_offset = 2014) {
    # This function converts a date to a decimal representation
    #
    # ex: "1976-03-01" -> 1976.163

    julian(date, origin = origin) / 365.25 + yr_offset
  }

  yearsToDate <- function(year_frac, origin = as.Date("2014-01-01"), yr_offset = 2014.0) {
    # This function is the inverse function of dateToYears; it takes
    # a decimal representation of a date and converts it into a Date.
    #
    # ex: 1976.163 -> "1976-03-01"

    as.Date((year_frac - yr_offset) * 365.25, origin = origin)
  }

  yearsToDateTime <- function(year_frac, origin = as.Date("2014-01-01"), yr_offset = 2014.0) {
    # Same as the function above, but a DateTime object rather than a Date
    # object.
    #
    # ex: 1976.163 -> "1976-03-01"
    as.POSIXct((year_frac - yr_offset) * 365.25 * 3600 * 24, origin = origin)
  }

  sims %>%
    tidyr::pivot_wider(  # Move rows to columns for department results
      data = .,
      id_cols = c(time, .id),
      names_from = unitname,
      values_from = c(cases, totInc)
    ) %>%
    dplyr::mutate(  # Sum all reported cholera infections
      ReportedAll = cases_Artibonite + cases_Centre +
        cases_Grande_Anse + cases_Nippes + cases_Nord +
        cases_Nord_Est + cases_Ouest + cases_Sud +
        cases_Sud_Est + cases_Nord_Ouest
    ) %>%
    dplyr::select(time, .id, ReportedAll) %>%
    dplyr::group_by(time) %>%
    dplyr::summarise(q05 = quantile(ReportedAll, 0.025, na.rm = T),
              mean = mean(ReportedAll, na.rm = T),
              q50 = quantile(ReportedAll, 0.5, na.rm = T),
              q95 = quantile(ReportedAll, 0.975, na.rm = T)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(date = yearsToDateTime(time)) %>%
    dplyr::mutate(date = as.Date(lubridate::round_date(date)))
}
