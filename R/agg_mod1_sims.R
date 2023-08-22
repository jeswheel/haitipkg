#' Aggregate Model 1 Simulations
#'
#' This function summarizes simulations from model 1.
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
#' @export
agg_mod1_sims <- function(sims) {

  sims |>
    dplyr::mutate(date = lubridate::ymd("2010-10-16") + lubridate::weeks(week)) |>
    dplyr::select(date, .id, cases) |>
    dplyr::group_by(date) |>
    dplyr::summarise(
      q05 = stats::quantile(cases, 0.025, na.rm = T),
      mean = mean(cases, na.rm = T),
      q50 = stats::quantile(cases, 0.5, na.rm = T),
      q95 = stats::quantile(cases, 0.975, na.rm = T)
    ) |>
    dplyr::ungroup()
}
