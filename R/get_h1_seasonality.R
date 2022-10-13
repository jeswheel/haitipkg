#' Get Haiti 1 Seasonality
#'
#' This function takes as input a pomp representation of Model 1 that contains
#' estimated parameter values, and returns a data.frame object with the
#' estimated seasonality over time.
#'
#' @param mod1 a pomp object representing Model 1. It must contain parameters.
#' @param times a vector containing the weeks that will be used to estimate
#'    the seasonality. By default, this is chosen as 431:800, which starts the
#'    week after available training data and covers at least the entirety of
#'    the following year. In this case, the linear trend in beta ends after the
#'    available training data.
#'
#' @return data.frame containing three columns: (1) week, which is week number
#'    corresponding to the measurement time (2) date, which is the date of the
#'    measurement and (3) trans_std, which contains the value of the Model 1
#'    seasonality component.
#'
#' @import dplyr
#' @export
get_h1_seasonality <- function(mod1, times = 431:800) {

  params <- coef(mod1)
  betas <- params[grepl("^beta[[:digit:]]+$", names(params))]
  betat <- params[grepl('^betat$', names(params))]

  covar_df <- as.data.frame(t(mod1@covar@table))
  covar_df$times <- mod1@covar@times

  results_df <- data.frame(
    'week' = times,
    'date' = lubridate::ymd("2010-10-16") + lubridate::weeks(times),
    'trans' = NA_real_
  )

  for (t in times) {
    si <- unlist(covar_df[covar_df$times == t, -7])

    if (t <= 430) {
      val <- as.numeric(exp(betas %*% si + betat * ((t - 215) / (430-215))))
    } else {
      val <- as.numeric(exp(betas %*% si + betat))
    }

    results_df[results_df$week == t, 'trans'] <- val
  }

  results_df %>%
    mutate(year = lubridate::year(date)) %>%
    filter(year == 2020) %>%
    mutate(trans_std = (trans - min(trans)) / (max(trans) - min(trans))) %>%
    select(week, date, trans_std)
}
