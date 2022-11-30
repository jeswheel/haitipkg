#' Project Rainfall
#'
#' This function projects rainfall in each of the Haitian departements. The
#' rainfall is used as a covariate in Model 3.
#'
#' The rainfall projection is made by assuming that the rainfall on certain
#' dates in the future will be similar to the rainfall on the same dates in
#' previous years. Therefore the projection is made by the following steps:
#' \enumerate{
#'    \item Starting from 1 day after the available data (either actual data
#'    or projected rainfall from previous steps), get the 14-day window of dates
#'    in which rainfall will be projected.
#'    \item With this 14-day window, randomly select a year in the past in which
#'    actual rainfall data is available.
#'    \item Project the rainfall on the selected 14-day window by using the
#'    observed rainfall measurements over the same dates in the randomly
#'    selected year.
#'    \item Repeat the process until a projection for all desired dates in the
#'    future is obtained.
#' }
#'
#' Note that this forecasting procedure is stochastic and that a different
#' rainfall projection will be obtained each time the function is called. While
#' it is unlikely that this function will provide precise estimates of future
#' rainfall, it will likely capture seasonal trends in rainfall that may affect
#' the spread of Cholera in Haiti.
#'
#' @param end_date: Date object representing when the projection of rainfall
#'  should stop.
#'
#' @importFrom magrittr %>%
#' @export

project_rain <- function(end_date = as.Date("2029-12-20"),
                         include_data = FALSE) {

  departments <- c(
    'Artibonite', 'Centre', 'Grande_Anse', 'Nippes', 'Nord',
    'Nord-Est', 'Nord-Ouest', 'Ouest', 'Sud', 'Sud-Est'
  )

  rf <- haitiRainfall %>%
    dplyr::filter(
      date >= as.Date("2010-10-23") - lubridate::days(8)
    )

  num_days <- 14
  dti <- dplyr::pull(rf[1, 'date'])
  dtf <- dplyr::pull(rf[nrow(rf), 'date'])
  rf_prj_index <- seq(dtf + 1, end_date, '1 day')

  rain_prj <- matrix(0, nrow = length(rf_prj_index),
                     ncol = 10)

  yrs <- seq(lubridate::year(dti) + 1, lubridate::year(dtf) - 2)
  date_range <- seq(dtf + 1, end_date, '14 day')  # 14 from num_days

  for (i in 1:length(date_range)) {
    dd <- lubridate::day(date_range[i])

    # Deal with leap years, as was done in original code.
    if (lubridate::month(date_range[i]) == 2 && dd == 29) {
      dd <- 28
    }

    # Pick random year, but use the same month and same day
    pick <- as.Date(
      ISOdate(
        year = sample(yrs, 1),  # This is where the randomness comes from.
        month = lubridate::month(date_range[i]),
        day = lubridate::day(date_range[i])
      )
    )

    # Ensure we are not sampling rainfall data from hurricane Matthew!
    while (lubridate::month(pick) == 10 && lubridate::year(pick) == 2016 && lubridate::day(pick) < 14) {
      # Pick random year, but use the same month and same day
      pick <- as.Date(
        ISOdate(
          year = sample(yrs, 1),  # This is where the randomness comes from.
          month = lubridate::month(date_range[i]),
          day = lubridate::day(date_range[i])
        )
      )
    }

    # After getting random day, get the 2 weeks following that day.
    pick_seq <- seq(
      pick, pick + num_days - 1, 'day'
    )

    # Get rain values for entire country during those two weeks.
    replacement <- as.matrix(rf[as.Date(rf$date) %in% pick_seq, -1])
    colnames(replacement) <- NULL
    rownames(replacement) <- NULL
    rain_prj[(num_days * (i - 1) + 1):(num_days * i), ] <- replacement
  }

  rain_prj <- data.frame(
    date = rf_prj_index, rain_prj
  )

  colnames(rain_prj) <- c('date', departments)

  if (include_data) {
    rf$is_data <- TRUE
    rain_prj$is_data <- FALSE
    all_rain <- dplyr::bind_rows(rf, rain_prj)
  } else {
    all_rain <- rain_prj
  }

  # The values were determined to be the same transformation that was used on
  # the "training" data, i.e., dividing by the maximum rainfall event in
  # each department from Oct 2010 to Jan 2019.
  result <- all_rain %>%
    dplyr::filter(date >= as.Date("2010-10-23") - 4) %>%
    mutate(
      Artibonite = Artibonite / 121.71200,
      Centre = Centre / 105.10750,
      Grande_Anse = Grande_Anse / 167.86139,
      Nippes = Nippes / 179.57455,
      Nord = Nord / 96.05684,
      `Nord-Est` = `Nord-Est` / 71.97000,
      `Nord-Ouest` = `Nord-Ouest` / 130.11900,
      Ouest = Ouest / 185.94521,
      Sud = Sud / 200.50558,
      `Sud-Est` = `Sud-Est` / 225.66032,
      time = lubridate::decimal_date(date)
    )

  colnames(result) <- c(
    "date",
    paste0(
      'rain_std', c(
        'Artibonite', 'Centre', 'Grande_Anse',
        'Nippes', 'Nord', 'Nord-Est', 'Nord-Ouest',
        'Ouest', 'Sud', 'Sud-Est'
      )
    ),
    'time'
  )

  result %>% dplyr::select(time, dplyr::starts_with("rain_std")) %>%
    as.matrix()
}
