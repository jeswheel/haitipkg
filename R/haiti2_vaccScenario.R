#' Model 2: Project Vaccination Scenario
#'
#' This function projects the following vaccination scenarios under Model 2:
#'   \describe{
#'     \item{V0}{The default vaccination scenario: no additional vaccinations are deployed.}
#'     \item{V1}{Two department (Centre and Artibonite) over two years.}
#'     \item{V2}{Three department (Centre, Artibonite, and Ouest) over two years.}
#'     \item{V3}{Countrywide vaccination over 5 years.}
#'     \item{V4}{Countrywide vaccination over 2 years.}
#'   }
#'
#' @param h2_params: a vector of parameters for the model
#'
#' @importFrom magrittr %>%
#' @return List containing two objects: \describe{
#'   \item{tot_inc}{Vector estimating the total cholera incidence from Feb 2019-Feb 2024 for each scenario.}
#'   \item{all_trajs}{data.frame containing trajectories for each of the vaccination scenarios.}
#' }
#'
#' @examples
#' h2 <- fit_haiti2()
#' results <- haiti2_vaccScenario(h2_params = h2$h2_params)
#'
#' @export
haiti2_vaccScenario <- function(h2_params) {

  # Create joint model to simulate the trajectory
  h2 <- haiti2(cutoff = 10000, measure = "log")

  if (!setequal(names(h2_params), names(pomp::coef(h2)))) {
    stop("Names of input parameters do not match the model.")
  }

  # Set model coefficients to input values
  pomp::coef(h2) <- h2_params

  # Get epidemic trajectory
  epi_traj <- pomp::trajectory(
    h2,
    params = pomp::coef(h2),
    format = 'data.frame',
    times = pomp::time(h2)
  )

  # # Get endemic trajectory
  # end_traj <- pomp::trajectory(
  #   h2_end,
  #   params = pomp::coef(h2_end),
  #   format = 'data.frame'
  # )

  # Get reported cases
  epi_traj$Ctotal <- rowSums(epi_traj[, paste0("C", 1:10)]) * h2_params['Rho']
  # end_traj$Ctotal <- rowSums(end_traj[, paste0("C", 1:10)]) * end_params['Rho']

  # Create vector of forcast times
  time_forecast <- lubridate::decimal_date(
    seq(
      as.Date(lubridate::round_date(lubridate::date_decimal(min(h2@times)), unit = 'day')),
      lubridate::ymd("2029-12-21"),
      by = "1 week")
  )

  # Set vacination scenario parameters to endemic model parameters
  V0_params <- h2_params
  V1_params <- h2_params
  V2_params <- h2_params
  V3_params <- h2_params
  V4_params <- h2_params
  V1_params['scenario'] <- 1
  V2_params['scenario'] <- 2
  V3_params['scenario'] <- 3
  V4_params['scenario'] <- 4

  # Project under V0
  V0_traj <- pomp::trajectory(
    h2,
    params = V0_params,
    times = time_forecast,
    format = "data.frame"
  ) %>%
    dplyr::filter(
     year >= max(h2@times)
    )
  V0_traj$scenario <- "V0"
  V0_traj$Itotal <- rowSums(V0_traj[, paste0("C", 1:10)])
  V0_traj$Atotal <- (1 - V0_params['k']) / (V0_params['k']) * V0_traj$Itotal
  V0_traj$totInc <- V0_traj$Itotal + V0_traj$Atotal
  V0_traj$ReportedCases <- V0_traj$Itotal * V0_params["Rho"]

  # Project under V1
  V1_traj <- pomp::trajectory(
    h2,
    params = V1_params,
    times = time_forecast,
    format = "data.frame"
  ) %>%
    dplyr::filter(
      year >= max(h2@times)
    )
  V1_traj$scenario <- "V1"
  V1_traj$Itotal <- rowSums(V1_traj[, paste0("C", 1:10)])
  V1_traj$Atotal <- (1 - V1_params['k']) / (V1_params['k']) * V1_traj$Itotal
  V1_traj$totInc <- V1_traj$Itotal + V1_traj$Atotal
  V1_traj$ReportedCases <- V1_traj$Itotal * V1_params["Rho"]

  # Project under V2
  V2_traj <- pomp::trajectory(
    h2,
    params = V2_params,
    times = time_forecast,
    format = "data.frame"
  ) %>%
    dplyr::filter(
      year >= max(h2@times)
    )
  V2_traj$scenario <- "V2"
  V2_traj$Itotal <- rowSums(V2_traj[, paste0("C", 1:10)])
  V2_traj$Atotal <- (1 - V2_params['k']) / (V2_params['k']) * V2_traj$Itotal
  V2_traj$totInc <- V2_traj$Itotal + V2_traj$Atotal
  V2_traj$ReportedCases <- V2_traj$Itotal * V2_params["Rho"]

  # Project under V3
  V3_traj <- pomp::trajectory(
    h2,
    params = V3_params,
    times = time_forecast,
    format = "data.frame"
  ) %>%
    dplyr::filter(
      year >= max(h2@times)
    )
  V3_traj$scenario <- "V3"
  V3_traj$Itotal <- rowSums(V3_traj[, paste0("C", 1:10)])
  V3_traj$Atotal <- (1 - V3_params['k']) / (V3_params['k']) * V3_traj$Itotal
  V3_traj$totInc <- V3_traj$Itotal + V3_traj$Atotal
  V3_traj$ReportedCases <- V3_traj$Itotal * V3_params["Rho"]

  # Project under V4
  V4_traj <- pomp::trajectory(
    h2,
    params = V4_params,
    times = time_forecast,
    format = "data.frame"
  ) %>%
    dplyr::filter(
      year >= max(h2@times)
    )
  V4_traj$scenario <- "V4"
  V4_traj$Itotal <- rowSums(V4_traj[, paste0("C", 1:10)])
  V4_traj$Atotal <- (1 - V4_params['k']) / (V4_params['k']) * V4_traj$Itotal
  V4_traj$totInc <- V4_traj$Itotal + V4_traj$Atotal
  V4_traj$ReportedCases <- V4_traj$Itotal * V4_params["Rho"]

  # Get the interesting results in a single matrix
  mod2_all_trajs <- rbind(
    V0_traj[, c('scenario', 'year', 'Itotal', 'Atotal', 'totInc', 'ReportedCases')],
    V1_traj[, c('scenario', 'year', 'Itotal', 'Atotal', 'totInc', 'ReportedCases')],
    V2_traj[, c('scenario', 'year', 'Itotal', 'Atotal', 'totInc', 'ReportedCases')],
    V3_traj[, c('scenario', 'year', 'Itotal', 'Atotal', 'totInc', 'ReportedCases')],
    V4_traj[, c('scenario', 'year', 'Itotal', 'Atotal', 'totInc', 'ReportedCases')]
  ) %>%
    as.data.frame() %>%  # Convert to data.frame
    dplyr::mutate(
      date = lubridate::date_decimal(year),
      date = lubridate::round_date(date, unit = 'day')  # Convert to type Date
    )

  # Create total incidence vector
  mod2_tot_Inc <- c(
    "V0" = round(
      sum(
        dplyr::filter(
          V0_traj,
          year >= lubridate::decimal_date(as.Date("2019-02-01")),
          year <= lubridate::decimal_date(as.Date("2024-02-01")))$totInc
      )
    ),
    "V1" = round(
      sum(
        dplyr::filter(
          V1_traj,
          year >= lubridate::decimal_date(as.Date("2019-02-01")),
          year <= lubridate::decimal_date(as.Date("2024-02-01")))$totInc
      )
    ),
    "V2" = round(
      sum(
        dplyr::filter(
          V2_traj,
          year >= lubridate::decimal_date(as.Date("2019-02-01")),
          year <= lubridate::decimal_date(as.Date("2024-02-01")))$totInc
      )
    ),
    "V3" = round(
      sum(
        dplyr::filter(
          V3_traj,
          year >= lubridate::decimal_date(as.Date("2019-02-01")),
          year <= lubridate::decimal_date(as.Date("2024-02-01")))$totInc
      )
    ),
    "V4" = round(
      sum(
        dplyr::filter(
          V4_traj,
          year >= lubridate::decimal_date(as.Date("2019-02-01")),
          year <= lubridate::decimal_date(as.Date("2024-02-01")))$totInc
      )
    )
  )

  # Create list to return results
  results <- list()
  results[['tot_inc']] <- mod2_tot_Inc
  results[['all_trajs']] <- mod2_all_trajs

  return(results)
}
