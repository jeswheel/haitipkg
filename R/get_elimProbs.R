#' Get Elimination Probabilities
#'
#' This function takes simulations from a model as input, and outputs a list
#' containing the probability of elimination, Cumulative Infections,
#' and probability of elimination over time.
#'
#' The probability of elimination is calculated as:
#' \itemize{
#'   \item{Probability of Elimination}{The proportion of simulations that
#'   achieved less than one infection with V \emph{cholerae} (including reported
#'   and unreported infections) over at least 52 consecutive weeks.}
#'   \item{Cumulative Infections}{Median estimate with 95% CI of cumulative
#'   infections from February 2019, to February 2024.}
#'   \item{Probability of elimination over time}{he proportion of simulations that
#'   achieved less than one infection with V \emph{cholerae} (including reported
#'   and unreported infections) over at least 52 consecutive weeks, over time.}
#' }
#'
#' @param sims data.frame containing simulations from one of the vaccination scenarios
#' @param model numeric in c(1, 3) indicating which model the simulations are from.
#'
#' @return A list containing the Probability of Elimination (probElim),
#' Cumulative Infections (cumInf), and Probability of Elimination over time
#' (ElimTime).
#'
#' @import dplyr
#' @export
get_elimProbs <- function(sims, model) {

  if (model == 3) {
    # Here, we need to add cases together to get observed cases.

    yearsToDate <- function(year_frac, origin = as.Date("2014-01-01"), yr_offset = 2014.0) {
      # This function is the inverse function of dateToYears; it takes
      # a decimal representation of a date and converts it into a Date.
      #
      # ex: 1976.163 -> "1976-03-01"

      as.Date((year_frac - yr_offset) * 365.25, origin = origin)
    }


    sims <- sims %>%
      tidyr::pivot_wider(
        id_cols = c(time, .id),
        names_from = unitname,
        values_from = c(cases, totInc)
      ) %>%
      mutate(
        ReportedAll = cases_Artibonite + cases_Centre +
          cases_Grande_Anse + cases_Nippes + cases_Nord +
          cases_Nord_Est + cases_Ouest + cases_Sud +
          cases_Sud_Est + cases_Nord_Ouest,
        IncidenceAll = totInc_Artibonite + totInc_Centre +
          totInc_Grande_Anse + totInc_Nippes + totInc_Nord +
          totInc_Nord_Est + totInc_Ouest + totInc_Sud +
          totInc_Sud_Est + totInc_Nord_Ouest
      ) %>%
      select(time, .id, ReportedAll, IncidenceAll) %>%
      group_by(.id) %>%
      mutate(
        cumIncidence = IncidenceAll,
        cumReported = cumsum(ReportedAll),
        IncidenceIncrease52 = cumIncidence - lag(cumIncidence, 52),
        ReportedIncrease52 = cumReported - lag(cumReported, 52)
      ) %>%
      ungroup()

    # # here, cases<Departement> represents the simulated reported cases in
    # # <Departement>
    # #
    # # Note that CasesAll is total number of Infectious Cases, and IncidenceAll
    # # is total number of all infections.
    # sims <- sims %>%
    #   mutate(
    #     ReportedAll = casesArtibonite + casesCentre +
    #       casesGrande_Anse + casesNippes + casesNord +
    #       casesNord_Est + casesOuest + casesSud +
    #       casesSud_Est + casesNord_Ouest,
    #     time = yearsToDate(time)
    #   ) %>%
    #   select(time, .id, ReportedAll, CasesAll, IncidenceAll, DosesAll) %>%
    #   group_by(.id) %>%
    #   mutate(
    #     cumIncidence = cumsum(IncidenceAll),  # Get total sum of incidence (A + I)
    #     cumReported = cumsum(ReportedAll),  # Get total sum of reported cases
    #     # cumCasesAll = cumsum(CasesAll),  # Get total sum of symptomatic infections (I)
    #     IncidenceIncrease52 = cumIncidence - lag(cumIncidence, 52),  # See if cases have increased in last 52 weeks.
    #     # CasesIncrease52 = cumCasesAll - lag(cumCasesAll, 52),
    #     ReportedIncrease52 = cumReported - lag(cumReported, 52)
    #   ) %>%
    #   ungroup()

    gc()
  } else if (model == 1) {
    sims <- sims %>%
      mutate(time = lubridate::ymd("2010-10-16") + lubridate::weeks(week)) %>%
      select(time, .id, cases, incid) %>%
      group_by(.id) %>%
      mutate(
        cumIncidence = cumsum(incid),  # Get total sum of incidence (A + I)
        cumReported = cumsum(cases),  # Get total sum of reported cases
        IncidenceIncrease52 = cumIncidence - lag(cumIncidence, 52),  # See if cases have increased in last 52 weeks.
        ReportedIncrease52 = cumReported - lag(cumReported, 52)
      ) %>%
      ungroup()
  } else if (model == 2){
    stop("Model 2 Not Yet Implemented")
  } else {
    stop("Model must be an integer value in c(1, 2, 3).")
  }

  probElim <- sims %>%
    group_by(.id) %>%
    filter(0 %in% IncidenceIncrease52) %>%  # Find simulations when Cumulative incidence didn't increase for 52 weeks
    ungroup() %>%
    summarize(n_distinct(.id) / length(unique(sims$.id))) %>%
    pull()

  # Now I want to know if any resurgence after elimination occured:
  # results %>%
  #   group_by(.id) %>%
  #   filter(0 %in% IncidenceIncrease52) %>%  # Only look at groups were elimination occured
  #   slice(min(which((IncidenceIncrease52 == 0))):n()) %>% # Only look at data after elimination occured
  #   filter(IncidenceIncrease52 > 0) # Only look at groups that had cases after elimination

  # Cumulative Number of infections after 5 years
  cumInf <- sims %>%
    filter(time >= "2019-02-01", time <= "2024-02-01") %>%
    group_by(.id) %>%
    summarize(fiveYearIncrease = cumIncidence[time == max(time)] - cumIncidence[time == min(time)]) %>%
    ungroup() %>%
    summarize(q025 = stats::quantile(fiveYearIncrease, probs = 0.025),
              q50 = stats::quantile(fiveYearIncrease, probs = 0.5),
              q975 = stats::quantile(fiveYearIncrease, probs = 0.975))

  sims0 <- sims %>%
    group_by(.id) %>%
    filter(0 %in% IncidenceIncrease52)

  if (nrow(sims0) == 0) {
    results_df <- data.frame('time' = unique(sims$time))
    results_df$elim_prob = 0
  } else {
    ElimTime <- sims0 %>%  # Find only sims that eliminated cholera
      slice(min(which((IncidenceIncrease52 == 0)))) %>%  # Find when cholera was eliminated
      ungroup() %>%
      count(time) %>%  # get count of elimination times
      arrange(time) %>%  # Important because of cumsum below
      mutate(elim_prob = cumsum(n) / length(unique(sims$.id))) %>%
      select(-n)

    results_df <- data.frame('time' = unique(sims$time))

    results_df <- dplyr::left_join(results_df, ElimTime, by = 'time') %>%
      tidyr::fill(elim_prob, .direction = 'down')
  }

  final <- list()
  final[['probElim']] <- probElim
  final[['cumInf']] <- cumInf
  final[['ElimTime']] <- results_df

  final

}
