#' Get Susceptible Quantiles for Model 3
#'
#' This function gets quantiles for the Susceptible compartment from the saved
#' states of a block particle filter of Model 3. This function is only used in
#' the supplement material.
#'
#' @param bpf A bpfiltered object with saved states.
#' @param i the time point \eqn{t_i} to get Susceptible qunatiles.
#' @return A data.frame containing quantiles of the susceptible compartment at
#'    time t_i.
#'
#' @importFrom magrittr %>%
#' @export

get_h3_S_quants <- function(bpf, i) {
  states <- bpf@saved.states[[i]][paste0(
    rep(c("S", "VSd", "VSdd", "VSd_alt", "VSdd_alt"), each = 10), 1:10
  ),
  ]

  new_states <- as.data.frame(t(states)) %>%
    dplyr::transmute(
      S_1  = S1  + VSd1  + VSdd1  + VSd_alt1  + VSdd_alt1,
      S_2  = S2  + VSd2  + VSdd2  + VSd_alt2  + VSdd_alt2,
      S_3  = S3  + VSd3  + VSdd3  + VSd_alt3  + VSdd_alt3,
      S_4  = S4  + VSd4  + VSdd4  + VSd_alt4  + VSdd_alt4,
      S_5  = S5  + VSd5  + VSdd5  + VSd_alt5  + VSdd_alt5,
      S_6  = S6  + VSd6  + VSdd6  + VSd_alt6  + VSdd_alt6,
      S_7  = S7  + VSd7  + VSdd7  + VSd_alt7  + VSdd_alt7,
      S_8  = S8  + VSd8  + VSdd8  + VSd_alt8  + VSdd_alt8,
      S_9  = S9  + VSd9  + VSdd9  + VSd_alt9  + VSdd_alt9,
      S_10 = S10 + VSd10 + VSdd10 + VSd_alt10 + VSdd_alt10
    )

  tidyr::pivot_longer(
    new_states,
    cols = dplyr::everything(),
    names_to = 'state',
    values_to = 'S',
    names_prefix = "S_"
  ) %>%
    dplyr::mutate(dep = as.numeric(state)) %>%
    dplyr::group_by(state) %>%
    dplyr::summarize(
      Q025 = stats::quantile(S, probs = 0.025),
      Q50  = stats::quantile(S, probs = 0.5),
      Q975 = stats::quantile(S, probs = 0.975)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      Q025 = Q025 / bpf@params[paste0("H", state)],
      Q50  = Q50  / bpf@params[paste0("H", state)],
      Q975 = Q975 / bpf@params[paste0("H", state)],
      dep = unit_names(bpf)[as.numeric(state)],
      time = bpf@times[i]
    ) %>%
    dplyr::select(time, dep, Q025, Q50, Q975)
}
