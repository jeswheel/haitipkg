#' Project from Filtering Distribution
#'
#' This function simulates from a Haiti-cholera POMP model into the future
#' starting from filtering distribution at the the last time point available.
#'
#' Because we are primarily interested in projecting the number of Cholera
#' cases in the future under various vaccination scenarios, it makes sense to
#' use as much information as possible, as simulating from the start of the
#' time series will ultimately disregard the known truth at the end of the
#' available date. The state of the system at the end of the time series
#' can be accounted for by simulating the filtering distribution:
#'
#' \deqn{Xk | Y1 = y1*, \ldots, Yk*}
#'
#' As currently implemented, the efficiency could be greatly improved for
#' projecting models that do not have covariates.
#'
#' @param mod: POMP model that you would like to simulate from.
#' @param PF: pfilter object, the result of calling
#'    pfilter(..., save.states = TRUE). Note that the states must be saved in
#'    this object for the function to work.
#' @param covarGen: function that is used to generate covariates, if the model
#'   needs covariates to run rprocess.
#' @param nsims: Integer number of simulations desired.
#' @param seed: seed for the random number generator.
#' @param cores: Integer number of cores to be used. If this number is
#'    greater than 1, the simulations will be done in parallel.
#'
#' @importFrom magrittr %>%
#' @importFrom foreach %dopar%
#' @export

project_from_filter <- function(mod, PF, covarGen = NULL,
                                nsims = 100, seed = 123, cores = 1) {

  set.seed(seed)
  max_cores <- parallel::detectCores()

  # Define helper functions
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

  std_rain <- function(x) {
    # This function simply standardizes the rain for us.
    x / max(x)
  }

  # Check if cores > 1, if so, register the cores and
  if (cores != 1 && cores > 0) {
    doParallel::registerDoParallel(as.integer(min(max_cores, cores)))
    doRNG::registerDoRNG(sample(1:99999, 1))
  }

  # Check if there are saved states, otherwise everything below fails
  if (length(saved.states(PF)) == 0) {
    stop("pfilter object must have saved.states.")
  }

  # Get the saved states from filtering distribution.
  ss <- pomp::saved.states(PF)

  end_states <- as.data.frame(t(ss[[length(ss)]]))
  state_names <- colnames(end_states)

  # Save a vector of the times we want to simulate in the future. We
  # start one week after the last day of available data, because the
  # last day of available data becomes t0.
  new_times <- dateToYears(
    seq.Date(
      yearsToDate(max(mod@times)) + 7,
      as.Date("2029-12-20"),
      by = "1 week"
    )
  )

  foreach::foreach(
    i = 1:nsims,
    .combine = dplyr::bind_rows
  ) %dopar% {

    # Sample from end states
    samp_end_state <- end_states[sample(nrow(end_states), 1), ]

    if (!is.null(covarGen)) {
      # Generate covariates
      cov_df <- covarGen()

      # Add covariates to the model, and then run rprocess from end-state:
      proc_sim <- mod %>%
        pomp::pomp(
          covar = pomp::covariate_table(cov_df, times = 'time')
        ) %>%
        pomp::rprocess(
          .,
          x0 = samp_end_state,
          t0 = max(mod@times),
          times = new_times,
          params = mod@params
        )

      measures <- mod %>%
        pomp::pomp(
          covar = pomp::covariate_table(cov_df, times = 'time')
        ) %>%
        rmeasure(
          object = .,
          x = proc_sim,
          times = new_times,
          params = mod@params
        ) %>%
        drop() %>% t() %>% as.data.frame()
    } else {

      # run rprocess from end-state:
      proc_sim <- mod %>%
        pomp::rprocess(
          .,
          x0 = samp_end_state,
          t0 = max(mod@times),
          times = new_times,
          params = mod@params
        )

      measures <- mod %>%
        rmeasure(
          object = .,
          x = proc_sim,
          times = new_times,
          params = mod@params
        ) %>%
        drop() %>% t() %>% as.data.frame()
    }


    # One of the dimensions in the array is empty, so we drop it. After, the
    # rows are the states, and columns are the times, but it is more natural
    # to have this transposed.
    sim_df <- as.data.frame(t(drop(proc_sim)))
    sim_df$time <- new_times
    sim_df$.id  <- i
    results <- dplyr::bind_cols(sim_df, measures)
    results %>%
      dplyr::select(time, .id, dplyr::everything())
  }
}
