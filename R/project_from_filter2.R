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
#'
#' @importFrom magrittr %>%
#' @importFrom foreach %dopar%
#' @import progress
#' @export

project_from_filter2 <- function(mod, PF, covarGen = NULL,
                                 nsims = 100, seed = 123) {

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

  departments <- c(
    'Artibonite', 'Centre', 'Grande_Anse', 'Nippes', 'Nord',
    'Nord-Est', 'Nord-Ouest', 'Ouest', 'Sud', 'Sud-Est'
  )

  # Check if there are saved states, otherwise everything below fails
  if (class(PF) == "bpfilterd_spatPomp") {
    if (length(PF@saved.states) == 0) stop("bpfilterd_spatPomp object must have saved.states.")
  } else if (length(saved.states(PF)) == 0) {
    stop("pfilter object must have saved.states.")
  }

  # Get the saved states from filtering distribution.
  if (class(PF) == 'bpfilterd_spatPomp') {
    ss <- PF@saved.states
  } else {
    ss <- pomp::saved.states(PF)  # This function doesn't yet exist for bpfiltered_spatPomp
  }

  end_states <- ss[[length(ss)]]

  rm(ss, PF)
  gc()

  state_names <- rownames(end_states)

  # Save a vector of the times we want to simulate in the future. We
  # start one week after the last day of available data, because the
  # last day of available data becomes t0.
  if (mod@timename == "week") {
    new_times <- (max(mod@times) + 1):(max(mod@times) + 570)
  } else {
    new_times <- dateToYears(
      seq.Date(
        yearsToDate(max(mod@times)) + 7,
        as.Date("2029-12-20"),
        by = "1 week"
      )
    )
  }


  if (class(mod) == 'spatPomp') {


    N_future_covar <- nrow(covarGen(include_data = FALSE))

    # Get the observed rainfall first
    covar_data <- haitiRainfall %>%
      dplyr::filter(date >= as.Date("2010-10-23") - 7) %>%
      dplyr::summarize(
        date = date, dplyr::across(Artibonite:`Sud-Est`, std_rain)
      ) %>%
      dplyr::mutate(
        time = dateToYears(date)
      )

    colnames(covar_data) <- c(
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

    covar_data <- covar_data %>% dplyr::select(time, dplyr::starts_with("rain_std")) %>%
      as.matrix()

    pb_covar <- progress_bar$new(format = "Creating Covariates: [:bar] :percent [:elapsedfull]",
                                 total = nsims - 1,
                                 complete = "=",   # Completion bar character
                                 incomplete = "-", # Incomplete bar character
                                 current = ">",    # Current bar character
                                 clear = FALSE,    # If TRUE, clears the bar when finish
                                 width = 75)      # Width of the progress bar

    covar_matrix <- covarGen(include_data = FALSE)
    for (i in 1:(nsims - 1)) {
      pb_covar$tick()
      covar_matrix <- rbind(covar_matrix, covarGen(include_data = FALSE))
    }

    pb_sims <- progress_bar$new(format = "Simulating Model   : [:bar] :percent [:elapsedfull]",
                                total = nsims,
                                complete = "=",   # Completion bar character
                                incomplete = "-", # Incomplete bar character
                                current = ">",    # Current bar character
                                clear = FALSE,    # If TRUE, clears the bar when finish
                                width = 75)      # Width of the progress bar

    foreach::foreach(
      i = 1:nsims,
      .combine = dplyr::bind_rows
    ) %do% {

      pb_sims$tick()

      future_covar <- covar_matrix[((i-1)*N_future_covar + 1):(i*N_future_covar), ]
      covar_df <- as.data.frame(rbind(covar_data, future_covar)) %>%
        tidyr::pivot_longer(
          data = .,
          cols = 2:11,
          names_to = "departement",
          values_to = "rain_std",
          names_prefix = "rain_std"
        ) %>%
        dplyr::arrange(time, departement) %>%
        dplyr::mutate(
          departement = gsub('-', '_', departement)
        )

      proc_sim <- mod %>%
        spatPomp::spatPomp(
          covar = as.data.frame(covar_df)
        ) %>%
        pomp::rprocess(
          .,
          x0 = end_states[, sample(ncol(end_states), size = 1)],
          t0 = max(mod@times),
          times = new_times,
          params = mod@params
        )

      measures <- mod %>%
        spatPomp::spatPomp(
          covar = as.data.frame(covar_df)
        ) %>%
        rmeasure(
          object = .,
          x = proc_sim,
          times = new_times,
          params = mod@params
        ) %>%
        drop() %>% t() %>% as.data.frame()


      # One of the dimensions in the array is empty, so we drop it. After, the
      # rows are the states, and columns are the times, but it is more natural
      # to have this transposed.
      sim_df <- as.data.frame(t(drop(proc_sim)))
      sim_df$time <- new_times
      sim_df$.id  <- i
      results <- dplyr::bind_cols(sim_df, measures)

      rm(
        proc_sim, measures, sim_df
      )
      gc()

      get_unit_index_from_statename <- function(statename){
        stringr::str_extract(statename, "[[:digit:]]+$")
      }
      get_state_obs_type_from_statename <- function(statename) {
        gsub("[[:digit:]]+$", "", statename)
      }

      results %>%
        dplyr::select(time, .id, dplyr::everything()) %>%
        tidyr::pivot_longer(
          cols = -c(time, .id),
          names_to = 'stateobs',
          values_to = 'val'
        ) %>%
        dplyr::mutate(
          ui = get_unit_index_from_statename(stateobs),
          unit = unit_names(mod)[as.integer(ui)]
        ) %>%
        dplyr::select(
          time, .id, unit, stateobs, val
        ) %>%
        dplyr::arrange(.id, time, unit, stateobs) %>%
        dplyr::mutate(stateobstype = get_state_obs_type_from_statename(stateobs)) %>%
        dplyr::select(-stateobs) %>%
        tidyr::pivot_wider(names_from = stateobstype, values_from = 'val') %>%
        dplyr::rename(unitname = unit)
    }
  } else if (!is.null(covarGen)) { ################################### POMP ##################

    N_future_covar <- nrow(covarGen(include_data = FALSE))

    # Get the observed rainfall first
    covar_data <- haitiRainfall %>%
      dplyr::filter(date >= as.Date("2010-10-23") - 7) %>%
      dplyr::summarize(
        date = date, dplyr::across(Artibonite:`Sud-Est`, std_rain)
      ) %>%
      dplyr::mutate(
        time = dateToYears(date)
      )

    colnames(covar_data) <- c(
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

    covar_data <- covar_data %>% dplyr::select(time, dplyr::starts_with("rain_std")) %>%
      as.matrix()

    pb_covar <- progress_bar$new(format = "Creating Covariates: [:bar] :percent [:elapsedfull]",
                                 total = nsims - 1,
                                 complete = "=",   # Completion bar character
                                 incomplete = "-", # Incomplete bar character
                                 current = ">",    # Current bar character
                                 clear = FALSE,    # If TRUE, clears the bar when finish
                                 width = 75)      # Width of the progress bar

    covar_matrix <- covarGen(include_data = FALSE)
    for (i in 1:(nsims - 1)) {
      pb_covar$tick()
      covar_matrix <- rbind(covar_matrix, covarGen(include_data = FALSE))
    }

    pb_sims <- progress_bar$new(format = "Simulating Model   : [:bar] :percent [:elapsedfull]",
                                total = nsims,
                                complete = "=",   # Completion bar character
                                incomplete = "-", # Incomplete bar character
                                current = ">",    # Current bar character
                                clear = FALSE,    # If TRUE, clears the bar when finish
                                width = 75)      # Width of the progress bar

    foreach::foreach(
      i = 1:nsims,
      .combine = dplyr::bind_rows
    ) %do% {

      pb_sims$tick()

      future_covar <- covar_matrix[((i-1)*N_future_covar + 1):(i*N_future_covar), ]

      proc_sim <- mod %>%
        pomp::pomp(
          covar = pomp::covariate_table(rbind(covar_data, future_covar), times = 'time')
        ) %>%
        pomp::rprocess(
          .,
          x0 = end_states[, sample(ncol(end_states), size = 1)],
          t0 = max(mod@times),
          times = new_times,
          params = mod@params
        )

      measures <- mod %>%
        pomp::pomp(
          covar = pomp::covariate_table(rbind(covar_data, future_covar), times = 'time')
        ) %>%
        rmeasure(
          object = .,
          x = proc_sim,
          times = new_times,
          params = mod@params
        ) %>%
        drop() %>% t() %>% as.data.frame()


      # One of the dimensions in the array is empty, so we drop it. After, the
      # rows are the states, and columns are the times, but it is more natural
      # to have this transposed.
      sim_df <- as.data.frame(t(drop(proc_sim)))
      sim_df$time <- new_times
      sim_df$.id  <- i
      results <- dplyr::bind_cols(sim_df, measures)

      rm(
        proc_sim, measures, sim_df
      )
      gc()

      results %>%
        dplyr::select(time, .id, dplyr::everything())
    }
  } else {

    sample_numbers <- sample(ncol(end_states), size = nsims, replace = TRUE)
    x_filter <- end_states[, sample_numbers]

    proc_sim <- mod %>%
      pomp::rprocess(
        .,
        x0 = x_filter,
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
      )

    if (nsims != dim(proc_sim)[2]) {
      stop("Wrong Dimensions in nsims in proc_sim.")
    } else if (length(new_times) != dim(proc_sim)[3]) {
      stop('Wrong Dimensions in times in proc_sim')
    }

    ntimes <- length(new_times)
    simnames <- rownames(proc_sim)

    final_out <- matrix(nrow = nsims * ntimes, ncol = length(state_names))
    measures_temp <- matrix(nrow = nsims * ntimes, ncol = dim(measures)[1])

    for (i in 1:nsims) {
      final_out[((i - 1) * ntimes + 1):((i) * ntimes), ] <- t(proc_sim[, i, ])
      measures_temp[((i - 1) * ntimes + 1):((i) * ntimes), ] <- measures[, i, ]
    }

    final_out <- cbind(rep(new_times, nsims), final_out, measures_temp, rep(1:nsims, each = ntimes))
    colnames(final_out) <- c(mod@timename, state_names, rownames(mod@data), ".id")


    # One of the dimensions in the array is empty, so we drop it. After, the
    # rows are the states, and columns are the times, but it is more natural
    # to have this transposed.
    results <- as.data.frame(final_out)

    rm(
      proc_sim, measures
    )
    gc()

    results %>%
      dplyr::select(mod@timename, .id, dplyr::everything())
  }

}
