#' Set Model 3 hyper-parameters for model fitting
#'
#' This function sets hyper-parameters needed to fit Model 3. The primary
#' argument is RUN_LEVEL, which determines the computational effort used to
#' fit the model. The computation effort is primarily determined by the number
#' of particles, number of BPF iterations, and the number of replicates.
#'
#' Originally, this function did not exist and these hyper-parameters were
#' just set in the R code right before fitting the model. Ultimately this
#' function was created with the intention of making the code in the ms.Rnw
#' file that uses this package easier to read.
#'
#' @param RUN_LEVEL Integer value 1-3. This determines the computational effort
#'   used to fit the model. \code{RUN_LEVEL == 1} is used primarily for debugging
#'   purposes, \code{RUN_LEVEL == 2} uses sufficient computational effort to
#'   obtain preliminary results, and \code{RUN_LEVEL == 3} gives results of the
#'   full computation.
#' @param rw_1.sd Random walk standard deviation specification for the first
#'   search.
#' @param rw_2.sd Random walk standard deviation specification for the second
#'   search.
#' @param rw_3.sd Random walk standard deviation specification for the third
#'   search.
#'
#' @return A list containing hyper-parameters for fitting model 3.
#' @export
#'
#' @examples
#' rwsd <- rw_sd(epsilon = 0.02)
#' set_h3_fit_parms(RUN_LEVEL = 3, rw_1.sd = rwsd, rw_2.sd = rwsd, rw_3.sd = rwsd)
set_h3_fit_parms <- function(RUN_LEVEL, rw_1.sd, rw_2.sd, rw_3.sd) {


  if (RUN_LEVEL == 1){
    SEARCH1 = list(
      NBPF = 15,
      NP = 500,
      SPAT_REGRESSION = 0.05,
      NREPS = 36,
      NP_EVAL = 500,
      NREPS_EVAL = 3,
      RW_SD = rw_1.sd,
      COOLING = 0.5,
      KEEP_TRACES = TRUE
    )

    SEARCH2 = list(
      TOP_N = 9,
      NBPF = 10,
      NP = 500,
      SPAT_REGRESSION = 0.05,
      NREPS = 4,
      NP_EVAL = 500,
      NREPS_EVAL = 3,
      RW_SD = rw_2.sd,
      COOLING = 0.5,
      KEEP_TRACES = TRUE
    )

    SEARCH3 = list(
      TOP_N = 6,
      NBPF = 10,
      NP = 500,
      SPAT_REGRESSION = 0.05,
      NREPS = 6,
      NP_EVAL = 500,
      NREPS_EVAL = 3,
      RW_SD = rw_3.sd,
      COOLING = 0.5,
      KEEP_TRACES = TRUE
    )

    N_SEARCHES <- 3L
  } else if (RUN_LEVEL == 2) {
    SEARCH1 = list(
      NBPF = 100,
      NP = 1000,
      SPAT_REGRESSION = 0.05,
      NREPS = 72,
      NP_EVAL = 1000,
      NREPS_EVAL = 8,
      RW_SD = rw_1.sd,
      COOLING = 0.5,
      KEEP_TRACES = TRUE,
      KEEP_LIK_MAT = TRUE
    )

    SEARCH2 = list(
      TOP_N = 9,
      NBPF = 100,
      NP = 1000,
      SPAT_REGRESSION = 0.05,
      NREPS = 4,
      NP_EVAL = 1000,
      NREPS_EVAL = 9,
      RW_SD = rw_2.sd,
      COOLING = 0.5,
      KEEP_TRACES = TRUE,
      KEEP_LIK_MAT = TRUE
    )

    SEARCH3 = list(
      TOP_N = 6,
      NBPF = 100,
      NP = 1000,
      SPAT_REGRESSION = 0.05,
      NREPS = 6,
      NP_EVAL = 1000,
      NREPS_EVAL = 9,
      RW_SD = rw_3.sd,
      COOLING = 0.5,
      KEEP_TRACES = TRUE,
      KEEP_LIK_MAT = TRUE
    )

    N_SEARCHES <- 3L
  } else if (RUN_LEVEL == 3) {

    SEARCH1 = list(
      NBPF = 100,
      NP = 2000,
      SPAT_REGRESSION = 0.05,
      NREPS = 144,
      NP_EVAL = 2000,
      NREPS_EVAL = 18,
      RW_SD = rw_1.sd,
      COOLING = 0.5,
      KEEP_TRACES = FALSE,
      KEEP_LIKE_MAT = FALSE
    )

    SEARCH2 = list(
      TOP_N = 12,
      NBPF = 100,
      NP = 2000,
      SPAT_REGRESSION = 0.05,
      NREPS = 6,
      NP_EVAL = 2000,
      NREPS_EVAL = 18,
      RW_SD = rw_2.sd,
      COOLING = 0.5,
      KEEP_TRACES = FALSE,
      KEEP_LIK_MAT = FALSE
    )

    SEARCH3 = list(
      TOP_N = 9,
      NBPF = 100,
      NP = 2000,
      SPAT_REGRESSION = 0.05,
      NREPS = 8,
      NP_EVAL = 2000,
      NREPS_EVAL = 18,
      RW_SD = rw_3.sd,
      COOLING = 0.5,
      KEEP_TRACES = FALSE,
      KEEP_LIK_MAT = TRUE
    )

    N_SEARCHES <- 3L
  }

  out <- list()
  out$SEARCH1 <- SEARCH1
  out$SEARCH2 <- SEARCH2
  out$SEARCH3 <- SEARCH3
  out$N_SEARCHES <- N_SEARCHES

  out
}
