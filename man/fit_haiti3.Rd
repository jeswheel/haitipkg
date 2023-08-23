% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_haiti3.R
\name{fit_haiti3}
\alias{fit_haiti3}
\title{Fit Model 3 SpatPomp}
\usage{
fit_haiti3(
  search1 = list(NBPF = 5, NP = 50, SPAT_REGRESSION = 0.5, NREPS = 6, NP_EVAL = 100,
    NREPS_EVAL = 6, RW_SD = NULL, COOLING = 0.5, KEEP_TRACES = FALSE, KEEP_LIK_MAT =
    FALSE),
  search2 = NULL,
  search3 = NULL,
  ncores = 3,
  nsearches = 1,
  search_rho = FALSE,
  search_gamma = FALSE,
  search_hur = FALSE,
  search_Iinit = FALSE,
  start_date = "2010-11-20"
)
}
\arguments{
\item{search1}{list containing parameters used to fit the model. See details.}

\item{search2}{list containing parameters used to fit the model. See details.}

\item{search3}{list containing parameters used to fit the model. See details.}

\item{ncores}{Number of cores used to fit the model. The code is written
so that the optimal number of cores with \code{RUN_LEVEL = 3} is 36.}

\item{nsearches}{integer number of searches to conduct. See details below.}

\item{search_rho}{Boolean indicating whether or not rho should be estimated.}

\item{search_gamma}{Boolean indicating whether or not gamma should be estimated.}

\item{search_hur}{Boolean indicating whether or not coefficients related to
Hurricane Matthew (2016) should be perturbed in the global search.}

\item{search_Iinit}{Boolean indicating whether or not units with zero case
counts at time t1 should have the value \eqn{I_u(t0)} estimated.}

\item{start_date}{The starting date for the model. All of the data starts
at "2010-10-23", but the data is suspect until "2010-11-20". See supplement
for more detail.}
}
\description{
This function fits the spatPomp version of Model 3 to cholera incidence
data from Oct 2010 - Jan 2019 in Haiti. The data are measured at a
weekly timescale. The function allows a user to fit the model in up to 3
stages by setting \emph{nsearches} equal to 3. Any number larger than 3
is ignored.
}
\details{
Hyper-parameters to the IBPF algorithm can be adjusted as inputs into a list
for each search (e.g. search1 = list(TOP_N = 3, ...)). The possible list
options are:
\itemize{
\item{TOP_N: }{Number of results from the previous search to use as a
starting point for the next search. Note that this argument is ignored for
the first global search, and for searches after local unit searches.}
\item{NBPF: }{Number of IBPF iterations.}
\item{NP: }{Number of particles to use in in the IBPF algorithm.}
\item{SPAT_REGRESSION: }{The regression coefficient for shared parameters}
\item{NREPS: }{For each starting point, how many replicated searches should
be conducted. Here, we recommend that NREPS * TOP_N be a multiple of
ncores.}
\item{NP_EVAL: }{Number of particles to use in model evaluation.}
\item{NREPS_EVAL: }{Number of times to replicate the evaluation using a
particle filter.}
\item{RW_SD: }{Specification of the rw.sd used in the IBPF algorithm.}
\item{COOLING: }{Optional, defaults to 0.5. Specifies the
cooling.fraction.50 for the rw.sd used in the IBPF algorithm.}
\item{KEEP_TRACES: }{Optional, defaults to FALSE. Saves the traces of the
IBPF output. This is a large list, so is primarily saved for debugging and
testing purposes.}
\item{KEEP_LIK_MAT: }{Optional, defaults to FALSE. Saves a likelihood
matrix containing likelihood estimates for each unit. This is a large
matrix, so it is primarily saved for debugging and testing purposes.}
}

It is assumed that a parallel backend compatable with the foreach and
doRNG packages has been registered prior to calling this function
}
\examples{
\dontrun{fit_m3 <- fit_haiti3()}

}