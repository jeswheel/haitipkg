##' Initial parameter values for Model 1: haiti1_dep()
##'
##' Dataset containing the 2 columns: epidemic MIF2 results and endemic MIF2 results for Model 1, Departement specific
##'
##' @format A data frame with 26 rows and 10 columns (1 for each departement)
##' \describe{
##'     \item{rho_epi}{reporting rate for the epidemic period}
##'     \item{tau_epi}{dispersion in observation process for the epidemic period}
##'     \item{sig_sq_epi}{dispersion in latent process for the epidemic period}
##'     \item{rho_end}{reporting rate for the endemic period}
##'     \item{tau_end}{dispersion in observation process for the endemic period}
##'     \item{sig_sq_end}{dispersion in latent process for the endemic period}
##'     \item{beta1}{seasonality parameter for basis splines}
##'     \item{beta2}{seasonality parameter for basis splines}
##'     \item{beta3}{seasonality parameter for basis splines}
##'     \item{beta4}{seasonality parameter for basis splines}
##'     \item{beta5}{seasonality parameter for basis splines}
##'     \item{beta6}{seasonality parameter for basis splines}
##'     \item{nu}{mixing coefficient}
##'     \item{gamma}{infectious period (days)}
##'     \item{sigma}{latent period (days)}
##'     \item{theta0}{proportion of non-vac Exposed becoming Asymptomatic}
##'     \item{alpha}{mean years of natural immunity}
##'     \item{mu}{birth rate per 1000 per week}
##'     \item{delta}{natural death rate per 1000 per week}
##'     \item{S_0}{initial proportion of population Susceptible at start of time window}
##'     \item{E_0}{initial proportion of population Exposed at start of time window}
##'     \item{I_0}{initial proportion of population Infectious at start of time window}
##'     \item{A_0}{initial proportion of population Asymptomatic at start of time window}
##'     \item{R_0}{initial proportion of population Recovered at start of time window}
##'     \item{pop_0}{initial population at start of time window}
##' }
##'
"haiti1_dep_params"
