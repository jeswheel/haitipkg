##' Initial parameter values for Model 1: haiti1()
##'
##' Dataset containing the 2 columns: epidemic MIF2 results and endemic MIF2 results for Model 1
##'
##' @format A data frame with 2 rows
##' \describe{
##'     \item{rho}{reporting rate}
##'     \item{tau}{dispersion in observation process}
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
##'     \item{sig_sq}{dispersion in latent process}
##'     \item{S_0}{initial proportion of population Susceptible at start of time window}
##'     \item{E_0}{initial proportion of population Exposed at start of time window}
##'     \item{I_0}{initial proportion of population Infectious at start of time window}
##'     \item{A_0}{initial proportion of population Asymptomatic at start of time window}
##'     \item{R_0}{initial proportion of population Recovered at start of time window}
##'     \item{pop_0}{initial population at start of time window}
##' }
##'
"haiti1_adj_params"
