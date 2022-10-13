#' Get starting values for fitting Model 1
#'
#' This function get a set of starting parameters that are intended to
#' reproduce the results of Lee et al (2020). The starting parameter values
#' were not provided, so this function attempts to estimate them based on the
#' description of their distributions provided in the supplement material.
#'
#' @param n_params number of starting parameters to use.
#' @return data.frame containing
#'
#' @import dplyr
#' @export
get_h1_startVals <- function(n_params) {
  set.seed(141837)

  # Sample beta and nu
  starting_params <- sample_beta_nu(n_params = n_params)

  # Sample rho
  starting_params$rho <- sample_rho(
    beta = starting_params$beta, n_params = n_params
  )

  # Sample tau
  starting_params$tau <- sample_tau(n_params = n_params)

  # Sample IVPs
  POP <- 10911819
  E_0_probs <- dnorm(seq(1, 1752, 1), mean = 400, sd = (1752-400) / qnorm(0.99))
  I_0_probs <- dnorm(seq(1, 2271, 1), mean = 500, sd = (2271-500) / qnorm(0.99))
  starting_params$E_0 <- sample(1:1752, size = n_params, replace = TRUE, prob = E_0_probs) / POP
  starting_params$I_0 <- sample(1:2271, size = n_params, replace = TRUE, prob = I_0_probs) / POP
  starting_params$S_0 <- 1 - starting_params$E_0 - starting_params$I_0

  starting_params
}

#' Sample starting values for beta and nu
#'
#' This helper function samples starting values for beta and nu, using their
#' joint density as shown in the supplement material of Lee et al (2020) as
#' a guide.
#'
#' @param mu_beta mean of the marginal distribution of beta. Value based on
#'    visual inspection of figure in Lee et al (2020).
#' @param mu_nu mean of the marginal distribution of nu.
#'    Value based on visual inspection of figure in Lee et al (2020)
#' @param n_params number of starting parameters to use.
#' @noRd
sample_beta_nu <- function(mu_beta = 5.8, mu_nu = 0.96, n_params) {
  # Determine a good value for the marginal standard deviations
  sd_beta <- (12.1-mu_beta) / qnorm(0.995)
  sd_nu <- (0.9-mu_nu) / qnorm(0.005)

  # Create covariance matrix for the bivariate-normal distribution
  cov_mat <- rbind(
    c(sd_beta^2, -0.7 * sd_beta * sd_nu),
    c(-0.7 * sd_beta * sd_nu, sd_nu^2)
  )

  # Sample from the bivariate normal distribution
  X <- MASS::mvrnorm(n = n_params, mu = c(mu_beta, mu_nu), Sigma = cov_mat)
  colnames(X) <- c("beta1", "nu")

  # Identify values that fell outside of range, so they can be resampled
  bad_beta <- which(X[, 'beta1'] > 12.1 | X[, 'beta1'] < 0.04)
  bad_nu <- which(X[, 'nu'] > 1 | X[, 'nu'] < 0.9)
  n_resamp <- length(unique(c(bad_beta, bad_nu)))

  # Keep resampling values until all values fall within the range.
  while (n_resamp > 0) {
    to_resamp <- unique(c(bad_beta, bad_nu))
    X[to_resamp, ] <- MASS::mvrnorm(n = n_resamp, mu = c(mu_beta, mu_nu), Sigma = cov_mat)

    bad_beta <- which(X[, 'beta1'] > 12.1 | X[, 'beta1'] < 0.04)
    bad_nu <- which(X[, 'nu'] > 1 | X[, 'nu'] < 0.9)
    n_resamp <- length(unique(c(bad_beta, bad_nu)))
  }

  # Return results as a data.frame
  as.data.frame(X)
}

#' Sample starting values for tau
#'
#' Function used to sample tau. Tau is approximately independent of the other
#' parameters.
#'
#' @param mu mean of the truncated normal distribution
#' @param n_params number of starting parameters to use.
#' @return vector of length n_params, sampled values of tau from truncated
#'    normal distribution
#' @noRd
sample_tau <- function(mu = 4.1, n_params) {
  X <- rnorm(n_params, mean = mu, sd = (1.7-mu) / qnorm(0.001))

  # Get how many resamples are needed
  n_resamp <- sum(X > 5.9 | X < 1.7)

  # Keep resampling until all values fall within desired range
  while (n_resamp > 0) {
    X[X > 5.9 | X < 1.7] <- rnorm(n_resamp, mean = mu, sd = (1.7-mu) / qnorm(0.001))

    n_resamp <- sum(X > 5.9 | X < 1.7)
  }

  # Return values of tau
  X
}

#' Sample starting values for rho
#'
#' This function is used to get starting values for the parameter rho,
#'  based on it's strong correlation with beta1 and nu.
#'  Here, we first sample beta1 and nu, and then sample rho based on a
#'  regression on beta1. Additionally, it is noted that rho has a bi-modal
#'  distribution, with lower values of rho corresponding with high
#'  probability to higher values of beta1.
#'
#' @param beta previously sampled values of beta1 that will be used to get
#'    values of rho.
#' @param n_params number of starting parameters to use.
#' @return vector of length n_params that will have an appropriate relationship
#'    with the parameter beta1.
#' @noRd
sample_rho <- function(beta, n_params) {
  rhos <- rep(0, n_params)

  # Get slope and intercept for the regression
  slope <- (0.4/(5.1 - 6.75))
  intercpt <- -6.75 * slope

  # Determine which mode should the value of rho belong to, based on regression
  case_when(
    beta > 6.75 ~ 0,  # Large beta --> group "0"
    # Small beta --> more likely to be group "1", with increasing probability as beta gets smaller
    beta <= 6.75 & beta > (0.95 - intercpt) / slope ~ slope * beta + intercpt,
    TRUE ~ 0.95
  ) -> probs

  # Sample the groups that rho belongs to
  rho_group1 <- rbernoulli(n_params, probs)
  rho_group2 <- !rho_group1

  # Keep track of how many values are in each group, used later for resampling
  n_group1 <- sum(rho_group1)
  n_group2 <- sum(rho_group2)

  # Resample values of rho that are larger than 1
  group1 <- rnorm(n_group1, mean = 1, sd = 0.05)
  n_bad1 <- sum(group1 > 1)
  while(n_bad1 > 0) {
    group1[group1 > 1] <- rnorm(n_bad1, mean = 1, sd = 0.175)
    n_bad1 <- sum(group1 > 1)
  }

  # get standard deviation for regression
  std_beta1 <- (beta[rho_group2] - mean(beta[rho_group2])) / sd(beta[rho_group2])

  group2 <- 0.25 - 0.025 * std_beta1 - (rlnorm(n_group2, sdlog = 0.2) - 1)

  # Resample values of rho in group2 as needed
  which_bad2 <- group2 < .13
  while(sum(which_bad2) > 0) {
    group2[group2 < 0.13] <- 0.25 - 0.1 * std_beta1[which_bad2] - (rlnorm(sum(which_bad2), sdlog = 0.2) - 1)

    which_bad2 <- group2 < .13
  }

  rhos[rho_group1] <- group1
  rhos[rho_group2] <- group2

  # Return the values of rho
  rhos
}


