#' Autoregressive Negative Binomial
#'
#' Fits a simple AR1 negative binomial model. Let \eqn{Y_1, \ldots, Y_N} be
#' the \eqn{N} random variables representing observed data
#' \eqn{y_1^*, \ldots y_N^*}, respectively. The model assumes that
#' \eqn{Y_n|Y_{n-1}} for some \eqn{n \in 1:N} has the following distribution:
#' \deqn{Y_n | Y_{n - 1} = y_{n-1}* \sim NB(a + by_{n-1}^*, \psi)}
#' Where \eqn{a, b} and \eqn{\psi} are model parameters that need to be
#' estimated.
#'
#' Missing values are currently ignored.
#'
#' @param data The data to fit, numeric vector.
#' @param init Initial value
#' @param start_vals Starting values for the parameters
#' @param transform Boolean, whether or not to use variable transformations in the optimization
#' @param a Only used if `transform == TRUE`. `a` is a numeric representing the maximum value
#'    for the AR1 coefficient, and `-a` is the minimum value.
#' @param ... other parameters to be passed into the `optim` function
#'
#' @return An AR1_NegBinom model
#' @export
#'
#' @examples ar1_NegBinom(round(sunspot.year), init = 5)
ar1_NegBinom <- function(data, init, start_vals = c(10, 0.5, 1), transform = TRUE, a = 2, ...) {

  # Fit the model
  #   - Done by creating a new function that calls conditional ll, and takes the sum.

  neg_ll <- function(pars) {
    -sum(.ar1_NegBinom_cond_ll(theta = pars, init = init, data = data, transform = transform))
  }

  if (transform) start_vals <- .undoTrans(start_vals, a = 2)

  out <- optim(start_vals, neg_ll, ...)


  res <- list()

  if (transform) {
    res$theta <- .doTrans(out$par, a = 2)
  } else {
    res$theta <- out$par
  }

  names(res$theta) <- c("Intercept", "AR1", "Size")

  res$cond_ll <- .ar1_NegBinom_cond_ll(theta = res$theta, init = init, data = data, transform = FALSE)
  res$ll <- -out$value

  class(res) <- 'AR1_NegBinom'
  res
}

#' Conditional log likelihood of an AR1 Negative Binomial model
#'
#' The log likelihood of the model is equal to the sum of the conditional log
#' likelihoods; because we are also interested at times in the conditional log
#' likelihoods, the sum of this function is what is optimized by `optim`, and
#' this function is later called when the `AR1_NegBinom` object is created.
#'
#' @param theta parameter values: (a, b, psi)
#' @param init initial value Y0
#' @param data data used to calculate the conditional log likelihood
#' @param transform whether or not to transform the parameters
#'
#' @return vector of conditional log-likelihoods.
#' @noRd
.ar1_NegBinom_cond_ll <- function(theta, init, data, transform) {

  lls <- numeric(length(data))

  if (transform) theta <- .doTrans(theta)

  a <- theta[1]
  b <- theta[2]
  psi <- theta[3]
  # xinit <- as.integer(round(theta[4]))

  x0 <- c(init, data)
  for (i in 1:(length(data))) {

    if (!is.na(x0[i + 1])) { # Obs not missing

      if (is.na(x0[i])) { # Previous obs is missing, go back two obs.
        mu <- a + b * x0[i - 1]
        lls[i] <- dnbinom(x0[i + 1], mu = mu, size = psi, log = TRUE)
      } else { # Previous obs isn't missing, and current obs ins't missing.
        mu <- a + b * x0[i]
        lls[i] <- dnbinom(x0[i + 1], mu = mu, size = psi, log = TRUE)
      }
    }

    # cat(paste0('i: ', i, "; logLik = ", log_lik, "\n"))
  }

  lls
}

#' Transform all reals to (-a, a)
#'
#' @param x real number to be transformed
#'
#' @return value in (-a, a)
#' @noRd
.trans_neg_a_a <- function(x, a = 2) {
  if (x < 0) {
    (a * exp(x) - a) / (exp(x) + 1)
  } else {
    (a - a * exp(-x)) / (1 + exp(-x))
  }
}

#' Transform from (-a, a) to the real number line
#'
#' @param x real number in (-a, a) to transform
#'
#' @return real number
#' @noRd
.undoTrans_neg_a_a <- function(x, tol = 1e-6, a = 2) {
  log((a + x) / (a - (x - tol)))
}

#' Performs transformation from R^3 -> R+ x (-a, a) x R+
#'
#' @param theta vector of three real numbers
#'
#' @return vector of real numbers. The first and last elements will be positive,
#'    the second element will be in (-a, a).
#' @noRd
.doTrans <- function(theta, a = 2) {
  c(exp(theta[1]), .trans_neg_a_a(theta[2], a = a), exp(theta[3]))
}

#' Performs transformation from R+ x (-a, a) x R+ -> R^3
#'
#' @param theta vector of three real numbers
#'
#' @return vector of real numbers.
#' @noRd
.undoTrans <- function(theta, a = 2) {
  c(log(theta[1]), .undoTrans_neg_a_a(theta[2], a = a), log(theta[3]))
}

#' @method print AR1_NegBinom
#' @export
print.AR1_NegBinom <- function(x, ...) {
  cat("Coefficients: \n\n")
  print.default(round(x$theta, 4), print.gap = 2)

  cat(paste0("\nLog likelihood = ", round(x$ll, 4)))
}
