#' Make covariate table for Model 1
#'
#' Generate a data frame containing weekly values for periodic basis splines for the seasonal transmission term, \beta.
#'
#' @param tmin first observation time
#' @param tmax last observation time
#' @param byt number of steps between each time
#' @param nbasis number of basis spline terms
#' @param degree degrees of freedom for the basis splines
#' @param per period for the basis splines
#' @return A data frame with tmax rows and nbasis + 1 columns
#' @examples
#' covar <- covars(tmin = 0, tmax = 100, nbasis = 6, degree = 6, per = 52.14)
#' @export

covars <- function(tmin, tmax, byt = 1, nbasis = 6, degree = 6, per = 52.14) {
  tbasis <- seq(from = tmin, to = tmax, by = byt)
  covar_tab <- data.frame(cbind(time = tbasis,
                                pomp::periodic.bspline.basis(
                                  x = tbasis,
                                  nbasis = nbasis,
                                  degree = degree,
                                  period = per,
                                  names = "seas%d"
                                )))
  return(covar_tab)
}
