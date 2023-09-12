#' Probability Mass Function of the Discrete Weibull Distribution
#'
#'
#' @param x A vector of values corresponding to realizations of the random variable.
#' @param a Scale parameter of the discrete weibull distribution. See details.
#' @param b Shape parameter of the discrete weibull distribution. See details.
#'
#'
#' @details
#' For this function, the parameterization for the discrete weibull distribution is:
#' \deqn{exp[-(\frac{x}{a})^b]-exp[-(\frac{x+1}{a})^b]}
#'
#'
#' @return A vector of values corresponding to the probability masses of the discrete weibull distribution.
#' @export
#'
#' @examples d.weibull(1:10, a=10, b=1)
d.weibull <- function(x, a, b) {

  val = exp(-(x/a)^b) - exp(-((x+1)/a)^b)

  return(val)
}




