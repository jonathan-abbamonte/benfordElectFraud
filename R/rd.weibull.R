#' Random Samples from Discrete Weibull Distribution
#'
#' @param n A scalar correpsonding to the number of draws from discrete weibull mass function.
#' @param a Scale parameter of the discrete weibull distribution. See details.
#' @param b Shape parameter of the discrete weibull distribution. See details.
#'
#'#' @details
#' For this function, the parameterization for the discrete weibull distribution is:
#' \deqn{exp[-(\frac{x}{a})^b]-exp[-(\frac{x+1}{a})^b]}
#'
#' @return A vector of values corresponding to quantiles from the discrete weibull distribution.
#' @export
#'
#' @examples rd.weibull(10, a=5, b=1)
rd.weibull <- function(n, a, b){
 val_vec = runif(n, 0, 1)
 vec = a*(-log(1-val_vec))^(1/b)-1
 return(vec)
}
