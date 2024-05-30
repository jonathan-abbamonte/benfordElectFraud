#' Random Samples from Discrete Weibull Distribution
#'
#' @param n A scalar correpsonding to the number of draws from discrete weibull mass function.
#' @param a Scale parameter of the discrete weibull distribution. See details.
#' @param b Shape parameter of the discrete weibull distribution. See details.
#'
#' @details
#' For this function, the parameterization for the discrete weibull distribution is:
#' \deqn{exp[-(\frac{x}{a})^b]-exp[-(\frac{x+1}{a})^b]}
#'
#' @return A vector of values corresponding to quantiles from the discrete weibull distribution.
#' @export
#'
#' @examples rd.weibull(10, a=5, b=1)
rd.weibull <- function(n, a, b){

  # x = 0
  # cum_prob <- 1 - exp(-((x+1)/a)^b)
  #
  # while (tail(cum_prob, n=1) < 1) {
  #   x = x + 1
  #   cum_prob_i = 1 - exp(-((x+1)/a)^b)
  #   cum_prob <- append(cum_prob, cum_prob_i)
  # }
  # probs = c(diff(cum_prob), 0)
  #
  # random_cum_probs <- sample(cum_prob, size=n, replace=TRUE, prob=probs)
  #
  # random_variates <- round( (a * ((-log(1-random_cum_probs))^(1/b))) - 1 )

  val_vec <- runif(n,0,1)
  vec = (a * ((-log(1-val_vec))^(1/b))) - 1
  random_variates <- round(vec)
  random_variates[random_variates < 0] <- 0
  return(random_variates)
}


