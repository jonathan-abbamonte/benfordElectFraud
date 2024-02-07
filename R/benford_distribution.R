#' Distribution of Probabilities Predicted by Benford's Law
#'
#' @description
#' Compute the distribution expected probabilities for first- or second-digit Benford's Law for any base.
#'
#'
#' @param digit digit from significand to use for calculating distribution of probabilities. Default is \code{"first"}
#'  with option to also select \code{"second"}.
#' @param base the number base to use.
#'
#' @return a vector of probabilities corresponding to each digit.
#' @export
#'
#' @references Anderson K.M., Dayaratna K., Gonshorowski D., & Miller S.J. (2022). "A New Benford Test for Clustered Data with Applications to American Elections." \emph{Stats}, 5(3), 841-855. https://doi.org/10.3390/stats5030049
#'
#' @examples benford_distribution(base=10)
benford_distribution <- function(digit="first", base) {

  bf <- c()

  if (digit == "first") {
    for (i in 1:(base-1)) {
      bf <- append(bf, log(1 + (1/i), base=base))
    }
    names(bf) <- seq(1, base-1, 1)
  }
  else if (digit == "second") {
    for (i in 0:(base-1)) {
      bfk <- c()
      for (k in 1:(base-1)) {
        bfk <- append(bfk, log(1 + (1/(base*k + i)), base=base))
      }
      bf[i+1] <- sum(bfk)
    }
    names(bf) <- seq(0, base-1, 1)
  }

  return (bf)
}
