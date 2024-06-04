#' Distribution of Probabilities Predicted by Benford's Law
#'
#' @description
#' Compute the distribution expected probabilities for any-digit Benford's Law for any base.
#'
#'
#' @param digit digit from significand to use for calculating distribution of probabilities. Default is \code{1}, and any positive integer is accepted.
#' @param base the number base to use. Default is \code{10}.
#'
#' @return a vector of probabilities corresponding to each digit.
#' @export
#'
#' @references Anderson K.M., Dayaratna K., Gonshorowski D., & Miller S.J. (2022). "A New Benford Test for Clustered Data with Applications to American Elections." \emph{Stats}, 5(3), 841-855. https://doi.org/10.3390/stats5030049
#'
#' @examples benford_distribution(base=10)
benford_distribution = function(digit=1, base=10) {
    if (digit==1){
      return(setNames(sapply(seq(1,base-1), function(x) {log((1+1/x), base=base)}), seq(1,base-1)))
    } else {
      return (setNames(sapply(seq(0,base-1), function(x){sum(sapply(seq(base^(digit-2),base^(digit-1)-1), function(k) {log(1+(1/(base*k+x)), base=base)}))}), seq(0,base-1)))
    }
  }
