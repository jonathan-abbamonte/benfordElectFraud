#' Extract a Digit from a Number
#'
#' @description Extract a digit (or digits) from a vector of numbers.
#'
#' @param x A vector of numbers to extract digits from.
#' @param digit An integer corresponding to the starting digit of the numbers in x to be extracted.
#' @param number.of.digits Number of digits to extract. Default is 1.
#' @param last If \code{TRUE}, function will extract digits from left to right, i.e. starting from the last digit.
#' If \code{FALSE} (default), function will extract digits from right to left, i.e. starting from the first digit.
#'
#' @return An integer vector of extracted digits.
#' @export
#'
#' @examples extract_digit(1:20, digit=1)
#' extract_digit(90:120, digit=2)
#' extract_digit(90:120, digit=1, number.of.digits=2)
#' extract_digit(1000:1030, number.of.digits=2, last=TRUE)
extract_digit <- function(x, digit, number.of.digits=1, last=FALSE){
  x <- abs(x)
  ifelse(number.of.digits > nchar(x),
         {stop('number.of.digits cannot be greater than the number of digits in x')},
         w <- 2)

  if (last==FALSE) {
    as.integer(as.numeric(substr(x,
                                 start = nchar(x)-(nchar(x)-digit),
                                 stop = nchar(x)-(nchar(x)-digit)+number.of.digits-1)))
  } else if (last == TRUE) {
    as.integer(as.numeric(substr(x,
                                 start = nchar(x)-number.of.digits+1,
                                 stop = nchar(x))))
  } else {stop("last must be either TRUE or FALSE")}
}

