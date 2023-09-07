#' Extract a Digit from a Number
#'
#' @description Extract a digit (or digits) from a vector of numbers.
#'
#' @param x a vector of numbers to extract digits from.
#' @param digit an integer corresponding to the digit of the numbers in x to be extracted.
#' @param number.of.digits number of digits to extract. Default is 1.
#'
#' @return an integer vector of extracted digits.
#' @export
#'
#' @examples extract_digit(1:20, digit=1)
#' extract_digit(90:120, digit=2)
#' extract_digit(90:120, digit=1, number.of.digits=2)
extract_digit <- function(x, digit, number.of.digits=1){
  x <- abs(x)
  ifelse(number.of.digits > nchar(x),
         {stop('number.of.digits cannot be greater than the number of characters in x')},
         w <- 2)

  as.integer(as.numeric(substr(x, nchar(x)-(nchar(x)-digit), nchar(x)-(nchar(x)-digit)+number.of.digits-1)))
}

