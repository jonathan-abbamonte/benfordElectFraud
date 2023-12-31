% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract_digit.R
\name{extract_digit}
\alias{extract_digit}
\title{Extract a Digit from a Number}
\usage{
extract_digit(x, digit, number.of.digits = 1)
}
\arguments{
\item{x}{a vector of numbers to extract digits from.}

\item{digit}{an integer corresponding to the digit of the numbers in x to be extracted.}

\item{number.of.digits}{number of digits to extract. Default is 1.}
}
\value{
an integer vector of extracted digits.
}
\description{
Extract a digit (or digits) from a vector of numbers.
}
\examples{
extract_digit(1:20, digit=1)
extract_digit(90:120, digit=2)
extract_digit(90:120, digit=1, number.of.digits=2)
}
