% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/benford_test.R
\name{benford_test}
\alias{benford_test}
\title{Benford's Law Test}
\usage{
benford_test(
  x,
  digit = c("first", "second"),
  base,
  len = 100,
  group_id = NULL,
  id = NULL,
  na.rm = TRUE,
  conf.int = 0.95,
  method = c("chisq", "multinom"),
  p.adj = c("BH", "holm", "hochberg", "hommel", "bonferroni", "BY", "fdr", "none"),
  labs = NULL,
  verbose = 0
)
}
\arguments{
\item{x}{a data frame or an object coercible to a data frame such as a matrix containing values to test for Benfordness
and column corresponding to group names. Each column pertains to a response. Columns must contain only positive values.
Decimal values will be truncated to zero decimal places and coerced to integer.}

\item{digit}{a string corresponding to the digit of the values in x to test for Benfordness.
User can select either \code{"first"} or \code{"second"}.}

\item{base}{the target base to convert all non-negative integer values in x to. Must be 3 or greater. Required.}

\item{len}{the length of the character string used to convert numerals in base 10 to another base.
Function calls \code{oro.dicom:dec2base}. It is highly recommended to set the \code{len} > 100 as
there is a bug in \code{dec2base} which sets some values with leading 1's to '0'.
A higher \code{len} may be required for data with values of higher magnitude.}

\item{group_id}{a character variable to group Benford's test on. Groups must contain a minimum of 10 observations.
Any group with fewer than 10 observations will be excluded. Observations within groups will
be tested for Benfordness using the test specified in \code{method}. Required.}

\item{id}{a character vector specifying the id for observations.}

\item{na.rm}{how to handle missing values and zeros. If \code{na.rm=TRUE}, NAs and zeros will be
excluded from analysis. If \code{na.rm=FALSE} then the user must remove NAs \code{x} and observations
with zeros will be automatically dropped with a warning.}

\item{conf.int}{confidence level for the test selected in method. Default is \code{conf.int=0.90}.}

\item{method}{method used to calculate the test statistic and p-value for the Benford's Law test. User can
select "chisq" to compute using Pearson's Chi-Square test or "multinom" to compute using a multinomial likelihood
ratio test.}

\item{p.adj}{string specifying method for calculating adjusted p-value for multiple comparisons.
Function calls \code{stats::p.adjust}. Default method for is Benjamini & Hochberg (1995) ("BH") but any method
can be selected. See \code{?stats::p.adjust} for options.}

\item{labs}{named vector of labels to apply to each column in x. If no named vector is given, the names of the columns of x
(other than the column(s) assigned to group_id and id) are used as labels in the output.}
}
\value{
a data frame showing all grouped observations that are statistically different from the
expected counts predicted by Benford's Law at the confidence level set in \code{conf.int} and calculated using the
method specified in \code{method}. The group_id, label, test statistic, p-value, and group size are reported in the output.
}
\description{
Identify grouped observations that are statistically inconsistent with Benford's Law.
}
\examples{
benford_test(ohio2016, base=3, group_id="county_name", id="precinct", method="chisq")
}
\references{
Anderson K.M., Dayaratna K., Gonshorowski D., & Miller S.J. (2022). "A New Benford Test for Clustered Data with Applications to American Elections." \emph{Stats}, 5(3), 841-855. https://doi.org/10.3390/stats5030049
}
