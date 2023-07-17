#' Benford's Law Test
#'
#'@description
#'Identify grouped observations that are statically inconsistent with Benford's Law.
#'
#'
#' @param x a data frame or an object coercible to a data frame such as a matrix containing values to test for Benfordness
#' and column corresponding to group names. Each column pertains to a response. Columns must contain only positive values.
#' Decimal values will be truncated to zero decimal places and coerced to integer.
#'
#' @param digit a string corresponding to the digit of the values in x to test for Benfordness.
#' User can select either \code{"first"} or \code{"second"}.
#'
#' @param base the target base to convert all non-negative integer values in x to. Must be 3 or greater. Required.
#'
#' @param len the length of the character string used to convert numerals in base 10 to another base.
#' Function calls \code{oro.dicom:dec2base}. It is highly recommended to set the \code{len} > 100 as
#' there is a bug in \code{dec2base} which sets some values with leading 1's to '0'.
#' A higher \code{len} may be required for data with values of higher magnitude.
#'
#' @param group_id a character variable to group Benford's test on. Groups must contain a minimum of 10 observations.
#' Any group with fewer than 10 observations will be excluded. Observations within groups will
#' be tested for Benfordness using the test specified in \code{method}. Required.
#'
#' @param id a character vector specifying the id for observations.
#'
#' @param conf.int confidence level for the test selected in method. Default is \code{conf.int=0.90}.
#'
#' @param method method used to calculate the test statistic and p-value for the Benford's Law test. User can
#' select "chisq" to compute using Pearson's Chi-Square test or "multinom" to compute using a multinomial likelihood
#' ratio test.
#'
#' @param p.adj string specifying method for calculating adjusted p-value for multiple comparisons.
#' Function calls \code{stats::p.adjust}. Default method for is Benjamini & Hochberg (1995) ("BH") but any method
#' can be selected. See \code{?stats::p.adjust} for options.
#'
#' @param na.rm how to handle missing values and zeros. If \code{na.rm=TRUE}, NAs and zeros will be
#' excluded from analysis. If \code{na.rm=FALSE} then the user must remove NAs \code{x} and observations
#' with zeros will be automatically dropped with a warning.
#'
#' @param labs named vector of labels to apply to each column in x. If no named vector is given, the names of the columns of x
#' (other than the column(s) assigned to group_id and id) are used as labels in the output.
#'
#' @return a data frame showing all grouped observations that are statistically different from the
#' expected counts predicted by Benford's Law at the confidence level set in \code{conf.int} and calculated using the
#' method specified in \code{method}. The group_id, label, test statistic, p-value, and group size are reported in the output.
#'
#' @export
#' @references Anderson K.M., Dayaratna K., Gonshorowski D., & Miller S.J. (2022). "A New Benford Test for Clustered Data with Applications to American Elections." \emph{Stats}, 5(3), 841-855. https://doi.org/10.3390/stats5030049
#'
#' @examples
#' benford_test(ohio2020, base=3, group_id="county_name", id="precinct", method="chisq")
benford_test <- function(x, digit=c("first", "second"), base, len=100, group_id=NULL, id=NULL, na.rm=TRUE,
                         conf.int=0.90, method=c("chisq", "multinom"),
                         p.adj=c("BH", "holm", "hochberg", "hommel", "bonferroni", "BY", "fdr", "none"),
                         labs=NULL) {

  df <- x

  if (!(is.character(digit[1])) | !(digit[1] %in% c('first','second'))) {stop("digit must be either 'first' or 'second' ")}
  if (digit[1] == 'first') {dgt <- 1} else {dgt <- 2}

  if (!(is.numeric(base))) {stop("must provide base")}
  if (base < 0) {stop("base must be a positive number")}
  if (base < 3) {stop("base must be 3 or greater")}

  if (!is.data.frame(df)) {try(df <- as.data.frame(df))}

  if (is.null(group_id)) {stop("must assign column name to group_id")}
  if (is.character(group_id)) {
    df <- df %>% relocate(all_of(group_id))
  }

  if (!is.null(id)) {
    if (is.character(id)) {
      df <- df %>% relocate(all_of(id), .after=all_of(group_id))
    }
    if (class(as.character(df[,2]))=="character") {
      df[,2] <- as.character(df[,2])
    } else {stop("id must be a character vector")}
  }

  if (is.null(id)) {
    if (any(is.numeric(df[,-1])) | any(is.character(df[,-1])))
      {warning("numeric values being truncated to zero decimal places and coerced to integer")}
    #if (any(!is.numeric(as.numeric(df[,-1])))) {stop("values cannot be coerced to integer")}
    df <- df %>% mutate_at(-1, function(z) as.integer(trunc(as.numeric(z))))
    if (any(df[,-1] < 0)) {stop("values must be positive")}
    col_id <- grep(group_id, colnames(df))
  } else {
    if (any(is.numeric(df[,-c(1,2)])) | any(is.character(df[,-c(1,2)])))
      {warning("numeric values being truncated to zero decimal places and coerced to integer")}
    #if (any(!is.numeric(as.numeric(df[,-c(1,2)])))) {stop("values cannot be coerced to integer")}
    df <- df %>% mutate_at(-c(1:2), function(z) as.integer(trunc(as.numeric(z))))
    if (any(df[,-c(1:2)] < 0)) {stop("values must be positive")}
    col_id <- c(grep(group_id, colnames(df)), grep(id, colnames(df)))
  }

  pH0 = 0.5     ##
  pHA = 1 - pH0 ##
  benford_probs <- benford_distribution(digit=digit[1], base=base)

  groupVec <- unique(df[,1])
  FinalMat <- data.frame(matrix(nrow=0, ncol=6))
  #BaseConvert <- data.frame(matrix(nrow=0, ncol=ncol(df[,-col_id])))

  for (i in 1:length(groupVec)) {
    y <- subset(df, df[,1] == groupVec[i])

    if (nrow(y) > 10) {
      Names <- y[,col_id]
      y <- y[,-col_id]

      if (na.rm==TRUE) {
        y <- y[complete.cases(y),]
      } else if (any(is.na(y))) {stop("Data frame cannot have missing values. Remove or select na.rm=TRUE")}

      if (any(y==0)) {
        warning("All values must be greater than 0. Dropping cases with zeros.")
        y <- y %>% filter_all(all_vars(. > 0))
      }

      if (is.null(labs)) {
        labs <- colnames(y[,-col_id])
      } else {
        y <- y %>% rename(all_of(labs))
        }

      BaseConvert_i <- data.frame(matrix(nrow=nrow(y), ncol=ncol(y)))
      Digits = data.frame(matrix(nrow=nrow(y), ncol=ncol(y)))
      Exp = benford_probs * nrow(y)

      for (j in 1:ncol(y)) {
        for (k in 1:nrow(y)) {
          BaseConvert_i[k,j] <- as.integer(dec2base(y[k,j], base=base, len=len))
        }

        if (any(BaseConvert_i[,j] == 0))
        {stop("Set len param higher. Nonzero values being converted to numbers with all zeros.")}
      }

      #BaseConvert <- rbind(BaseConvert, BaseConvert_i)
      FinalMat0 <- data.frame(matrix(nrow=ncol(y), ncol=6))

      for (h in 1:ncol(y)) {
        Digits[,(h+dgt)] = extract_digit(BaseConvert_i[,h], digit=dgt)

        Obs <- rep(0, times=base-2+dgt)
        names(Obs) <- seq(from=(2-dgt), to=(base-2+dgt))
        for (k in (2-dgt):length(Obs)) {
            Obs[k] = length(which(Digits[,(h+dgt)]==k))
          }


        if (method == "chisq") {
          test_stat = sum( (Obs - Exp)^2 / Exp )
          pval = pchisq(test_stat, df=base-3+dgt, lower.tail=FALSE)
        } else if (method == "multinom") {
          obs_probs <- Obs / nrow(y)
          test_stat = -2 * sum( Obs[obs_probs!=0]*log(benford_probs[obs_probs!=0]/obs_probs[obs_probs!=0]) )
          cor_test_stat = test_stat / (1 + base/(6*nrow(y)) + ((base-1)^2)/(6*nrow(y)^2) ) # Smith (1981)
          pval = pchisq(cor_test_stat, df=base-3+dgt, lower.tail=FALSE)

          # library(XNomial)
          # if (nrow(y) > 100) {
          #   pval = xmonte(Obs, Exp)$pLLR
          # } else {
          #   pval = xmulti(Obs, Exp)$pLLR}
        }

        BF = Benford_calc(theta = Exp/sum(Exp),
                          x = Obs,
                          alpha = rep(1,length(Exp)))

        FinalMat0[ ,1] = groupVec[i]
        FinalMat0[h,2] = colnames(y)[h]
        FinalMat0[h,3] = test_stat
        FinalMat0[h,4] = pval
        FinalMat0[h,5] = BF
        FinalMat0[h,6] = nrow(y)
      }

      FinalMat <- FinalMat %>% rbind(FinalMat0)
    }
  }


  if (nrow(FinalMat) > 0) {
    names(FinalMat) <- c("Group_ID", "Label", "Test.Statistic", "p.val", "BF", "n")

    FinalMat$pval_adj = p.adjust(FinalMat[,4], method=p.adj[1])
    FinalMat$PostProb = (1 + (pHA/( (pH0)^(1/dim(FinalMat)[1]) )) * (1/exp(FinalMat[,5])) )^-1
    FinalMat$MinPostProb = 1 / (1 + (-exp(1)*FinalMat$pval_adj*log(FinalMat$pval_adj))^-1)

  } else {stop("no groups with 10 or more cases. Test cannot be calculated.")}

  result <- subset(FinalMat, (FinalMat$PostProb < (1 - conf.int) ))

  return(result)

}



# digit = "first"
# base = 3
# current_base = 10
# group_id = "county_name"
# id = "precinct"
# len = 100
# conf.int=0.90
# method = "chisq"
# p.adj = "BH"
# labs <- c("trump"="votes.r", "clinton"="votes.d")
# na.rm=TRUE
#
#
# benford_test(x, base=3, group_id="county_name", id="precinct", method="chisq")

# ohio2020 <- df %>%
#   rename(labs)


# y <- subset(df, df[,1] == groupVec[1])



# x <- x[,c(1,2,4,6)]
#
# is.numeric(x$votes.r)
# v <- x[,c(1,4,6)]
# v$votes.r[1] <- -1*v$votes.r[1]
# x$precinct <- as.factor(x$precinct)
# y[3,1] <- 0

# dec2base(243, base=3)
# max.len <- as.integer(as.character(max(round(log(max(244, 1)) / log(3), digits=50) + 1, 0)))
# pwr <- rep(1, length(243)) * 3^((max.len-1):0)
# jl <- 243 * rep(1, max.len)
# di5 <- floor((jl %% (3 * pwr)) / pwr)
# paste(c(as.character(0:9), LETTERS)[di5 + 1], collapse="")








