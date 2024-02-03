#' Benford's Law Test
#'
#'@description
#'Identify grouped observations that are statistically inconsistent with Benford's Law.
#'
#'
#' @param x A data frame or an object coercible to a data frame such as a matrix containing values to test for Benfordness
#' and column corresponding to group names. Each column pertains to a response. Columns must contain only positive values.
#' Decimal values will be truncated to zero decimal places and coerced to integer.
#'
#' @param digit A string corresponding to the digit of the values in x to test for Benfordness.
#' User can select either \code{"first"}, \code{"second"}, or \code{"last two"}. If \code{"second"} or \code{"last two"},
#' values with less than 2 digits will automatically be dropped.Currently, only base 10 is supported for option \code{"last two"}.
#'
#' @param base The target base to convert all non-negative integer values in x to. Must be 3 or greater. Required.
#'
#' @param group_id A character variable to group Benford's test on. Groups must contain a minimum of 10 observations.
#' Any group with fewer than 10 observations will be excluded. Observations within groups will
#' be tested for Benfordness using the test specified in \code{method}. Required.
#'
#' @param id A character vector specifying the id for observations.
#'
#' @param method Method used to calculate the test statistic and p-value for the Benford's Law test. User can
#' select "chisq" to compute using Pearson's Chi-Square test or "multinom" to compute using a multinomial likelihood
#' ratio test.
#'
#' @param n.sims Number of Monte Carlo simulations to run for testing group-wise data conformity with discrete weibull
#' distribution and Benford's Law.
#'
#' @param conf.level Confidence level for the test selected in method. Default is \code{conf.level=0.95}.
#'
#' @param p.adj String specifying method for calculating adjusted p-value for multiple comparisons.
#' Function calls \code{stats::p.adjust}. Default method for is Benjamini & Hochberg (1995) ("BH") but any method
#' can be selected. See \code{?stats::p.adjust} for options.
#'
#' @param na.rm How to handle missing values and zeros. If \code{na.rm=TRUE}, NAs and zeros will be
#' excluded from analysis. If \code{na.rm=FALSE} then the user must first remove NAs but observations
#' with zeros will be automatically dropped.
#'
#' @param labs Named vector of labels to apply to each column in x. If no named vector is given, the names of the columns of x
#' (other than the column(s) assigned to group_id and id) are used as labels in the output.
#'
#' @param len The length of the character string used to convert numerals in base 10 to another base.
#' Function calls \code{oro.dicom:dec2base}. It is highly recommended to set the \code{len} > 100 as
#' there is a bug in \code{dec2base} which sets some values with leading 1's to '0'.
#' A higher \code{len} may be required for data with values of higher magnitude.
#'
#' @param verbose \code{verbose=1} prints warnings while \code{verbose=0} silences them.
#'
#'
#' @return A data frame showing all grouped observations that are statistically different from the
#' expected counts predicted by Benford's Law at the confidence level set in \code{conf.level} and calculated using the
#' method specified in \code{method}. The group_id, label, test statistic, p-value, and group size are reported in the output.
#'
#' @import dplyr
#' @import oro.dicom
#'
#' @export
#'
#' @references Anderson K.M., Dayaratna K., Gonshorowski D., & Miller S.J. (2022). "A New Benford Test for Clustered Data with Applications to American Elections." \emph{Stats}, 5(3), 841-855. https://doi.org/10.3390/stats5030049.
#' @references Smith PJ, Rae DS, Manderscheid RW, Silbergeld S. Approximating the moments and distribution of the likelihood ratio statistic for multinomial goodness of fit. Journal of the American Statistical Association. 1981 Sep 1;76(375):737-40.
#'
#' @examples
#' benford_test(ohio2016, base=3, group_id="county_name", id="precinct", method="chisq")
#'
benford_test <- function(x, digit=c("first", "second", "last two"), base, group_id=NULL, id=NULL, na.rm=TRUE,
                         method=c("chisq", "multinom"), n.sims=100, conf.level=0.95,
                         p.adj=c("BH", "holm", "hochberg", "hommel", "bonferroni", "BY", "fdr", "none"),
                         labs=NULL, len=100, verbose=0) {

  df <- x

  if (!(is.character(digit[1])) | !(digit[1] %in% c('first','second','last two'))) {stop("digit must be either 'first', 'second', or 'last two' ")}
  if (digit[1] == 'first') {dgt <- 1} else {dgt <- 2}

  if (!(is.numeric(base))) {stop("must provide base")}
  if (base < 0) {stop("base must be a positive number")}
  if (base < 3) {stop("base must be 3 or greater")}
  if (digit[1]=='last two') {if (base!=10) {stop("Only base 10 is currently supported for the last two digits")}}

  if (!is.data.frame(df)) {try(df <- as.data.frame(df))}

  if (is.null(group_id)) {stop("must assign column name to group_id")}
  if (is.character(group_id)) {
    df <- df %>% relocate(all_of(group_id))
  }

  if (!is.null(id)) {
    if (is.character(id)) {
      df <- df %>% relocate(all_of(id), .after=all_of(group_id))
    }
    if (!class(eval(parse(text=paste0('df$', id))))=="character") {
      if (class(as.character(df[,2]))=="character") {
        df[,2] <- as.character(df[,2])
      } else {stop("id must be a character vector")}
    }
  }

  if (is.null(id)) {
    if ((any(is.numeric(df[,-1])) | any(is.character(df[,-1]))) & verbose==1)
      {warning("numeric values being truncated to zero decimal places and coerced to integer")}
    #if (any(!is.numeric(as.numeric(df[,-1])))) {stop("values cannot be coerced to integer")}
    df <- df %>% mutate_at(-1, function(z) as.integer(trunc(as.numeric(z))))
    if (any(df[,-1] < 0)) {stop("values must be positive")}
    col_id <- grep(group_id, colnames(df))
  } else {
    if ((any(is.numeric(df[,-c(1,2)])) | any(is.character(df[,-c(1,2)]))) & verbose==1)
      {warning("numeric values being truncated to zero decimal places and coerced to integer")}
    #if (any(!is.numeric(as.numeric(df[,-c(1,2)])))) {stop("values cannot be coerced to integer")}
    df <- df %>% mutate_at(-c(1:2), function(z) as.integer(trunc(as.numeric(z))))
    if (any(df[,-c(1:2)] < 0)) {stop("values must be positive")}
    col_id <- c(grep(group_id, colnames(df)), grep(id, colnames(df)))
  }

  pH0 = 0.5     ##
  pHA = 1 - pH0 ##

  if (digit[1] %in% c('first','second')) {
    benford_probs <- benford_distribution(digit=digit[1], base=base)
  } else if (digit[1]=='last two') {
    benford_probs <- c(10/99, 89/99)
    names(benford_probs) <- c('A', 'N')
  }

  groupVec <- as.vector(unlist(unique(df[,1])))
  FinalMat <- data.frame(matrix(nrow=0, ncol=6))

  for (i in 1:length(groupVec)) {
    y <- subset(df, df[,1] == groupVec[i])

    if (nrow(y) > 10) {
      Names <- y[,col_id]
      y <- y[,-col_id]

      if (na.rm==TRUE) {
        y <- y[complete.cases(y),]
      } else if (any(is.na(y))) {stop("Data frame cannot have missing values. Remove or select na.rm=TRUE")}

      if (any(y==0)) {
        if (verbose==1) {warning("All values must be greater than 0. Dropping cases with zeros.")}
        y <- y %>% filter_all(all_vars(. > 0))
      }

      if (digit[1] %in% c('second','last two')) {
        if (verbose==1) {warning("All values must have 2 digits or more. Dropping ids where at least one observation has less than 2 digits.")}
        y <- y %>% filter_all(all_vars(nchar(.) > 1))
      }

      if (!is.null(labs)) {y <- y %>% rename(all_of(labs))}

      BaseConvert_i <- data.frame(matrix(nrow=nrow(y), ncol=ncol(y)))
      Digits = data.frame(matrix(nrow=nrow(y), ncol=ncol(y)))

      Exp = diff(c(0, round(cumsum(benford_probs * nrow(y)))))   #### Should we round expected values?
      if (sum(Exp) != nrow(y)) {stop("Sum of rounded expected values from Benford's disribution not equal to length of group_id")}

      for (j in 1:ncol(y)) {
        for (k in 1:nrow(y)) {
          BaseConvert_i[k,j] <- as.integer(dec2base(as.integer(y[k,j]), base=base, len=len))
        }

        if (any(BaseConvert_i[,j] == 0))
        {stop("Set len param higher. Nonzero values being converted to numbers with all zeros.")}
      }


      FinalMat0 <- data.frame(matrix(nrow=ncol(y), ncol=6))

      for (h in 1:ncol(y)) {
        if (digit[1] %in% c('first','second')) {
          Digits[,h] = extract_digit(BaseConvert_i[,h], digit=dgt)
        } else if (digit[1]=='last two') {
          Digits[,h] = extract_digit(BaseConvert_i[,h], number.of.digits=2, last=TRUE)
        }

        Obs <- rep(0, times=length(benford_probs))
        names(Obs) <- names(benford_probs)

        if (digit[1] %in% c('first','second')) {
          for (g in (2-dgt):length(Obs)) {
            Obs[g] = length(which(Digits[,h]==g))}
        } else if (digit[1]=='last two') {
            Obs[1] = length(which(Digits[,h] %in% c(0, seq(11, 99, by=11))))
            Obs[2] = length(Digits[,h]) - Obs[1]
        }

        if (digit[1] %in% c('first','second')) {dgts <- dgt} else if (digit[1]=='last two') {dgts = -base + 4}


        if (method == "chisq") {
          test_stat = sum( (Obs - Exp)^2 / Exp )
          pval = stats::pchisq(test_stat, df=base-3+dgts, lower.tail=FALSE)
        } else if (method == "multinom") {
          obs_probs <- Obs / nrow(y)
          test_stat <- 2 * sum(Obs[Obs!=0] * log(Obs[Obs!=0] / (nrow(y)*benford_probs[Obs!=0])))
          pval = stats::pchisq(test_stat, df=base-3+dgts, lower.tail=FALSE)

          # cor_test_stat = test_stat / (1 + base/(6*nrow(y)) + ((base-1)^2)/(6*nrow(y)^2) ) # Smith (1981)

          # library(XNomial)
          # if (nrow(y) > 100) {
          #   pval = xmonte(Obs, Exp)$pLLR
          # } else {
          #   pval = xmulti(Obs, Exp)$pLLR}
        }

        BF = Benford_calc(theta = Exp/sum(Exp),
                          x = Obs,
                          alpha = rep(1,length(Exp)))

        FinalMat0[h,1] = groupVec[i]
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

    FinalMat$pval_adj = stats::p.adjust(FinalMat[,4], method=p.adj[1])
    FinalMat$PostProb = (1 + (pHA/( (pH0)^(1/dim(FinalMat)[1]) )) * (1/exp(FinalMat[,5])) )^-1
    FinalMat$MinPostProb = 1 / (1 + (-exp(1)*FinalMat$pval_adj*log(FinalMat$pval_adj))^-1)

  } else {stop("no groups with 10 or more cases. Test cannot be calculated.")}

  result <- subset(FinalMat, (FinalMat$pval_adj < (1 - conf.level) ))




  anom_counties <- result %>%
    dplyr::select(Group_ID, Label)

  df2 <- df %>% filter(eval(parse(text=group_id)) %in% anom_counties$Group_ID)
  bf_valid_test <- data.frame(matrix(nrow=0, ncol=7,
                                     dimnames = list(NULL,
                                                     c('Group_ID', 'Label', 'Scale', 'Shape', 'N', 'D.weibull_KS_pval', 'BF_pval'))))
  bf_valid_test[,1:2] <- lapply(bf_valid_test[,1:2], as.character)
  bf_valid_test[,3:7] <- lapply(bf_valid_test[,1:2], as.numeric)


  for (z in 1:length(unique(anom_counties$Group_ID))) {
    group_id_i <- unique(anom_counties$Group_ID)[z]

    cty_df <- df2 %>%
      rename(group_id = all_of(group_id)) %>%
      filter(group_id == group_id_i)

    if (!is.null(id)) {
      cty_df <- cty_df %>%
        rename(id = all_of(id)) %>%
        group_by(id) %>%
        #filter(!grepl("ABS", mode, ignore.case=TRUE)) %>%
        filter(!grepl("TRANSFER", id, ignore.case=TRUE)) %>%
        ungroup() %>%
        group_by(group_id, id) %>%
        summarise_all(sum, na.rm=TRUE) %>%
        ungroup()
    }

    for (w in 2:(ncol(df2)-1)) {

      if (!is.null(id)) {
        cty_df_i <- cty_df %>%
          dplyr::select(group_id, id, w+1)
      } else {
        cty_df_i <- cty_df %>%
          dplyr::select(group_id, w)
      }

      cty_df_i <- cty_df_i[complete.cases(cty_df_i),]

      if (any(cty_df_i==0)) {
        if (verbose==1) {warning("All values must be greater than 0. Dropping cases with zeros.")}
        cty_df_i <- cty_df_i %>% filter_all(all_vars(. > 0))
      }

      if (!is.null(id)) {data_idx <- 3} else {data_idx <- 2}

      datavec <- as.numeric(unlist(cty_df_i[,data_idx]))
      datavec <- datavec[!is.na(datavec)]

      if (digit[1] %in% c('second','last two')) {
        if (any(nchar(datavec) < 2)) {
          if (verbose==1) {warning("All values must have 2 digits or more. Dropping cases with less than 2 digits.")}
        }
        datavec <- datavec[nchar(datavec) > 1]
      }

      n_precincts <- length(datavec)


      ###########################
      # Estimate Parameters
      ###########################
      if (length(datavec) > 1) {

        #source("lldweibull_Kevin.R")
        #---------------------------------------------------------------------------------------------
        loglike <- function(p) suppressWarnings({-sum(log(d.weibull(datavec, p[1], p[2])))})
        mle <- stats::optim(c(2000,2), loglike, control=list(maxit=1e5))  # NB: Changed the initial value for alpha #********

        ks_pval <- rep(0, n.sims)
        simvec <- list()

        for (k in 1:n.sims) {
          simvec[[k]] = rd.weibull(n_precincts, abs(mle$par[1]), abs(mle$par[2]))
          ks_result = suppressWarnings(stats::ks.test(datavec, simvec[[k]]))
          ks_pval[k] = ks_result$p.value
        }

        ks_pval_avg = mean(ks_pval)
        #---------------------------------------------------------------------------------------------

      } else {
        mle <- list(par=c(NA, NA), convergence=NA)
        ks_pval_avg = NA
      }


      ##################################################################
      # Test MLE of Individual-Level Discrete Weibull for Benfordness
      ##################################################################
      if (length(datavec) > 1) {
        Exp2 = diff(c(0, round(cumsum(benford_probs * n_precincts))))
        if (sum(Exp2) != n_precincts) {stop("Sum of rounded expected values from Benford's disribution not equal to length of group_id")}

        BF_pval <- c()

        for (v in 1:n.sims) {
          BaseConvert_2 <- c()

          for (l in 1:n_precincts) {
            BaseConvert_2[l] <- as.integer(dec2base(as.integer(simvec[[v]][l]), base=base, len=len))
          }

          BaseConvert_2 <- BaseConvert_2[BaseConvert_2 > 0]

          if (digit[1] %in% c("second","last two")) {
            BaseConvert_2 <- BaseConvert_2[BaseConvert_2 > 9]
          }

          if (digit[1] %in% c('first','second')) {
            Digits_2 = extract_digit(BaseConvert_2, digit=dgt)
          } else if (digit[1]=='last two') {
            Digits_2 = extract_digit(BaseConvert_2, number.of.digits=2, last=TRUE)
          }

          Obs2 <- rep(0, times=length(benford_probs))
          names(Obs2) <- names(benford_probs)

          if (digit[1] %in% c('first','second')) {
            for (g in (2-dgt):length(Obs2)) {
              Obs2[g] = length(which(Digits_2==g))}
          } else if (digit[1]=='last two') {
            Obs2[1] = length(which(Digits_2 %in% c(0, seq(11, 99, by=11))))
            Obs2[2] = length(Digits_2) - Obs2[1]
          }


          if (method == "chisq") {
            test_stat = sum( (Obs2 - Exp2)^2 / Exp2 )
            BF_pval[v] = stats::pchisq(test_stat, df=base-3+dgts, lower.tail=FALSE)

          } else if (method == "multinom") {
            test_stat2 <- 2 * sum(Obs2[Obs2!=0] * log(Obs2[Obs2!=0] / (n_precincts*benford_probs[Obs2!=0])))
            BF_pval[v] = stats::pchisq(test_stat2, df=base-3+dgts, lower.tail=FALSE)

            # cor_test_stat = test_stat / (1 + base/(6*n_obs) + ((base-1)^2)/(6*n_obs^2) ) # Smith (1981)
          }
        }

        BF_pval_avg = mean(BF_pval)

        bf_valid_test_row <- data.frame(Group_ID = group_id_i,
                                        Label = utils::tail(colnames(cty_df_i), n=1),
                                        Scale = mle$par[1],
                                        Shape = mle$par[2],
                                        N = n_precincts,
                                        D.weibull_KS_pval = ks_pval_avg,
                                        BF_pval = BF_pval_avg)

        bf_valid_test <- rbind(bf_valid_test, bf_valid_test_row)
      }
    }
  }

  bf_valid_test$BF_adj_pval <- stats::p.adjust(bf_valid_test$BF_pval, method=p.adj[1])

  bf_valid_test$conformity <- ifelse((bf_valid_test$BF_adj_pval > 0.05) & (bf_valid_test$D.weibull_KS_pval > 0.05), '*', '')

  result <- result %>%
    left_join(bf_valid_test, by=c("Group_ID", "Label"))

  return(result)

}


# round(bf_valid_test[,7], digits=3)
# table(round(p.adjust(bf_valid_test[,7], method=p.adj[1]), digits=3))


#
# x <- ohio2016
# digit = "first"
# base = 4
# group_id = "county_name"
# id = "precinct"
# len = 100
# conf.level=0.95
# method = "chisq"
# n.sims=100
# p.adj = "BH"
# na.rm=TRUE
# verbose = 0
# labs = NULL


#

# benford_test(x, base=3, group_id="county_name", id="precinct", method="chisq",
#              labs = c("trump"="votes.r", "clinton"="votes.d"))
#
# # ohio2020 <- df %>%
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








