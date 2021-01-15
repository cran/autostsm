#' Detect if log transformation is best
#'
#' @param y an object created from stsm_detect_frequency
#' @param freq Frequency of the data
#' @param sig_level Significance level to determine statistically significant seasonal frequencies
#' @param prior A data table created by stsm_prior
#' @import data.table
#' @return a logical indicating if the model should be multiplicative or not
#' @examples
#' \dontrun{
#' #GDP Not seasonally adjusted
#' library(autostsm)
#' data("NA000334Q", package = "autostsm") #From FRED
#' NA000334Q = data.table(NA000334Q, keep.rownames = TRUE)
#' colnames(NA000334Q) = c("date", "y")
#' NA000334Q[, "date" := as.Date(date)]
#' NA000334Q[, "y" := as.numeric(y)]
#' NA000334Q = NA000334Q[date >= "1990-01-01", ]
#' multiplicative = stsm_detect_multiplicative(y = NA000334Q$y, freq = 4)
#' }
#' @export
stsm_detect_multiplicative = function(y, freq, sig_level = 0.01, prior = NULL){
  multiplicative = FALSE
  
  if(all(stats::na.omit(y) > 0)){
    if(any(is.na(y))){
      y = suppressWarnings(imputeTS::na_kalman(y))
    }
    if(is.null(prior)){
      prior = stsm_prior(y, freq) 
    }else{
      prior = copy(prior)
    }
    
    if(!all(prior$seasonal == 0)){
      #Test for increasing/decreasing seasonal amplitude
      multiplicative = (multiplicative | tsutils::coxstuart(prior$seasonal, type = "deviation")$p.value <= sig_level)
    }
    
    #Test for linear vs exponential trend
    prior[, "t" := 1:.N]
    seasonal_adj = NULL
    if(any(prior$seasonal_adj < 0)){
      prior[, "seasonal_adj" := seasonal_adj + abs(min(seasonal_adj, na.rm = TRUE)) + 1]
    }
    prior = prior[!is.na(seasonal_adj), ]
    lm_lin = stats::lm(seasonal_adj ~ t, prior)
    lm_log = stats::lm(log(seasonal_adj) ~ t, prior)
    lm_log = stats::lm(prior$seasonal_adj ~ offset(exp(lm_log$fitted.values)))
    #Compare LnL (reduced form AIC because same number of parameters)
    multiplicative = (multiplicative | (stats::logLik(lm_log) > stats::logLik(lm_lin))) 
  }
  return(multiplicative)
}
