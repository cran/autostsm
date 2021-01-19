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
  
  seasonal_adj = trend = NULL
  if(all(stats::na.omit(y) > 0)){
    if(any(is.na(y))){
      y = suppressWarnings(imputeTS::na_kalman(y))
    }
    if(is.null(prior)){
      prior = stsm_prior(y, freq) 
    }else{
      prior = copy(prior)
    }
    prior[, "seasonal" := y - trend]
    
    if(!all(prior$seasonal == 0)){
      #Test for increasing/decreasing seasonal amplitude
      multiplicative = (multiplicative | tsutils::coxstuart(prior$seasonal, type = "dispersion")$p.value <= sig_level)
    }
    
    #PE test for non-nested models and functional form choice: linear vs log model
    prior[, "t" := 1:.N]
    if(any(prior$seasonal_adj < 0)){
      prior[, "seasonal_adj" := seasonal_adj + abs(min(seasonal_adj, na.rm = TRUE)) + 1]
    }
    prior = prior[!is.na(seasonal_adj), ]
    lm_lin = stats::lm(seasonal_adj ~ t, prior)
    lm_log = stats::update(lm_lin, log(seasonal_adj) ~ t)
    stat = lmtest::petest(lm_lin, lm_log)$`Pr(>|t|)`
    #To select log model, the linear model must be rejected while the log model must fail to be rejected
    multiplicative = (multiplicative | (stat[1] <= sig_level & stat[2] > sig_level))
  }
  return(multiplicative)
}
