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
  #Bind data.table variables to the global environment
  seasonal_adj = trend = seasonal = seasonal2 = cycle = cycle2 = NULL
  multiplicative = FALSE
  
  if(all(stats::na.omit(y) > 0)){
    if(is.null(prior)){
      prior = stsm_prior(y, freq) 
    }else{
      prior = copy(prior)
    }
      
    if(!all(prior[!is.na(seasonal), ]$seasonal == 0)){
      #Test for increasing/decreasing seasonal amplitude
      #Remove outliers for sensitivity
      prior[, "seasonal2" := y - trend - cycle]
      prior[, "seasonal2" := forecast::tsclean(seasonal2, replace.missing = FALSE)]
      multiplicative = (multiplicative | round(stsm_coxstuart(stats::na.omit(prior$seasonal2), type = "deviation")$p.value, 2) <= sig_level)
      prior[, "seasonal2" := NULL]
    }
    
    if(!all(prior[!is.na(cycle), ]$cycle == 0)){
      #Test for increasing/decreasing cycle amplitude
      #Remove outliers for sensitivity
      prior[, "cycle2" := y - trend - seasonal]
      prior[, "cycle2" := forecast::tsclean(cycle2, replace.missing = FALSE)]
      multiplicative = (multiplicative | round(stsm_coxstuart(stats::na.omit(prior$cycle2), type = "deviation")$p.value, 2) <= sig_level)
      prior[, "cycle2" := NULL]
    }
    
    #PE test for non-nested models and functional form choice: linear vs log model
    prior[, "t" := 1:.N][, "seasonal_adj" := y - seasonal - cycle]
    if(any(prior[!is.na(seasonal_adj), ]$seasonal_adj < 0)){
      prior[, "seasonal_adj" := seasonal_adj + abs(min(seasonal_adj, na.rm = TRUE)) + 1]
    }
    #Remove outliers for sensitivity
    prior[, "seasonal_adj" := forecast::tsclean(seasonal_adj, replace.missing = FALSE)]
    lm_lin = stats::lm(seasonal_adj ~ t, prior[!is.na(seasonal_adj), ])
    lm_log = stats::update(lm_lin, log(seasonal_adj) ~ t)
    stat = lmtest::petest(lm_lin, lm_log, vcov. = sandwich::vcovHAC)$`Pr(>|t|)`
    #To select log model, the linear model must be rejected while the log model must fail to be rejected
    multiplicative = (multiplicative | (round(stat[1], 2) <= sig_level & round(stat[2], 2) > sig_level))
    prior[, "seasonal_adj" := NULL]
  }
  return(multiplicative)
}
