#' Detect trend type
#'
#' @param y Univariate time series of data values. May also be a 2 column data frame containing a date column.
#' @param freq Frequency of the data (1 (yearly), 4 (quarterly), 12 (monthly), 365.25/7 (weekly), 365.25 (daily))
#' @param decomp Decomposition model ("tend-cycle-seasonal", "trend-seasonal", "trend-cycle", "trend-noise")
#' @param sig_level Significance level to determine statistically significant seasonal frequencies
#' @param prior A data table created by stsm_prior
#' @param seasons The seasonal periods
#' @param cycle The cycle period
#' @import data.table
#' @return list with trend type and logical flag for deterministic trend if the trend is determined to have 0 differencing
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
#' trend = stsm_detect_trend(y = NA000334Q$y, freq = 4)
#' }
#' @export
stsm_detect_trend = function(y, freq, decomp = "", sig_level = 0.01, prior = NULL, seasons = NULL, cycle = NULL){
  #Bind data.table variables to the global environment
  N = d = test = pval = sig = trend = remainder = seasonal = seasonal_adj = NULL
  
  #Set the prior
  if(is.null(prior)){
    prior = stsm_prior(y, freq, decomp , seasons, cycle) 
  }else{
    prior = copy(prior)
  }
  prior[, "seasonal_adj" := y - seasonal - cycle]
  
  ##### Find the number of differences to make it stationary #####
  #Replace outliers for the trend test
  prior[, "seasonal_adj" := forecast::tsclean(seasonal_adj, replace.missing = FALSE)]
  #Test for a trend
  trend_test = (stsm_coxstuart(stats::na.omit(prior$seasonal_adj), type = "trend")$p.value <= sig_level)
  #Calculate unit root tests for differences 0 to 2
  if(trend_test == FALSE){
    ndiffs = c(forecast::ndiffs(prior$seasonal_adj, test = "adf", type = "level", alpha = sig_level, max.d = 2), 
             forecast::ndiffs(prior$seasonal_adj, test = "pp", type = "level", alpha = sig_level, max.d = 2), 
             forecast::ndiffs(prior$seasonal_adj, test = "kpss", type = "level", alpha = sig_level, max.d = 2))
    ndiffs = round(mean(ndiffs, na.rm = T))
  }else{
    ndiffs = c(forecast::ndiffs(diff(prior$seasonal_adj), test = "adf", type = "level", alpha = sig_level, max.d = 1), 
                forecast::ndiffs(diff(prior$seasonal_adj), test = "pp", type = "level", alpha = sig_level, max.d = 1), 
                forecast::ndiffs(diff(prior$seasonal_adj), test = "kpss", type = "level", alpha = sig_level, max.d = 1))
    ndiffs = round(mean(ndiffs, na.rm = T)) + 1
  }
  
  if(ndiffs >= 2){
    trend = "double-random-walk"
  }else if(ndiffs == 1 & trend_test == TRUE){
    trend = "random-walk-drift"
  }else if(ndiffs == 1 & trend_test == FALSE){
    trend = "random-walk"
  }else if(ndiffs == 0 & trend_test == FALSE){
    trend = "random-walk"
  }else if(ndiffs == 0 & trend_test == TRUE){
    trend = "random-walk-drift"
  }
  
  return(list(trend = trend, det_trend = (ndiffs == 0)))
}
