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
  N = d = test = pval = sig = trend = remainder = NULL
  
  #Set the prior
  if(is.null(prior)){
    prior = stsm_prior(y, freq, decomp, seasons, cycle) 
  }else{
    prior = copy(prior)
  }
  prior[, "seasonal_adj" := trend + remainder]
  
  ##### Find the number of differences to make it stationary #####
  #Replace outliers for the trend test
  ol = forecast::tsoutliers(stats::ts(prior$seasonal_adj, frequency = freq))
  prior[ol$index, "seasonal_adj" := ol$replacements]
  #Test for a trend
  trend_test = (tsutils::coxstuart(stats::na.omit(prior$seasonal_adj), type = "trend")$p.value <= sig_level)
  #Calculate unit root tests for differences 0 to 2
  ndiffs = data.table(d = rep(0:2, 3), test = unlist(lapply(c("adf", "pp", "kpss"), rep, 3)))
  for(nd in ndiffs$d){
    if(nd == 0){
      x = stats::na.omit(prior$seasonal_adj)
    }else{
      x = stats::na.omit(diff(prior$seasonal_adj, differences = nd))
    }
    suppressWarnings(ndiffs[d == nd & test == "adf", "pval" := tseries::adf.test(x, alternative = "stationary")$p.value])
    suppressWarnings(ndiffs[d == nd & test == "pp", "pval" := tseries::pp.test(x, alternative = "stationary")$p.value])
    suppressWarnings(ndiffs[d == nd & test == "kpss", "pval" := tseries::kpss.test(x, null = "Level")$p.value])
  }
  ndiffs[test %in% c("adf", "pp"), "sig" := pval <= sig_level]
  ndiffs[test == "kpss", "sig" := pval > sig_level]
  ndiffs = ndiffs[ndiffs[sig == TRUE, .I[which.min(pval)], by = "test"]$V1, ]
  if(nrow(ndiffs) == 0){
    ndiffs = 2
  }else{
    #Find the most agreed upon differences
    ndiffs = ndiffs[, .N, by = "d"][N == max(N), ]
    #If there is a tie, default to 1 difference if exists, then 2, then the most frequent which is likely to be 0
    ndiffs = ifelse(1 %in% ndiffs$d, 1, ifelse(2 %in% ndiffs$d, 2, ndiffs[1, ]$d))
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
