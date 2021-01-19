#' Detect seasonality from the data
#'
#' @param y Univariate time series of data values.
#' @param freq Frequency of the data (1 (yearly), 4 (quarterly), 12 (monthly), 365.25/7 (weekly), 365.25 (daily))
#' @param sig_level Significance level to determine statistically significant seasonal frequencies
#' @param full Whether to check for full spectrum of seasonality rather than weekly and yearly suggestion
#' @param spectrum, Vector of seasonal periods to check
#' @param prior A data table created from stsm_prior
#' @import data.table
#' @return Numeric vector of seasonal periodicities
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
#' seasonality = stsm_detect_seasonality(y = NA000334Q$y, freq = 4)
#' }
#' @export
stsm_detect_seasonality = function(y, freq, sig_level = 0.0001, prior = NULL, full = F, spectrum = NULL){
  #Get naive model
  if(is.null(prior)){
    prior = stsm_prior(y, freq) 
  }else{
    prior = copy(prior)
  }
  
  #Adjust the seasonal value
  s = round = diff = trend = cycle = rn = `t value` = power = df = pval = . = NULL
  prior[, "seasonal" := y - trend - cycle]
  if(tsutils::coxstuart(stats::na.omit(prior$seasonal), type = "dispersion")$p.value <= 0.01){
    prior[, "cycle" := cycle/trend + 1]
    prior[, "seasonal" := y/(cycle*trend)]
  }
  prior[, "t" := 1:.N]
  
  #Wavelet analysis for seasonality using forward stepwise regression starting from harmonic 1 to floor(freq/2)
  if(full == F & is.null(spectrum)){
    if(floor(freq) == 365){
      #weekend, weekday, weekly, yearly
      spectrum = c(2, 5, 7, freq)
    }else if(floor(freq) == 260){
      #weekly, yearly
      spectrum = c(5, freq)
    }else if(floor(freq) == 52){
      #yearly
      spectrum = c(freq)
    }else if(floor(freq) == 12){
      #yearly
      spectrum = c(freq) 
    }else if(floor(freq) == 4){
      #semi-yearly, yearly
      spectrum = c(freq) 
    }
  }else{
    spectrum = c(2:floor(freq), freq)
  }
  for(i in rev(spectrum)){
    #Build the ith harmonic
    prior[, paste0("sin", freq/i) := sin(2*pi*1/i*t)]
    prior[, paste0("cos", freq/i) := cos(2*pi*1/i*t)]
    
    #Linear regression for cycle ~ harmonics
    lm = stats::lm(seasonal ~ ., prior[, c("seasonal", colnames(prior)[grepl("sin|cos", colnames(prior))]), with = FALSE])
    
    #Standardize the coefficient value for the Chi-square test
    #X = B_1^2 + B_2^2 ~ X^2(2)
    coefs = data.table(summary(lm)$coef, keep.rownames = TRUE)[rn != "(Intercept)", ]
    coefs = coefs[grepl("sin|cos", rn), c("rn", "t value")]
    coefs[, "harmonics" := gsub("sin|cos", "", rn)]
    wave = coefs[, .(power = sqrt(sum(`t value`^2)), df = .N), by = "harmonics"][harmonics != "(Intercept)", ]
    wave[, "pval" := stats::pchisq(power^2, df = df, lower.tail = FALSE)]
    if(nrow(coefs[harmonics %in% wave[pval > sig_level, ]$harmonics, ]) > 0){
      prior[, coefs[harmonics %in% wave[pval > sig_level, ]$harmonics, ]$rn := NULL]
    }
  }
  harmonics = unique(as.numeric(gsub("sin|cos", "", colnames(prior)[grepl("sin|cos", colnames(prior))])))
  suppressWarnings(prior[, colnames(prior)[grepl("sin|cos", colnames(prior))] := NULL])
  
  if(freq > 1){
    #Double check seasonality at the maximum frequency. Sometimes the Fourier analysis will not detect
    #seasonality at the frequency of the data if there is not enough observations
    if(length(y) > 2*freq){
      if(seastests::wo(stats::ts(suppressWarnings(imputeTS::na_kalman(y)), frequency = floor(freq)))$stat){
        seasonal.periods = unique(c(freq, freq/harmonics))
      }else{
        seasonal.periods = unique(freq/harmonics)
      }
    }else{
      seasonal.periods = unique(freq/harmonics)
    }
  }
  
  #Reduce the dimensionality
  seasonal.periods = data.table(s = seasonal.periods)[, "round" := round(s)][, "diff" := abs(s - round)]
  seasonal.periods = seasonal.periods[seasonal.periods[, .I[which.min(diff)], by = "round"]$V1, ]$s
  
  return(seasons = unique(seasonal.periods))
}
