#' Detect seasonality from the data
#'
#' @param y Univariate time series of data values.
#' @param freq Frequency of the data (1 (yearly), 4 (quarterly), 12 (monthly), 365.25/7 (weekly), 365.25 (daily))
#' @param sig_level Significance level to determine statistically significant seasonal frequencies
#' @param full_spectrum check the full spectrum of seasonal frequencies for seasonality
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
stsm_detect_seasonality = function(y, freq, sig_level = 0.0001, prior = NULL, full_spectrum = FALSE){
  #Bind data.table variables to the global environment
  s = seasons = round = diff = trend = cycle = rn = `t value` = power = df = pval = . = NULL
  
  #Return no seasons if the frequency is yearly
  if(freq == 1){
    return(seasons = numeric(0))
  }
  
  #Get the prior
  if(is.null(prior)){
    prior = stsm_prior(y, freq) 
  }else{
    prior = copy(prior)
  }
  
  #Adjust the seasonal value
  prior[, "seasonal" := y - trend]
  if(tsutils::coxstuart(stats::na.omit(prior$seasonal), type = "dispersion")$p.value <= 0.01){
    prior[, "seasonal" := y/trend]
  }
  prior[, "t" := 1:.N]
  
  #Wavelet analysis for seasonality using forward stepwise regression
  if(full_spectrum == FALSE){
    if(floor(freq) == floor(60*60*24*365.25)){
      #Secondly frequency
      #minutely, hourly, daily, weekend, weekday, weekly, monthly, quarterly, semi-yearly, yearly
      spectrum = rev(c(60, 60*60, 60*60*24, 60*60*24*2, 60*60*24*5, 60*60*24*7, freq/12, freq/4, freq/2, freq))
    }else if(floor(freq) == floor(60*60*24*365.25*5/7)){
      #Secondly frequency, weekday only
      #minutely, hourly, daily, weekday, monthly, quarterly, semi-yearly, yearly
      spectrum = rev(c(60, 60*60, 60*60*24, 60*60*24*5, freq/12, freq/4, freq/2, freq))
    }else if(floor(freq) == floor(60*24*365.25)){
      #Minutely frequency
      #hourly, daily, weekend, weekday, weekly, monthly, quarterly, semi-yearly, yearly
      spectrum = rev(c(60, 60*24, 60*24*2, 60*24*5, 60*24*7, freq/12, freq/4, freq/2, freq))
    }else if(floor(freq) == floor(60*24*365.25*5/7)){
      #Minutely frequency, weekday only
      #hourly, daily, weekday, monthly, quarterly, semi-yearly, yearly
      spectrum = rev(c(60, 60*24, 60*24*5, freq/12, freq/4, freq/2, freq)) 
    }else if(floor(freq) == floor(24*365.25)){
      #Hourly frequency
      #daily, weekend, weekday, weekly, monthly, quarterly, semi-yearly, yearly
      spectrum = rev(c(24, 24*2, 24*5, 24*7, freq/12, freq/4, freq/2, freq))
    }else if(floor(freq) == floor(24*365.25*5/7)){
      #Hourly frequency, weekday only
      #daily, weekday, monthly, quarterly, semi-yearly, yearly
      spectrum = rev(c(24, 24*5, freq/12, freq/4, freq/2, freq))
    }else if(floor(freq) == 365){
      #Daily frequency
      #weekend, weekday, weekly, monthly, quarterly, semi-yearly, yearly
      spectrum = rev(c(2, 5, 7, freq/12, freq/4, freq/2, freq))
    }else if(floor(freq) == floor(365.25*5/7)){
      #Daily frequency, weekday only
      #weekly, monthly, quarterly, semi-yearly, yearly
      spectrum = rev(c(5, freq/12, freq/4, freq/2, freq))
    }else if(floor(freq) == 52){
      #Weekly frequency
      #monthly, quarterly, semi-yearly
      spectrum = rev(c(freq/12, freq/4, freq/2, freq))
    }else if(floor(freq) == 12){
      #Monthly frequency
      #quarterly, semi-yearly, yearly
      spectrum = rev(c(freq/4, freq/2, freq))
    }else if(floor(freq) == 4){
      #Quarterly frequency
      #semi-yearly, yearly
      spectrum = rev(c(freq/2, freq))
    }else{
      spectrum = floor(freq):2
      spectrum = c(freq, spectrum[2:length(spectrum)])
      #Must have at least 3 cycles for identification of a seasonal pattern
      spectrum = spectrum[spectrum <= nrow(prior)/3]
    }
  }else{
    spectrum = floor(freq):2
    spectrum = c(freq, spectrum[2:length(spectrum)])
    #Must have at least 3 cycles for identification of a seasonal pattern
    spectrum = spectrum[spectrum <= nrow(prior)/3]
  }
  for(i in unique(spectrum)){
    #Build the ith harmonic
    prior[, paste0("sin", i) := sin(2*pi*1/i*t)]
    prior[, paste0("cos", i) := cos(2*pi*1/i*t)]
  
    #Linear regression for seasonal ~ fourier series
    lm = stats::lm(seasonal ~ ., prior[, c("seasonal", colnames(prior)[grepl("sin|cos", colnames(prior))]), with = FALSE])
    
    #Standardize the coefficient value for the Chi-square test
    #X = B_1^2 + B_2^2 ~ X^2(2)
    coefs = data.table(summary(lm)$coef, keep.rownames = TRUE)[rn != "(Intercept)", ]
    coefs = coefs[grepl("sin|cos", rn), c("rn", "t value")]
    coefs[, "seasons" := gsub("sin|cos", "", rn)]
    wave = coefs[, .(power = sqrt(sum(`t value`^2)), df = .N), by = "seasons"][seasons != "(Intercept)", ]
    wave[, "pval" := stats::pchisq(power^2, df = df, lower.tail = FALSE)]
    if(nrow(coefs[seasons %in% wave[pval > sig_level, ]$seasons, ]) > 0){
      prior[, coefs[seasons %in% wave[pval > sig_level, ]$seasons, ]$rn := NULL]
    }
  }
  seasons = unique(as.numeric(gsub("sin|cos", "", colnames(prior)[grepl("sin|cos", colnames(prior))])))
  
  #Check that the regression F stat is statistically significant
  if(length(colnames(prior)[grepl("sin|cos", colnames(prior))]) > 0){
    lm = stats::lm(seasonal ~ ., prior[, c("seasonal", colnames(prior)[grepl("sin|cos", colnames(prior))]), with = FALSE])
    if(1 - stats::pf(summary(lm)$fstatistic["value"],  df1 = summary(lm)$fstatistic["numdf"], df2 = summary(lm)$fstatistic["dendf"]) <= sig_level){
      seasons = unique(as.numeric(gsub("sin|cos", "", colnames(prior)[grepl("sin|cos", colnames(prior))])))
    }else{
      seasons = numeric(0)
    }
  }else{
    seasons = numeric(0)
  }
  
  return(seasons)
}
