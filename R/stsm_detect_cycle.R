#' Detect cycle from the data
#'
#' @param y Univariate time series of data values.
#' @param freq Frequency of the data (1 (yearly), 4 (quarterly), 12 (monthly), 365.25/7 (weekly), 365.25 (daily))
#' @param sig_level Significance level to determine statistically significant seasonal frequencies
#' @param prior A data table created by stsm_prior
#' @import data.table
#' @return Numeric value of cycle periodicity
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
#' cycle = stsm_detect_cycle(y = NA000334Q$y, freq = 4)
#' }
#' @export
stsm_detect_cycle = function(y, freq, sig_level = 0.0001, prior = NULL){
  #Get naive seasonal model
  if(is.null(prior)){
    prior = stsm_prior(y, freq) 
  }else{
    prior = copy(prior)
  }
  
  #Redefine the cycle value
  seasonal = trend = power = harmonic = slope = pval = . = NULL
  prior[, "cycle" := y - seasonal - trend]
  if(tsutils::coxstuart(stats::na.omit(prior$cycle), type = "dispersion")$p.value <= 0.01){
    if(!all(prior$seasonal == 0)){
      #Convert the seasonal component to a multiplicative factor
      prior[, "seasonal" := seasonal/trend + 1]
      prior[, "cycle" := y/(seasonal*trend)]
    }else{
      prior[, "cycle" := y/trend]
    }
  }
  prior[, "t" := 1:.N]
  
  #Wavelet analysis for the cycle
  wave = data.table()
  for(i in c(seq(0.01, 0.99, by = 0.01))){
    #Build the ith harmonic
    prior[, paste0("sin", i) := sin(2*pi*i/freq*t)]
    prior[, paste0("cos", i) := cos(2*pi*i/freq*t)]
    
    #Linear regression for cycle ~ harmonics
    lm = stats::lm(cycle ~ ., prior[, c("cycle", colnames(prior)[grepl("sin|cos", colnames(prior))]), with = FALSE])
    
    #Standardize the coefficient value for the Chi-square test
    #X = B_1^2 + B_2^2 ~ X^2(2)
    coefs = summary(lm)$coef
    coefs = coefs[rownames(coefs) %in% c(paste0("sin", i), paste0("cos", i)), "t value"]
    
    if(all(!is.na(coefs))){
      wave = rbind(wave, data.table(harmonic = i, power = sqrt(sum(coefs^2)), 
                                    pval = stats::pchisq(sum(coefs^2), df = length(coefs), lower.tail = FALSE)))
    }
    suppressWarnings(prior[, c(paste0("sin", i), paste0("cos", i)) := NULL])
  }
  #Find the maxima
  wave[, "slope" := (shift(power, type = "lead", n = 1) - shift(power, type = "lag", n = 1))/(shift(harmonic, type = "lead", n = 1) - shift(harmonic, type = "lag", n = 1))]
  #Find significant maxima
  harmonics = wave[(shift(power, type = "lag", n = 1) < power & shift(power, type = "lead", n = 1) < power) &
                     (shift(sign(slope), type = "lag", n = 1) == 1 & shift(sign(slope), type = "lead", n = 1) == -1), ][
                       pval <= sig_level, ]$harmonic
  #Cycle must be 3 years or longer (use 2.5 for estimation uncertainty) and must have at enough data to observe 2.5 cycles
  cycle.period = wave[harmonic %in% harmonics, ][harmonic > 2.5*freq/length(y) & harmonic < 1/2.5, ][which.max(power), ]$harmonic
  
  return(cycle = freq/cycle.period)
}
