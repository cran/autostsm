#' Return a naive model prior decomposition
#'
#' @param y an object created from stsm_detect_frequency
#' @param freq Frequency of the data
#' @param harmonics The harmonics to split the seasonality into
#' @param decomp decomposition string
#' @import data.table
#' @return data table containing a naive decomposition using STL
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
#' prior = stsm_prior(y = NA000334Q$y, freq = 4)
#' }
#' @export
stsm_prior = function(y, freq, decomp = "", harmonics = NULL){
  if(any(is.na(y))){
    y = suppressWarnings(imputeTS::na_kalman(y))
  }
  
  trend = seasonal = cycle = remainder = NULL
  if(length(y) > 2*freq){
    if(seastests::wo(stats::ts(y, frequency = floor(freq)))$stat == TRUE | grepl("seasonal", decomp)){
      prior = data.table(stats::stl(stats::ts(y, frequency = floor(freq)), s.window = 7, s.degree = 1, l.degree = 1, t.degree = 1)$time.series)
    }else{
      prior = data.table(trend = stats::predict(stats::loess(y ~ t, data.table(y = y)[, "t" := 1:.N])))[, "seasonal" := 0][, "remainder" := y - trend]
    }
  }else{
    prior = data.table(trend = stats::predict(stats::loess(y ~ t, data.table(y = y)[, "t" := 1:.N])))[, "seasonal" := 0][, "remainder" := y - trend]
  }
  if(all(prior$seasonal == 0) & ifelse(!is.null(harmonics), length(harmonics) > 0, FALSE)){
    prior[, "seasonal" := y - trend][, "remainder" := y - trend - seasonal][, "seasonal_adj" := y - seasonal]
  }
  prior[, "seasonal_adj" := y - seasonal]
  if(grepl("cycle", decomp) | decomp == ""){
    if(!all(prior$seasonal == 0)){
      new_trend = stats::predict(stats::loess(y ~ t, data.table(y = y)[, "t" := 1:.N]))
      prior[, "cycle" := trend - new_trend]
      prior[, "trend" := new_trend]
    }else{
      prior[, "cycle" := stats::predict(stats::smooth.spline(remainder))$y]
      prior[, "remainder" := y - seasonal - cycle]
    }
  }else{
    prior[, "cycle" := 0]
  }
  prior[, "drift" := trend - shift(trend, type = "lag", n = 1)]
  prior[1, "drift" := 0]
  prior[, "remainder" := y - trend - seasonal - cycle]
  if(ifelse(!is.null(harmonics), length(harmonics) > 0, FALSE)){
    ff = data.table(seasonal = prior$seasonal, t = 1:nrow(prior))
    for(j in harmonics){
      ff[, paste0("sin", j) := sin(2*pi*j/freq*t)]
      ff[, paste0("cos", j) := cos(2*pi*j/freq*t)]
    }
    ff[, "t" := NULL]
    lm_s = stats::lm(seasonal ~ ., data = ff)
    prior = cbind(prior, do.call("cbind", lapply(harmonics, function(x){
      matrix(stats::coef(lm_s)["(Intercept)"]/length(harmonics) + 
               c(as.matrix(ff[, paste0(c("sin", "cos"), x), with = FALSE]) %*% 
                   matrix(stats::coef(lm_s)[paste0(c("sin", "cos"), x)], ncol = 1)) + 
               lm_s$residuals/length(harmonics), 
             ncol = 1)
    })))
    colnames(prior)[(ncol(prior) - length(harmonics) + 1):ncol(prior)] = paste0("seasonal", floor(harmonics))
  }
  return(prior)
}
