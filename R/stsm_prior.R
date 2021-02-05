#' Return a naive model prior decomposition
#'
#' @param y an object created from stsm_detect_frequency
#' @param freq Frequency of the data
#' @param seasons The seasonal periods to split the seasonality into
#' @param decomp decomposition string
#' @param cycle The cycle periods
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
stsm_prior = function(y, freq, decomp = "", seasons = NULL, cycle = NULL){
  #Bind data.table global variables
  trend = seasonal = cycle = remainder = seasonalcycle = NULL
  
  #Impute missing values for using the Kalman filter for prior decompoosition
  if(any(is.na(y))){
    y = suppressWarnings(imputeTS::na_kalman(y))
  }
  
  #Split the data into trend + seasonal-cycle + remainder
  prior = data.table(trend = stats::predict(stats::loess(y ~ t, data.table(y = y)[, "t" := 1:.N])))[, "remainder" := y - trend]
  prior[, "seasonalcycle" := stats::predict(stats::smooth.spline(remainder))$y][, "remainder" := y - trend - seasonalcycle]
  prior[, "seasonalcycle" := seasonalcycle + remainder*2/3][, "remainder" := remainder*1/3]
  prior[, "seasonal" := 0][, "cycle" := 0]
  
  #Get the drift and 
  prior[, "drift" := trend - shift(trend, type = "lag", n = 1)]
  prior[1, "drift" := 0]
  
  #Get a smooth cycle estimate
  if(grepl("cycle", decomp) | decomp == "" | 
     (ifelse(!is.null(cycle), length(cycle) > 0, FALSE) & ifelse(!is.null(seasons), length(seasons) == 0, TRUE))){
    prior[, "cycle" := stats::predict(stats::smooth.spline(seasonalcycle))$y]
    prior[, "remainder" := y - trend - seasonal - cycle]
  }
  
  #Estimate the seasonal and cycle components
  if(ifelse(!is.null(seasons), length(seasons) > 0, FALSE)){
    ff = data.table(seasonalcycle = prior$seasonalcycle, t = 1:nrow(prior))
    for(j in c(seasons, cycle)){
      ff[, paste0("sin", j) := sin(2*pi*1/j*t)]
      ff[, paste0("cos", j) := cos(2*pi*1/j*t)]
    }
    ff[, "t" := NULL]
    lm_s = stats::lm(seasonalcycle ~ ., data = ff)
    prior = cbind(prior, do.call("cbind", lapply(seasons, function(x){
      matrix(stats::coef(lm_s)["(Intercept)"]/(length(seasons) + 1) +
               c(as.matrix(ff[, paste0(c("sin", "cos"), x), with = FALSE]) %*%
                   matrix(stats::coef(lm_s)[paste0(c("sin", "cos"), x)], ncol = 1)),
             ncol = 1)
    })))
    colnames(prior)[(ncol(prior) - length(seasons) + 1):ncol(prior)] = paste0("seasonal", seasons)
    prior[, "seasonal" := rowSums(prior[, paste0("seasonal", seasons), with = FALSE])]
    prior[, "cycle" := seasonalcycle - seasonal - lm_s$residuals]
    prior[, "remainder" := y - trend - seasonal - cycle]
  }
  
  prior[, "seasonal_adj" := trend + remainder]
  return(prior)
}
