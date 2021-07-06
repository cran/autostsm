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
  
  #Impute missing values for using the Kalman filter for prior decomposition
  if(any(is.na(y))){
    y = suppressWarnings(stsm_na_kalman(y))
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
    temp = forecast::mstl(forecast::msts(y - prior$trend, seasonal.periods = seasons, ts.frequency = freq))
    colnames(temp)[grepl("Seasonal", colnames(temp))] = sapply(colnames(temp)[grepl("Seasonal", colnames(temp))], function(x){
      s = as.numeric(gsub("Seasonal", "", x))
      return(paste0("seasonal", seasons[which.min(abs(s - seasons))]))
    })
    prior = cbind(prior, temp[, grepl("seasonal", colnames(temp))])
    prior[, colnames(prior)[grepl("seasonal", colnames(prior))] := lapply(.SD, as.numeric), 
          .SDcols = colnames(prior)[grepl("seasonal", colnames(prior))]]
    prior[, "cycle" := as.numeric(temp[, "Trend"])]
    prior[, "seasonalcycle" := rowSums(temp[, !colnames(temp) %in% c("Data")])]
    prior[, "seasonal" := rowSums(prior[, paste0("seasonal", seasons), with = FALSE])]
    prior[, "remainder" := y - trend - seasonal - cycle]
  }
  suppressWarnings(prior[, "seasonalcycle" := NULL])
  return(prior)
}

