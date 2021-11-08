#' Create dates to interpolate
#'
#' @param y Univariate time series of data values.
#' @param dates Vector of date values for y
#' @param interpolate Character string of how to interpolate
#' @param exo Matrix of exogenous variables. Can be used to specify regression effects or other seasonal effects like holidays, etc.
#' @import data.table
#' @return List of the data, dates, and exo
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
#' dates_interp = stsm_dates_to_interpolate(y = NA000334Q$y, dates = NA000334Q$date, 
#' interpolate = "monthly")
#' }
stsm_dates_to_interpolate = function(y, dates, exo = NULL, interpolate){
  `%m+%` = lubridate::`%m+%`
  data = data.table(date = dates, y)
  if(!is.null(exo)){
    data = cbind(data, exo)
  }
  
  if(interpolate == "quarterly"){
    new_dates = seq(min(dates), max(dates), by = "quarter")
  }else if(interpolate == "monthly"){
    new_dates = seq(min(dates), max(dates), by = "month")
  }else if(interpolate == "weekly"){
    new_dates = seq(min(dates), max(dates), by = "week")
  }else if(interpolate == "daily"){
    new_dates = seq(min(dates), max(dates), by = "day")
  }
  data = merge.data.table(data, data.table(date = new_dates), by = "date", all = T)
  return(list(y = data$y, dates = data$date, exo = data$exo))
}
