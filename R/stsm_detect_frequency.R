#' Detect frequency and dates from the data
#'
#' @param y Univariate time series of data values. May also be a 2 column data frame containing a date column.
#' @param freq Initial setting for the frequency detection
#' @import data.table
#' @return List giving the dates and frequency of the data
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
#' freq = stsm_detect_frequency(y = NA000334Q)
#' }
#' @export
stsm_detect_frequency = function(y, freq = NULL){
  if(stats::is.ts(y)){
    if(ifelse(is.null(ncol(y)), FALSE, ncol(y) > 1)){
      stop("Data must be a univariate time series.")
    }
    dates = stats::time(y)
    freq = stats::frequency(y)
  }else{
    y = as.data.table(y)
    datecol = unlist(y[, lapply(.SD, class), .SDcols = colnames(y)])
    datecol = names(datecol[datecol %in% c("Date", "yearmon")])
    if(length(datecol) == 1){
      if(class(y[, c(datecol), with = FALSE][[1]]) == "yearmon"){
        y[, c(datecol) := as.Date(eval(parse(text = datecol)))]
        is.yearmon = TRUE
      }else{
        is.yearmon = FALSE
      }
      if(ncol(y[, colnames(y) != datecol, with = FALSE]) > 1){
        stop("Data must be a univariate time series.")
      }
      if(is.yearmon == FALSE){
        y[, "day" := weekdays(eval(parse(text = datecol)))]
      }
      
      datediffs = unique(diff(unlist(unique(y[, c(datecol), with = FALSE]))))
      dates = unique(y[, c(datecol), with = FALSE][[1]])
      if(is.yearmon == TRUE){
        dates = unique(zoo::as.yearmon(dates))
      }
      if(is.null(freq)){
        freq = datediffs[which.max(tabulate(match(diff(y[, c(datecol), with = FALSE][[1]]), datediffs)))]
        freq = c(365.25, 365.25/7, 12, 4, 1)[which.min(abs(freq -  c(1, 365.25/7, 365.25/12, 365.25/4, 365.25)))]
        if(freq == 365.25 & all(!unique(y$day) %in% c("Saturday", "Sunday"))){
          freq = 365.25/7*5
        }
      }
      suppressWarnings(y[, "day" := NULL])
      y = unlist(y[, colnames(y)[colnames(y) != datecol], with = FALSE])
      rm(datediffs, datecol)
    }else if(length(datecol) > 1){
      stop("Too many date columns. Include only 1 date column or set the frequency manually.")
    }else if(length(datecol) == 0 & is.null(freq)){
      stop("No date column detected. Include a date column or set the frequency.")
    }else{
      y = unlist(y)
      dates = stats::time(stats::ts(y, freq = freq))
    }
  }
  return(list(data = y, dates = dates, freq = freq))
}
