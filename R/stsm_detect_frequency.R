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
  #Bind data.table variables to the global environment
  N = . = days = name = NULL
  
  #Define default value for standard_freq
  standard_freq = TRUE
  
  #Define a calendar to calculate days in a month, quarter, years
  calendar = data.table(month = c("January", "February", "March", "April", "May", "June", 
                                  "July", "August", "September", "October", "November", "December")) 
  calendar[, "mo" := 1:.N]
  calendar[, "qtr" := sort(rep(1:4, 3))]
  calendar[, "days" := c(31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)]
  
  ##### Test for the frequency of the data #####
  if(stats::is.ts(y)){
    #If the data is a ts object, it must by univariate
    if(ifelse(is.null(ncol(y)), FALSE, ncol(y) > 1)){
      stop("Data must be a univariate time series.")
    }
    dates = stats::time(y)
    freq = tryCatch(as.Date(as.numeric(stats::frequency(y))), 
                    error = function(err){
                      as.Date(as.numeric(stats::time(y)), origin = "0000-01-01")
                    })
  }else{
    #Find the date column
    y = as.data.table(y)
    datecol = unlist(y[, lapply(.SD, class), .SDcols = colnames(y)])
    datecol = names(datecol[datecol %in% c("Date", "yearmon", "POSIXct", "POSIXt", "POSIXlt")])
    datecol = sapply(colnames(y), function(x){any(grepl(x, datecol))})
    datecol = names(which(datecol))
    if(length(datecol) == 1){
      #Convert a yearmon column to a date column
      if(any(class(y[, c(datecol), with = FALSE][[1]]) == "yearmon")){
        y[, c(datecol) := as.Date(eval(parse(text = datecol)))]
      }
      
      #Make sure there is only one numeric column
      if(ncol(y[, colnames(y) != datecol, with = FALSE]) > 1){
        stop("Data must be a univariate time series.")
      }
      
      #Get the day of week names
      y[, "day" := weekdays(eval(parse(text = datecol)))]
      
      #Get the dates
      dates = unique(y[, c(datecol), with = FALSE][[1]])
      
      #Find the frequency if it is not already defined
      if(is.null(freq)){
        #Calculate the day differences between subsequent dates and get the unique values
        y[, "diff" := difftime(eval(parse(text = datecol)), 
                               shift(eval(parse(text = datecol)), type = "lag", n = 1), units = "days")]
        datediffs = unique(y$diff)
        
        #Get the most frequent date difference
        freq = as.numeric(sort(y[, .N, by = "diff"][N == max(N), ]$diff))[1]
        
        #Set the numeric frequency
        if(any(freq %in% c(sum(floor(calendar$days)), sum(ceiling(calendar$days))))){
          freq = 1
          name = "yearly"
        }else if(any(freq %in% unique(unlist(calendar[, .(sum(floor(days)), sum(ceiling(days))), by = "qtr"][, c("V1", "V2")])))){
          freq = 4
          name = "quarterly"
        }else if(any(freq %in% unique(c(floor(calendar$days), ceiling(calendar$days))))){
          freq = 12
          name = "monthly"
        }else if(any(freq %in% c(1, 3))){
          freq = 365.25
          name = "daily"
        }else if(any(freq %in% c(52, 53))){
          freq = 365.25/7
          name = "weekly"
        }else if(any(freq %in% c(1/c(23, 24)))){
          freq = 365.25*24
          name = "hourly"
        }else if(any(freq %in% c(1/c(c(23, 24)*60)))){
          freq = 365.25*24*60
          name = "minutely"
        }else if(any(freq %in% c(1/c(c(23, 24)*60*60)))){
          freq = 365.25*24*60*60
          name = "secondly"
        }else{
          freq = nrow(y)
          standard_freq = FALSE
          warning("Standard frequency (secondly, minutely, hourly, daily, weekly, monthly, quarterly, yearly), not detected. Defaulting to the length of the data set.")
        }
        
        #Test if the data should be weekday (no weekends) only
        if(freq >= 365.25 & standard_freq == TRUE & all(!unique(y$day) %in% c("Saturday", "Sunday")) & length(unique(y$day)) >= 5){
          freq = freq*5/7
          name = paste0(name, ", weekday only")
        }
        y[, "diff" := NULL]
      }
      y[, "day" := NULL]
      y = unlist(y[, colnames(y)[colnames(y) != datecol], with = FALSE])
    }else if(length(datecol) > 1){
      stop("Too many date columns. Include only 1 date column or set the frequency manually.")
    }else if(length(datecol) == 0 & is.null(freq)){
      stop("No date column detected. Include a date column or set the frequency.")
    }
  }
  return(list(data = y, dates = dates, freq = freq, standard_freq = standard_freq, name = name))
}