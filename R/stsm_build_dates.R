#' Build the date sequence as a Date type
#'
#' @param y a list object created from stsm_detect_frequency
#' @import data.table
#' @return a list with the univariate time series and corrected dates
stsm_build_dates = function(y){
  `%m+%` = lubridate::`%m+%`
  
  #Get the dates and handle errors
  dates = y$dates
  
  #Make sure data is not skipping dates
  if(y$standard_freq == TRUE){
    if(floor(y$freq) == floor(365.25*24*60*60)){
      #secondly
      dates2 = seq(lubridate::ymd_hms(min(dates)), lubridate::ymd_hms(max(dates)), by = "sec")
    }else if(floor(y$freq) == floor(60*60*24*365.25*5/7)){
      #secondly, weekdays only
      dates2 = seq(lubridate::ymd_hms(min(dates)), lubridate::ymd_hms(max(dates)), by = "sec")
      dates2 = dates2[which(!weekdays(dates2) %in% c("Saturday", "Sunday"))]
    }else if(floor(y$freq) == floor(60*24*365.25)){
      #minutely
      dates2 = seq(lubridate::ymd_hms(min(dates)), lubridate::ymd_hms(max(dates)), by = "min")
    }else if(floor(y$freq) == floor(60*24*365.25*5/7)){
      #minutely, weekdays only
      dates2 = seq(lubridate::ymd_hms(min(dates)), lubridate::ymd_hms(max(dates)), by = "min")
      dates2 = dates2[which(!weekdays(dates2) %in% c("Saturday", "Sunday"))]
    }else if(floor(y$freq) == floor(24*365.25)){
      #hourly
      dates2 = seq(lubridate::ymd_hms(min(dates)), lubridate::ymd_hms(max(dates)), by = "hour")
    }else if(floor(y$freq) == floor(24*365.25*5/7)){
      #hourly
      dates2 = seq(lubridate::ymd_hms(min(dates)), lubridate::ymd_hms(max(dates)), by = "hour")
      dates2 = dates2[which(!weekdays(dates2) %in% c("Saturday", "Sunday"))]
    }else if(floor(y$freq) == 365){
      #daily
      dates2 = seq(lubridate::ymd(min(dates)), lubridate::ymd(max(dates)), by = "day")
    }else if(floor(y$freq) == floor(365*5/7)){
      #daily, weekdays only
      dates2 = seq(lubridate::ymd(min(dates)), lubridate::ymd(max(dates)), by = "day")
      dates2 = dates2[which(!weekdays(dates2) %in% c("Saturday", "Sunday"))]
    }else if(floor(y$freq) == 52){
      #weekly
      dates2 = seq(lubridate::ymd(min(dates)), lubridate::ymd(max(dates)), by = "week")
    }else if(floor(y$freq) == 12){
      #monthly
      dates2 = seq(lubridate::ymd(min(dates)), lubridate::ymd(max(dates)), by = "month")
    }else if(floor(y$freq) == 4){
      #quarterly
      dates2 = seq(lubridate::ymd(min(dates)), lubridate::ymd(max(dates)), by = "quarter")
    }else if(floor(y$freq) == 1){
      #yearly
      dates2 = seq(lubridate::ymd(min(dates)), lubridate::ymd(max(dates)), by = "year")
    }else{
      dates2 = dates
    }
  }else{
    dates2 = dates
  }
  
  #Combine date sequences create missing values
  if(length(dates) != length(dates2)){
    dt = data.table(date = dates, y = y$data)
    date_dt = data.table(date = dates2)
    dt = merge(date_dt, dt, all = TRUE, by = "date")
    if(floor(y$freq) %in% floor(c(60*60*24, 60*24, 24, 1)*365.25*5/7)){
      dt[, "day" := weekdays(date)]
      dt = dt[!date %in% c("Saturday", "Sunday"), ]
      dt[, "day" := NULL]
    }
    y$data = dt$y
    y$dates = dt$date
    rm(dt, date_dt)
  }
  
  return(y)
}