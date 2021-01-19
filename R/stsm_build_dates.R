#' Build the date sequence as a Date type
#'
#' @param y a list object created from stsm_detect_frequency
#' @import data.table
#' @return a list with the univariate time series and corrected dates
stsm_build_dates = function(y){
  `%m+%` = lubridate::`%m+%`
  dates = tryCatch(as.Date(y$dates),
                   error = function(err){
                     diff = mean(unique(diff(y$dates)))
                     years = as.numeric(floor(y$dates))
                     parts = as.numeric(y$dates - years)
                     if(floor(y$freq) == 1){
                       dates = as.Date(paste0(years, "-01-01")) %m+% lubridate::years(ceiling((parts/diff)))
                     }else if(floor(y$freq) == 4){
                       dates = as.Date(paste0(years, "-01-01")) %m+% months(ceiling((parts/diff))*3)
                     }else if(floor(y$freq) == 12){
                       dates = as.Date(paste0(years, "-01-01")) %m+% months(ceiling((parts/diff)))
                     }else if(floor(y$freq) == 52){
                       dates = as.Date(paste0(years, "-01-01")) %m+% lubridate::weeks(ceiling((parts/diff)))
                     }else if(floor(y$freq) == 365){
                       dates = as.Date(paste0(years, "-01-01")) %m+% lubridate::days(ceiling((parts/diff)))
                     }
                     return(dates)
                   })
  
  #Make sure data is not skipping dates
  if(floor(y$freq) == 365){
    dates2 = seq(lubridate::ymd(min(dates)), lubridate::ymd(max(dates)), by = "day")
  }else if(floor(y$freq) == 260){
    dates2 = seq(lubridate::ymd(min(dates)), lubridate::ymd(max(dates)), by = "day")
    dates2 = dates2[which(!weekdays(dates2) %in% c("Saturday", "Sunday"))]
  }else if(floor(y$freq) == 52){
    dates2 = seq(lubridate::ymd(min(dates)), lubridate::ymd(max(dates)), by = "week")
  }else if(floor(y$freq) == 12){
    dates2 = seq(lubridate::ymd(min(dates)), lubridate::ymd(max(dates)), by = "month")
  }else if(floor(y$freq) == 4){
    dates2 = seq(lubridate::ymd(min(dates)), lubridate::ymd(max(dates)), by = "quarter")
  }else if(floor(y$freq) == 1){
    dates2 = seq(lubridate::ymd(min(dates)), lubridate::ymd(max(dates)), by = "year")
  }
  
  if(length(dates) != length(dates2)){
    dt = data.table(date = dates, y = y$data)
    date_dt = data.table(date = dates2)
    dt = merge(date_dt, dt, all = TRUE, by = "date")
    if(floor(y$freq) == 260){
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