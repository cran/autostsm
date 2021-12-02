#' Detect trend type
#'
#' @param y Univariate time series of data values. May also be a 2 column data frame containing a date column.
#' @param freq Frequency of the data (1 (yearly), 4 (quarterly), 12 (monthly), 365.25/7 (weekly), 365.25 (daily))
#' @param decomp Decomposition model ("tend-cycle-seasonal", "trend-seasonal", "trend-cycle", "trend-noise")
#' @param sig_level Significance level to determine statistically significant seasonal frequencies
#' @param prior A data table created by stsm_prior
#' @param seasons The seasonal periods
#' @param cycle The cycle period
#' @param cores Number of cores to use
#' @param cl a parallel cluster object
#' @param verbose Logical whether to print messages or not
#' @import data.table
#' @return list with trend type and logical flag for deterministic trend if the trend is determined to have 0 differencing
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
#' trend = stsm_detect_trend(y = NA000334Q$y, freq = 4)
#' }
#' @export
stsm_detect_trend = function(y, freq, decomp = "", sig_level = 0.01, prior = NULL, seasons = NULL, cycle = NULL, 
                             cl = NULL, cores = NULL, verbose = FALSE){
  #Bind data.table variables to the global environment
  N = d = test = pval = sig = trend = remainder = seasonal = seasonal_adj = 
    drift = m = m.lead = NULL
  
  #Set the prior
  if(is.null(prior)){
    prior = stsm_prior(y, freq, decomp , seasons, cycle) 
  }else{
    prior = copy(prior)
  }
  prior[, "seasonal_adj" := y - seasonal - cycle]
  
  ##### Find the number of differences to make it stationary #####
  #Replace outliers for the trend test
  prior[, "seasonal_adj" := forecast::tsclean(seasonal_adj, replace.missing = FALSE)]
  #Test for a trend
  trend_test = (round(stsm_coxstuart(stats::na.omit(y - prior$seasonal), type = "trend")$p.value, 2) <= sig_level)
  #Calculate unit root tests for differences 0 to 2
  maxlags = max(c(trunc(12*(length(stats::na.omit(prior$seasonal_adj))/100)^(1/4)), #Schwert method
                  trunc((length(stats::na.omit(prior$seasonal_adj)) - 1)^(1/3)))) #arma method
  if(trend_test == FALSE){
    ndiffs = c(forecast::ndiffs(stats::na.omit(prior$seasonal_adj), test = "adf", type = "level", alpha = sig_level, max.d = 2, lags = maxlags, selectlags = "BIC"), 
             forecast::ndiffs(stats::na.omit(prior$seasonal_adj), test = "pp", type = "level", alpha = sig_level, max.d = 2, lags = "long"), 
             forecast::ndiffs(stats::na.omit(prior$seasonal_adj), test = "kpss", type = "level", alpha = sig_level, max.d = 2, lags = "long"))
    ndiffs = round(mean(ndiffs, na.rm = TRUE))
  }else{
    ndiffs = c(forecast::ndiffs(diff(stats::na.omit(prior$seasonal_adj)), test = "adf", type = "level", alpha = sig_level, max.d = 1, lags = maxlags, selectlags = "BIC"), 
                forecast::ndiffs(diff(stats::na.omit(prior$seasonal_adj)), test = "pp", type = "level", alpha = sig_level, max.d = 1, lags = "long"), 
                forecast::ndiffs(diff(stats::na.omit(prior$seasonal_adj)), test = "kpss", type = "level", alpha = sig_level, max.d = 1, lags = "long"))
    ndiffs = round(mean(ndiffs, na.rm = TRUE)) + 1
  }
  
  #Check for a structural break in the drift going from positive to negative
  if(ndiffs == 1 & min(prior$drift, na.rm = TRUE) < 0){
    if(verbose == TRUE){
      message("Checking for structural breaks in the drift...")
    }
    
    #Setup parallel computing
    if(is.null(cl) & ifelse(!is.null(cores), cores > 1, TRUE)){
      #if cluster is null and (cores is not null or > 0)
      cl = tryCatch(parallel::makeCluster(max(c(1, ifelse(is.null(cores), parallel::detectCores())))),
                    error = function(err){
                      message("Parallel setup failed. Using single core.")
                      return(NULL)
                    })
      hpc = ifelse(!is.null(cl), "foreach", "none")
      stop_cluster = TRUE
    }else{
      stop_cluster = FALSE
      hpc = ifelse(length(cl) > 1, "foreach", "none")
    }
    
    bp = suppressWarnings(strucchange::breakpoints(stats::na.omit(prior$drift) ~ 1, hpc = hpc))
    if(!all(is.na(bp$breakpoints))){
      breaks = data.table(b = bp$breakpoints, m = as.numeric(NA), m.lead = as.numeric(NA))
      if(length(bp$breakpoints) == 1){
        breaks$m = mean(prior[!is.na(drift), ][1:breaks$b[1], ]$drift, na.rm = TRUE)
        breaks$m.lead = mean(prior[!is.na(drift), ][(breaks$b + 1):.N, ]$drift, na.rm = TRUE)
      }else{
        for(i in 1:nrow(breaks)){
          if(breaks$b[i] == min(breaks$b)){
            breaks[i, ]$m = mean(prior[!is.na(drift), ][1:breaks$b[i], ]$drift, na.rm = TRUE)
            breaks[i, ]$m.lead = mean(prior[!is.na(drift), ][(breaks$b[i] + 1):breaks$b[i + 1], ]$drift, na.rm = TRUE)
          }else if(breaks$b[i] == max(breaks$b)){
            breaks[i, ]$m = mean(prior[!is.na(drift), ][breaks$b[i]:.N, ]$drift, na.rm = TRUE)
            breaks[i, ]$m.lead = mean(prior[!is.na(drift), ][(breaks$b[i] + 1):.N, ]$drift, na.rm = TRUE)
          }else{
            breaks[i, ]$m = mean(prior[!is.na(drift), ][breaks$b[i-1]:breaks$b[i], ]$drift, na.rm = TRUE)
            breaks[i, ]$m.lead = mean(prior[!is.na(drift), ][(breaks$b[i] + 1):.N, ]$drift, na.rm = TRUE)
          }
        }
      }
      
      breaks[, "test" := (sign(m) == 1 & sign(m.lead) == -1) | (sign(m) == -1 & sign(m.lead) == 1)]
      if(any(breaks$test == TRUE)){
        ndiffs = 2
      }
    }
    if(stop_cluster == TRUE){
      parallel::stopCluster(cl)
    }
  }
  
  if(ndiffs >= 2){
    trend = "double-random-walk"
  }else if(ndiffs == 1 & trend_test == TRUE){
    trend = "random-walk-drift"
  }else if(ndiffs == 1 & trend_test == FALSE){
    trend = "random-walk"
  }else if(ndiffs == 0 & trend_test == FALSE){
    trend = "random-walk"
  }else if(ndiffs == 0 & trend_test == TRUE){
    trend = "random-walk-drift"
  }
  
  return(list(trend = trend, det_trend = (ndiffs == 0)))
}
