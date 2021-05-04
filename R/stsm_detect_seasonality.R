#' Detect seasonality from the data
#'
#' @param y Univariate time series of data values.
#' @param freq Frequency of the data (1 (yearly), 4 (quarterly), 12 (monthly), 365.25/7 (weekly), 365.25 (daily))
#' @param sig_level Significance level to determine statistically significant seasonal frequencies
#' @param prior A data table created from stsm_prior
#' @param cores Number of cores to use
#' @param cl a parallel cluster object
#' @param show_progress Whether to show progress bar
#' @import data.table
#' @return Numeric vector of seasonal periodicities
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
#' seasonality = stsm_detect_seasonality(y = NA000334Q$y, freq = 4)
#' }
#' @export
stsm_detect_seasonality = function(y, freq, sig_level = 0.01, prior = NULL, 
                                   cl = NULL, cores = NULL, show_progress = FALSE){
  #Bind data.table variables to the global environment
  s = seasons = round = period = period2 = test = frequency = frequency2 = 
    trend = rn = df = pval = NULL
  
  #Return no seasons if the frequency is yearly
  if(freq == 1){
    return(seasons = numeric(0))
  }
  
  #Get the prior
  if(is.null(prior)){
    prior = stsm_prior(y, freq) 
  }else{
    prior = copy(prior)
  }
  
  #Adjust the seasonal value
  prior[, "seasonal" := y - trend]
  if(tsutils::coxstuart(stats::na.omit(prior$seasonal), type = "dispersion")$p.value <= sig_level){
    prior[, "seasonal" := y/trend]
  }
  prior[, "t" := 1:.N]
  
  #Wavelet analysis for seasonality using forward stepwise regression
  if(floor(freq) == floor(60*60*24*365.25)){
    #Secondly frequency
    #minutely, hourly, workday hours, 1/2 daily, daily, weekend, weekday, weekly, monthly, quarterly, semi-yearly, yearly
    spectrum = rev(c(60, 60*60, 60*60*8, 60*60*12, 60*60*24, 60*60*24*2, 60*60*24*5, 60*60*24*7, freq/12, freq/4, freq/2, freq))
  }else if(floor(freq) == floor(60*60*24*365.25*5/7)){
    #Secondly frequency, weekday only
    #minutely, hourly, workday hours, 1/2 daily, weekday, monthly, quarterly, semi-yearly, yearly
    spectrum = rev(c(60, 60*60, 60*60*8, 60*60*12, 60*60*24, 60*60*24*5, freq/12, freq/4, freq/2, freq))
  }else if(floor(freq) == floor(60*24*365.25)){
    #Minutely frequency
    #hourly, workday hours, 1/2 daily, daily, weekend, weekday, weekly, monthly, quarterly, semi-yearly, yearly
    spectrum = rev(c(60, 60*8, 60*12, 60*24, 60*24*2, 60*24*5, 60*24*7, freq/12, freq/4, freq/2, freq))
  }else if(floor(freq) == floor(60*24*365.25*5/7)){
    #Minutely frequency, weekday only
    #hourly, workday hours, 1/2 daily, daily, weekday, monthly, quarterly, semi-yearly, yearly
    spectrum = rev(c(60, 60*8, 60*12, 60*24, 60*24*5, freq/12, freq/4, freq/2, freq))
  }else if(floor(freq) == floor(24*365.25)){
    #Hourly frequency
    #workday hours, 1/2 daily, daily, weekend, weekday, weekly, monthly, quarterly, semi-yearly, yearly
    spectrum = rev(c(8, 12, 24, 24*2, 24*5, 24*7, freq/12, freq/4, freq/2, freq))
  }else if(floor(freq) == floor(24*365.25*5/7)){
    #Hourly frequency, weekday only
    #workday hours, 1/2 daily, daily, weekday, monthly, quarterly, semi-yearly, yearly
    spectrum = rev(c(8, 12, 24, 24*5, freq/12, freq/4, freq/2, freq))
  }else if(floor(freq) == 365){
    #Daily frequency
    #weekend, weekday, weekly, monthly, quarterly, semi-yearly, yearly
    spectrum = rev(c(2, 5, 7, freq/12, freq/4, freq/2, freq))
  }else if(floor(freq) == floor(365.25*5/7)){
    #Daily frequency, weekday only
    #weekly, monthly, quarterly, semi-yearly, yearly
    spectrum = rev(c(5, freq/12, freq/4, freq/2, freq))
  }else if(floor(freq) == 52){
    #Weekly frequency
    #monthly, quarterly, semi-yearly
    spectrum = rev(c(freq/12, freq/4, freq/2, freq))
  }else if(floor(freq) == 12){
    #Monthly frequency
    #quarterly, semi-yearly, yearly
    spectrum = rev(c(freq/4, freq/2, freq))
  }else if(floor(freq) == 4){
    #Quarterly frequency
    #semi-yearly, yearly
    spectrum = rev(c(freq/2, freq))
  }else{
    spectrum = NULL
  }
  
  if(!is.null(spectrum)){
    #Define the harmonics for the regressions
    harmonics = freq/spectrum
    #Create a linear space of 100 harmonics between those defined by the spectrum of interest
    # so that each harmonic represents 1% of the space between the spectrum
    harmonics = sort(unique(c(harmonics, unlist(lapply(1:(length(harmonics) - 1), function(x){
      s = seq(harmonics[x], harmonics[x + 1], length.out = 102)
      unique(round(s[!s %in% harmonics]))
    })))))
    #Filter out harmonics that are within 1 of the spectrum of interest
    harmonics = harmonics[!(harmonics %in% floor(freq/spectrum) |
                           harmonics %in% ceiling(freq/spectrum)) |
                            harmonics %in% c(freq/spectrum)]
  }else{
    harmonics = 1:floor(freq/2)
    #Check a maximum of 1000 harmonics
    harmonics = harmonics[seq(1, length(harmonics), length.out = 1000)]
    #Limit the harmonics to only those that have 3 full cycles in the data
    harmonics = harmonics[freq/harmonics <= nrow(prior)/3]
  }
  
  if(length(harmonics) > 0){
    #Setup parallel computing
    if(is.null(cl)){
      cl = parallel::makeCluster(max(c(1, ifelse(is.null(cores), min(c(length(harmonics), parallel::detectCores())), cores))))
      doSNOW::registerDoSNOW(cl)
      stop_cluster = TRUE
    }else{
      stop_cluster = FALSE
    }
    `%fun%` = foreach::`%dopar%`
    
    if(show_progress == TRUE){
      pb = progress::progress_bar$new(
        format = " [:bar] :percent complete |:elapsed elapsed |:eta remaining ",
        total = length(harmonics),  
        clear = FALSE, show_after = 0
      )
      invisible(pb$tick(0))
      progress = function(n){pb$tick()}
    }else{
      progress = NULL
    }
    
    #Wavelet analysis for the seasonality
    wave = foreach::foreach(i = harmonics, .combine = "rbind",
                            .packages = c("data.table", "stats"),
                            .options.snow = list(progress = progress)) %fun% {
      
      #Build the ith harmonic
      prior = copy(prior)                       
      prior[, colnames(prior)[grepl("sin|cos", colnames(prior))] := NULL]
      prior[, paste0("sin", i) := sin(2*pi*i/freq*t)]
      prior[, paste0("cos", i) := cos(2*pi*i/freq*t)]
      
      #Perform the harmonic regression
      lm = stats::lm(seasonal ~ ., prior[, c("seasonal", colnames(prior)[grepl("sin|cos", colnames(prior))]), with = FALSE])
      pval = stats::pf(summary(lm)$fstatistic["value"], lower.tail = FALSE,
                       df1 = summary(lm)$fstatistic["numdf"], 
                       df2 = summary(lm)$fstatistic["dendf"])
      return(data.table(frequency = i/freq, period = freq/i, fstat = summary(lm)$fstatistic["value"], pval = pval))
    }
    wave[, "df" := unique(signif(diff(frequency)))[1]]
    seasons = wave[pval <= sig_level, ]
    
    if(nrow(seasons) > 0){
      #Assign the closest value in the predefined spectrum
      if(!is.null(spectrum)){
        for(i in 1:nrow(seasons)){
          val = seasons$period[i]
          seasons[i, "period2" := spectrum[which.min(abs(val - spectrum))]]
        }
        seasons[, "frequency2" := 1/period2]
        
        #Set to t to the predefined spectrum t if its frequency falls within one
        #difference of the predefined spectrum's frequency
        seasons[, "test" := frequency >= frequency2 - df & frequency <= frequency2 + df]
        seasons[, "period" := ifelse(test == TRUE, period2, period)]
        seasons[test == FALSE, "period" := round(period)]
      }else{
        seasons[, "period" := round(period)]
      }
      seasons = unique(seasons[frequency <= freq, ]$period)
    }else{
      seasons = numeric(0)
    }
    
    #Final specification tests
    if(length(seasons) > 0){
      for(i in seasons){
        prior[, paste0("sin", i) := sin(2*pi*1/i*t)]
        prior[, paste0("cos", i) := cos(2*pi*1/i*t)]
      }
      
      #Backward stepwise regression using F-tests for paired coefficients
      #Only take seasonalities that are statistically significant
      go = TRUE
      while(go == TRUE){
        #Estimate the final harmonic regression
        lm = stats::lm(seasonal ~ ., prior[, c("seasonal", colnames(prior)[grepl("sin|cos", colnames(prior))]), with = FALSE])
        coef = summary(lm)$coefficients
        coef = data.table(coef, keep.rownames = TRUE)
        coef[, "group" := gsub("sin|cos", "", rn)]
        
        #Calculate F-test for joint significance of each seasonality (the sin and cos coefficients)
        f_tests = foreach::foreach(i = unique(coef[rn != "(Intercept)", ]$group), .combine = "rbind",
                                   .packages = c("data.table", "sandwich", "lmtest")) %fun% {
          
          #Estimate the restricted regression without seasonality i
          dat = prior[, c("seasonal", colnames(prior)[grepl("sin|cos", colnames(prior))]), with = FALSE]
          dat[, paste0(c("sin", "cos"), i) := NULL]
          lm_res = stats::lm(seasonal ~ ., dat)    
          
          #Calculate robust F-test for joint significance of each seasonality
          test = tryCatch(lmtest::waldtest(lm_res, lm, vcov = sandwich::vcovHAC(lm)), 
                          error = function(x){
                            list(`Pr(>F)` = list(NA, stats::pf(summary(lm)$fstatistic["value"], lower.tail = FALSE,
                                                               df1 = summary(lm)$fstatistic["numdf"], 
                                                               df2 = summary(lm)$fstatistic["dendf"])))
                          })
          return(data.table(season = i, pval = test$`Pr(>F)`[2]))
        }
        
        #Remove insignificant seasonalities
        if(nrow(f_tests[pval > sig_level, ]) > 0){
          cols = unlist(lapply(f_tests[pval > sig_level, ]$season, function(x){paste0(c("sin", "cos"), x)}))
          prior[, c(cols) := NULL]
        }
        
        #If no seasonalities are significant or all seasonalities are significant, stop the loop
        if(length(colnames(prior)[grepl("sin|cos", colnames(prior))]) == 0 |
           nrow(f_tests[pval > sig_level, ]) == 0){
          go = FALSE
        }
      }
      
      if(nrow(f_tests[pval <= sig_level, ]) > 0){    
        #Calculate robust F-test for the final specification
        lm = stats::lm(seasonal ~ ., prior[, c("seasonal", colnames(prior)[grepl("sin|cos", colnames(prior))]), with = FALSE])
        lm_res = stats::lm(seasonal ~ 1, prior)
        final_ftest = tryCatch(lmtest::waldtest(lm_res, lm, vcov = sandwich::vcovHAC(lm)), 
                        error = function(x){
                          list(`Pr(>F)` = list(NA, stats::pf(summary(lm)$fstatistic["value"], lower.tail = FALSE,
                                                    df1 = summary(lm)$fstatistic["numdf"], 
                                                    df2 = summary(lm)$fstatistic["dendf"])))
                        })
        
        if(final_ftest$`Pr(>F)`[2] <= sig_level){
          seasons = unique(as.numeric(f_tests[pval <= sig_level, ]$season))
        }else{
          seasons = numeric(0)
        }
      }else{
        seasons = numeric(0)
      }
    }else{
      seasons = numeric(0)
    }
  
    if(stop_cluster == TRUE){
      parallel::stopCluster(cl)
    }
  }else{
    seasons = numeric(0)
  }
  
  return(seasons)
}
