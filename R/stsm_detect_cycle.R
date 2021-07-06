#' Detect cycle from the data
#'
#' @param y Univariate time series of data values.
#' @param freq Frequency of the data (1 (yearly), 4 (quarterly), 12 (monthly), 365.25/7 (weekly), 365.25 (daily))
#' @param sig_level Significance level to determine statistically significant seasonal frequencies
#' @param cores Number of cores to use
#' @param show_progress Whether to show progress bar
#' @param prior A data table created by stsm_prior
#' @param cl a parallel cluster object
#' @import data.table
#' @return Numeric value of cycle periodicity
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
#' cycle = stsm_detect_cycle(y = NA000334Q$y, freq = 4)
#' }
#' @export
stsm_detect_cycle = function(y, freq, sig_level = 0.01, prior = NULL, 
                             cl = NULL, cores = NULL, show_progress = FALSE){
  #Bind data.table variables to the global environment
  seasonal = trend = fstat = pval = NULL
  
  #Get the prior
  if(is.null(prior)){
    prior = stsm_prior(y, freq) 
  }else{
    prior = copy(prior)
  }
  
  #Redefine the cycle value
  prior[, "cycle" := y - trend]
  if(stsm_coxstuart(stats::na.omit(forecast::tsclean(prior$cycle, replace.missing = FALSE)), 
                    type = "deviation")$p.value <= sig_level){
    prior[, "cycle" := y/trend]
  }
  prior[, "t" := 1:.N]
  
  #Setup parallel computing
  harmonics = seq(0.01, 0.99, by = 0.01)
  harmonics = harmonics[freq/harmonics >= 2.5*freq & harmonics < 1/2]
  harmonics = harmonics[freq/harmonics <= nrow(prior)]
  if(length(harmonics) > 0){
    if(is.null(cl) & ifelse(!is.null(cores), cores > 1, TRUE)){
      cl = tryCatch(parallel::makeCluster(max(c(1, ifelse(is.null(cores), min(c(length(harmonics), parallel::detectCores())), cores)))),
                    error = function(err){
                      message("Parallel setup failed. Using single core.")
                      return(NULL)
                    })
      if(!is.null(cl)){
        doSNOW::registerDoSNOW(cl)
        `%fun%` = foreach::`%dopar%`
        stop_cluster = TRUE
      }else{
        stop_cluster = FALSE
        `%fun%` = foreach::`%do%`
      }
    }else{
      stop_cluster = FALSE
      if(length(cl) > 1){
        `%fun%` = foreach::`%dopar%`
      }else{
        `%fun%` = foreach::`%do%`
      }
    }
    
    if(show_progress == TRUE){
      pb = progress::progress_bar$new(
        format = " [:bar] :percent complete |:elapsed elapsed |:eta remaining ",
        total = length(harmonics),  
        clear = FALSE, show_after = 0
      )
      invisible(pb$tick(0))
      progress = function(n){pb$tick()}
    }else{
      pb = progress = NULL
    }
    
    #Wavelet analysis for the cycle
    wave = foreach::foreach(i = harmonics, .combine = "rbind", 
                            .packages = c("data.table", "stats"),
                            .options.snow = list(progress = progress)) %fun% {
      #Build the ith harmonic
      prior = copy(prior)                       
      suppressWarnings(prior[, colnames(prior)[grepl("sin|cos", colnames(prior))] := NULL])
      prior[, paste0("sin", i) := sin(2*pi*i/freq*t)]
      prior[, paste0("cos", i) := cos(2*pi*i/freq*t)]
      
      #Perform the harmonic regression
      lm = stats::lm(cycle ~ ., prior[, c("cycle", colnames(prior)[grepl("sin|cos", colnames(prior))]), with = FALSE])
      pval = stats::pf(summary(lm)$fstatistic["value"], lower.tail = FALSE,
                       df1 = summary(lm)$fstatistic["numdf"], 
                       df2 = summary(lm)$fstatistic["dendf"])
      
      if(stop_cluster == FALSE & show_progress == TRUE){
        pb$tick()
      }
      return(data.table(frequency = i/freq, period = freq/i, fstat = summary(lm)$fstatistic["value"], pval = pval))
    }
    
    #Find the maxima
    cycle = wave[which.max(fstat), ][pval <= sig_level, ]$period
    
    #Check the final harmonic specification
    if(length(cycle) > 0){
      i = cycle
      prior[, paste0("sin", i) := sin(2*pi*1/i*t)]
      prior[, paste0("cos", i) := cos(2*pi*1/i*t)]
      
      #Estimate the final harmonic regression and the restricted regression with a constant only
      lm_res = stats::lm(cycle ~ 1, prior)
      lm = stats::lm(cycle ~ ., prior[, c("cycle", colnames(prior)[grepl("sin|cos", colnames(prior))]), with = FALSE])
      
      #Calculate robust F-test for joint significance of the cycle (the sin and cos coefficients)
      final_ftest = tryCatch(lmtest::waldtest(lm_res, lm, vcov = sandwich::vcovHAC(lm)), 
                      error = function(x){
                        list(`Pr(>F)` = list(NA, stats::pf(summary(lm)$fstatistic["value"], lower.tail = FALSE,
                                                           df1 = summary(lm)$fstatistic["numdf"], 
                                                           df2 = summary(lm)$fstatistic["dendf"])))
                      })
     if(final_ftest$`Pr(>F)`[2] > sig_level){
        cycle = numeric(0)
      }
    }else{
      cycle = numeric(0)
    }
    
    if(stop_cluster == TRUE){
      parallel::stopCluster(cl)
    }
  }else{
    cycle = numeric(0)
  }
  return(cycle)
}

