#' Get initial values for the Kalman filter
#'
#' @param y an object created from stsm_detect_frequency
#' @param par parameter values for the state space model
#' @param freq Frequency of the data
#' @param trend Trend specification ("random-walk", "random-walk-drift", "double-random-walk", "random-walk2"). 
#' @param decomp Decomposition model ("tend-cycle-seasonal", "trend-seasonal", "trend-cycle", "trend-noise")
#' @param seasons The seasonal periods to split the seasonality into
#' @param prior A data table created by stsm_prior
#' @param cycle The cycle period
#' @import data.table
#' @return list containing the initial values for the Kalman filter
stsm_init_vals = function(y, par, freq, trend, decomp = "", seasons = NULL, prior = NULL, cycle = NULL){
  #Bind data.table variables to the global environment
  drift = remainder = NULL
  
  #Build the prior
  if(is.null(prior)){
    prior = stsm_prior(y, freq, decomp, seasons, cycle) 
  }else{
    prior = copy(prior)
  }
  
  sp = stsm_ssm(par, y, decomp, trend)
  init = list(B0 = sp$B0, P0 = sp$P0)
  if("trend" %in% colnames(prior)){
    init[["B0"]][names(init[["B0"]]) == "Tt_0"] = prior[!is.na(trend), ]$trend[1]
  }
  if("cycle" %in% colnames(prior)){
    init[["B0"]][names(init[["B0"]]) %in% c("Ct_0", "Cts_0")] = prior[!is.na(cycle), ]$cycle[1]
  }
  if("drift" %in% colnames(prior)){
    init[["B0"]][names(init[["B0"]]) == "Dt_0"] = prior[!is.na(drift), ]$drift[1]
  }
  if("remainder" %in% colnames(prior)){
    init[["B0"]][names(init[["B0"]]) == "et_0"] = prior[!is.na(remainder), ]$remainder[1]
  }
  if(grepl("seasonal", decomp)){
    if(!is.null(seasons)){
      for(j in seasons){
        init[["B0"]][names(init[["B0"]]) %in% paste0(c("St", "Sts"), j, "_0")] = 
          prior[!is.na(eval(parse(text = paste0("seasonal", j)))), c(paste0("seasonal", j)), with = FALSE][[1]][1]
      }
    }
  }
  return(init)
}
