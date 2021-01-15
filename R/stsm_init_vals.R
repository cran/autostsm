#' Get initial values for the Kalman filter
#'
#' @param y an object created from stsm_detect_frequency
#' @param par parameter values for the state space model
#' @param freq Frequency of the data
#' @param trend Trend specification ("random-walk", "random-walk-drift", "double-random-walk", "random-walk2"). 
#' @param decomp Decomposition model ("tend-cycle-seasonal", "trend-seasonal", "trend-cycle", "trend-noise")
#' @param harmonics The harmonics to split the seasonality into
#' @param prior A data table created by stsm_prior
#' @import data.table
#' @return list containing the initial values for the Kalman filter
stsm_init_vals = function(y, par, freq, trend, decomp = "", harmonics = NULL, prior = NULL){
  if(is.null(prior)){
    prior = stsm_prior(y, freq, decomp, harmonics) 
  }else{
    prior = copy(prior)
  }
  
  cycle = drift = remainder = NULL
  sp = stsm_ssm(par = par, yt = y, freq = freq, decomp = decomp, trend = trend)
  init = list(B0 = sp$B0, P0 = sp$P0)
  if("trend" %in% colnames(prior)){
    init[["B0"]][names(init[["B0"]]) == "Tt0"] = prior[!is.na(trend), ]$trend[1]
  }
  if("cycle" %in% colnames(prior)){
    init[["B0"]][names(init[["B0"]]) %in% c("Ct", "Cts")] = prior[!is.na(cycle), ]$cycle[1]
  }
  if("drift" %in% colnames(prior)){
    init[["B0"]][names(init[["B0"]]) == "Mt0"] = prior[!is.na(drift), ]$drift[1]
  }
  if("remainder" %in% colnames(prior)){
    init[["B0"]][names(init[["B0"]]) == "et0"] = prior[!is.na(remainder), ]$remainder[1]
  }
  if(grepl("seasonal", decomp)){
    if(!is.null(harmonics)){
      for(j in harmonics){
        init[["B0"]][names(init[["B0"]]) %in% paste0(c("St", "Sts"), j)] = 
          prior[!is.na(eval(parse(text = paste0("seasonal", floor(j))))), c(paste0("seasonal", floor(j))), with = FALSE][[1]][1]
      }
    }
  }
  return(init)
}
