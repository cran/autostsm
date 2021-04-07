#' Get initial parameter estimates for estimation
#'
#' @param y an object created from stsm_detect_frequency
#' @param freq Frequency of the data
#' @param trend Trend specification ("random-walk", "random-walk-drift", "double-random-walk", "random-walk2"). 
#' @param cycle The period for the longer-term cycle
#' @param decomp Decomposition model ("tend-cycle-seasonal", "trend-seasonal", "trend-cycle", "trend-noise")
#' @param seasons The seasonal lengths to split the seasonality into
#' @param prior A data table created by stsm_prior
#' @import data.table
#' @return named vector containing the initial parameter estimates for estimation
stsm_init_pars = function(y, freq, trend, cycle, decomp = "", seasons = NULL, prior = NULL){
  if(any(is.na(y))){
    y = suppressWarnings(imputeTS::na_kalman(y))
  }
  
  if(is.null(prior)){
    prior = stsm_prior(y, freq, decomp, seasons) 
  }else{
    prior = copy(prior)
  }
  
  #Trend and drift parameters
  if(trend == "random-walk"){
    par = c(sig_t = stats::sd(diff(prior$trend), na.rm = TRUE))
  }else if(trend == "random-walk2"){
    par = c(sig_t = stats::sd(diff(diff(prior$trend)), na.rm = TRUE))
  }else if(trend %in% c("random-walk-drift", "double-random-walk")){
    if(trend == "random-walk-drift"){
      arima = tryCatch(forecast::Arima(prior$drift, order = c(1, 0, 0), include.constant = TRUE),
                       error = function(err){NULL})
      if(!is.null(arima)){
        arima = c(arima$coef, sigma = sqrt(arima$sigma2))
        #Keep initial parameters away from the feasible borders
        arima["ar1"] = ifelse(arima["ar1"] > 0.99, 0.75, 
                              ifelse(arima["ar1"] < 0, 0.25, arima["ar1"]))
      }else{
        arima = c(sigma = stats::sd(prior$drift, na.rm = TRUE), intercept = mean(prior$drift, na.rm = TRUE)*(1 - 0.75), ar1 = 0.75) 
      }
      par = c(sig_t = unname(sqrt(abs(stats::var(diff(prior$trend), na.rm = TRUE) - arima["sigma"]^2/(1 - arima["ar1"]^2)))), 
              sig_d = unname(arima["sigma"]), d = unname(arima["intercept"]), phi_d = unname(arima["ar1"]))
      rm(arima)
    }else if(trend == "double-random-walk"){
      par = c(sig_t = stats::sd(diff(prior$trend), na.rm = TRUE)/sqrt(2),
              sig_d = stats::sd(diff(prior$trend), na.rm = TRUE)/sqrt(2))
    }
  }
  
  div = max(c(1, sum(sapply(c("cycle", "seasonal"), function(x){grepl(x, decomp)}))))
  arima = tryCatch(forecast::Arima(stats::ts(prior$cycle + prior$remainder/div, frequency = freq), order = c(2, 0, 1)), 
                   error = function(err){NULL})
  #Cycle parameters
  if(grepl("cycle", decomp)){
    if(is.null(cycle)){
      cycle = freq*5
    }
    #Set phi_c to be related to the speed of convergence of cycle towards steady state: i.e. maximum eigenvalue of ARMA(2,1)
    if(!is.null(arima)){
      Fm = rbind(arima$coef[c("ar1", "ar2", "ma1")], 
                 c(1, 0, 0), 
                 c(0, 0, 0))
      phi_c = max(Mod(eigen(Fm)$values))
      phi_c = max(c(0.01, min(c(phi_c, 1))))
      sig_c = sqrt(arima$sigma2)
    }else{
      phi_c = 0.95
      sig_c = stats::sd(diff(prior$cycle + prior$remainder/div), na.rm = TRUE)
    }
    par = c(par, phi_c = phi_c, lambda = 2*pi/cycle, sig_c = sig_c)
  }else{
    if(!is.null(arima) & length(arima$coef) > 0){
      par = c(par, arima$coef[names(arima$coef) != "intercept"], sig_c = sqrt(arima$sigma2))
      names(par)[grepl("ar\\d+", names(par))] = gsub("ar", "phi_c.", names(par)[grepl("ar\\d+", names(par))])
      names(par)[grepl("ma\\d+", names(par))] = gsub("ma", "theta_c.", names(par)[grepl("ma\\d+", names(par))])
    }else{
      par = c(par, phi_c.1 = 1.25, phi_c.2 = -0.3, theta_c.1 = 0.75, sig_c = stats::sd(diff(prior$cycle + prior$remainder/div)))
    }
    par = par[names(par) != "phi_c.2"]
  }
  
  #Seasonal parameters
  if(grepl("seasonal", decomp)){
    par = c(par, sig_s = unname(rep(stats::sd(diff(prior$seasonal + prior$remainder/div), na.rm = TRUE)/sqrt(length(seasons)), length(seasons))))
    names(par)[grepl("sig_s", names(par))] = paste0("sig_s", seasons)
  }
  par = c(par, sig_e = stats::sd(prior$remainder/div, na.rm = TRUE))
  if(any(par[grepl("sig_", names(par))] <= 0)){
    par[grepl("sig_", names(par))] = ifelse(stats::sd(prior$remainder, na.rm = TRUE) > 0, stats::sd(prior$remainder, na.rm = TRUE), 0.01)
  }
  return(par)
}


