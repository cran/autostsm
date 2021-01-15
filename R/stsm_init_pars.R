#' Get initial parameter estimates for estimation
#'
#' @param y an object created from stsm_detect_frequency
#' @param freq Frequency of the data
#' @param trend Trend specification ("random-walk", "random-walk-drift", "double-random-walk", "random-walk2"). 
#' @param cycle The period for the longer-term cycle
#' @param decomp Decomposition model ("tend-cycle-seasonal", "trend-seasonal", "trend-cycle", "trend-noise")
#' @param harmonics The harmonics to split the seasonality into
#' @param prior A data table created by stsm_prior
#' @import data.table
#' @return named vector containing the initial parameter estimates for estimation
stsm_init_pars = function(y, freq, trend, cycle, decomp = "", harmonics = NULL, prior = NULL){
  if(is.null(prior)){
    prior = stsm_prior(y, freq, decomp, harmonics) 
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
              sig_m = unname(arima["sigma"]), d = unname(arima["intercept"]), phi = unname(arima["ar1"]))
      rm(arima)
    }else if(trend == "double-random-walk"){
      par = c(sig_t = stats::sd(diff(prior$trend), na.rm = TRUE)/sqrt(2),
              sig_m = stats::sd(diff(prior$trend), na.rm = TRUE)/sqrt(2))
    }
  }
  
  arima = tryCatch(forecast::Arima(stats::ts(prior$seasonal_adj - prior$trend, frequency = freq), order = c(2, 0, 1)), 
                   error = function(err){NULL})
  #Cycle parameters
  if(grepl("cycle", decomp)){
    if(is.null(cycle)){
      cycle = freq*5
    }
    #Set rho to be related to the speed of convergence of cycle towards steady state: i.e. maximum eigenvalue of ARMA(2,1)
    if(!is.null(arima)){
      Fm = rbind(arima$coef[c("ar1", "ar2", "ma1")], 
                 c(1, 0, 0), 
                 c(0, 0, 0))
      rho = max(Mod(eigen(Fm)$values))
      rho = max(c(0.01, min(c(rho, 1))))
    }else{
      rho = 0.95
    }
    par = c(par, rho = rho, lambda = 2*pi/cycle, sig_c = stats::sd(diff(prior$cycle), na.rm = TRUE))
  }else{
    if(!is.null(arima) & length(arima$coef) > 0){
      par = c(par, arima$coef, sig_c = sqrt(arima$sigma2))
    }else{
      par = c(par, ar1 = 1.25, ar2 = -0.3, ma1 = 0.75, sig_c = stats::sd(diff(prior$cycle)))
    }
    par = par[names(par) != "ar2"]
  }
  
  #Seasonal parameters
  if(grepl("seasonal", decomp)){
    par = c(par, sig_s = unname(rep(stats::sd(diff(prior$seasonal), na.rm = TRUE)/sqrt(length(harmonics)), length(harmonics))))
    names(par)[grepl("sig_s", names(par))] = paste0("sig_s", harmonics)
  }
  par = c(par, sig_e = stats::sd(prior$remainder, na.rm = TRUE))
  if(any(par[grepl("sig_", names(par))] <= 0)){
    par[grepl("sig_", names(par))] = ifelse(stats::sd(prior$remainder, na.rm = TRUE) > 0, stats::sd(prior$remainder, na.rm = TRUE), 0.01)
  }
  return(par)
}

# stsm_init_pars = function(y, freq, trend, cycle, decomp = "", harmonics = NULL, prior = NULL){
#   if(is.null(prior)){
#     prior = stsm_prior(y, freq, decomp, harmonics) 
#   }else{
#     prior = copy(prior)
#   }
#   
#   #Trend and drift parameters
#   if(trend == "random-walk"){
#     par = c(sig_t = sd(diff(prior$trend), na.rm = TRUE))
#   }else if(trend == "random-walk2"){
#     par = c(sig_t = sd(diff(diff(prior$trend)), na.rm = TRUE))
#   }else if(trend %in% c("random-walk-drift", "double-random-walk")){
#     if(trend == "random-walk-drift"){
#       arima = tryCatch(forecast::Arima(prior$drift, order = c(1, 0, 0), include.constant = TRUE),
#                        error = function(err){NULL})
#       if(!is.null(arima)){
#         arima = c(arima$coef, sigma = sqrt(arima$sigma2))
#         #Keep initial parameters away from the feasible borders
#         arima["ar1"] = ifelse(arima["ar1"] > 0.99, 0.95, 
#                               ifelse(arima["ar1"] < 0, 0.05, arima["ar1"]))
#         arima["intercept"] = (1 - arima["ar1"])*mean(prior$drift, na.rm = TRUE)
#       }else{
#         arima = c(intercept = mean(prior$drift, na.rm = TRUE), ar1 = 0.75) 
#         arima = c(arima, sigma = sqrt(var(prior$drift, na.rm = TRUE)*(1 - arima["ar1"]^2)))
#       }
#       par = c(sig_m = unname(arima["sigma"]), d = unname(arima["intercept"]), phi = unname(arima["ar1"]))
#       par["sig_t"] = sqrt(abs(var(diff(prior$trend), na.rm = TRUE) - par["sig_m"]^2/(1 - par["phi"]^2)))
#       rm(arima)
#     }else if(trend == "double-random-walk"){
#       par = c(sig_t = sd(diff(diff(prior$trend)), na.rm = TRUE)/sqrt(3),
#               sig_m = sd(diff(diff(prior$trend)), na.rm = TRUE)/sqrt(3))
#     }
#   }
#   
#   arima = tryCatch(forecast::Arima(stats::ts(prior$seasonal_adj - prior$trend, frequency = freq), order = c(2, 0, 1)), 
#                    error = function(err){NULL})
#   #Cycle parameters
#   if(grepl("cycle", decomp)){
#     if(is.null(cycle)){
#       cycle = freq*5
#     }
#     #Set rho to be related to the speed of convergence of cycle towards steady state: i.e. maximum eigenvalue of ARMA(2,1)
#     if(!is.null(arima)){
#       Fm = rbind(arima$coef[c("ar1", "ar2", "ma1")],
#                  c(1, 0, 0),
#                  c(0, 0, 0))
#       rho = max(Re(eigen(Fm)$values))
#       rho = max(c(0.01, min(c(rho, 1))))
#     }else{
#       rho = 0.95
#     }
#     par = c(par, rho = rho, lambda = 2*pi/cycle, 
#             sig_c = ifelse(!is.null(arima), sqrt(arima$sigma2), sd(diff(prior$cycle), na.rm = TRUE)))
#   }else{
#     if(!is.null(arima) & length(arima$coef) > 0){
#       par = c(par, arima$coef, sig_c = sqrt(arima$sigma2))
#     }else{
#       par = c(par, ar1 = 1.25, ar2 = -0.3, ma1 = 0.75, sig_c = sd(diff(prior$cycle), na.rm = TRUE))
#     }
#   }
#   
#   #Seasonal parameters
#   if(grepl("seasonal", decomp)){
#     par = c(par, sig_s = unname(rep(sd(diff(prior$seasonal), na.rm = TRUE)/sqrt(length(harmonics)), length(harmonics))))
#     names(par)[grepl("sig_s", names(par))] = paste0("sig_s", harmonics)
#   }
#   
#   #Capture any non-zero correlations in the error term
#   cov = cov(diff(stats::ts(prior[, c("trend", "cycle", "remainder", "seasonal")])), use = "all.obs")
#   par = c(par, sig_e = sqrt(abs(var(prior$remainder, na.rm = TRUE) + sum(cov[lower.tri(cov, diag = FALSE)]))))
#   if(any(par[grepl("sig_", names(par))] <= 0)){
#     par[grepl("sig_", names(par))] = ifelse(sd(prior$remainder, na.rm = TRUE) > 0, sd(prior$remainder, na.rm = TRUE), 0.01)
#   }
#   
#   return(par)
# }
