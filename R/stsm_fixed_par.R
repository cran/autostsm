#' Fixed parameter setting
#' 
#' @param par Initial parameters
#' @param y Vector of univariate time series
#' @param det_obs Set the observation equation error variance to 0 (deterministic observation equation)
#' If det_obs = TRUE then the error variance of the observation equation (sig_e) is set to 0
#' @param det_trend Set the trend error variance to 0 (deterministic trend)
#' If det_trend = TRUE then the error variance of the trend equation (sig_t) is set to 0 and 
#' is referred to as a smooth trend
#' @param det_seas Set the seasonality error variances to 0 (deterministic seasonality)
#' If det_seas = TRUE then the error variance all seasonality frequency j equations (sig_s) 
#' are set to 0 and is referred to as deterministic seasonality
#' @param det_cycle Set the cycle error variance to 0 (deterministic cycle)
#' If det_cycle = TRUE then the error variance of the cycle equation (sig_c) is set to 0 and 
#' is referred to as a deterministic cycle
#' @param det_drift Set the drift error variance to 0 (deterministic drift)
#' If det_drift = TRUE then the error variance of the drift equation (sig_d) is set to 0 and 
#' is refereed to as a deterministic drift
#' @param exo Matrix of exogenous variables. Can be used to specify regression effects or other seasonal effects like holidays, etc.
#' @param saturating_growth Force the growth rate to converge to 0 in the long term 
stsm_fixed_pars = function(par, y, det_obs = FALSE, det_trend = FALSE, det_drift = FALSE, 
                           det_cycle = FALSE, det_seas = FALSE, 
                           saturating_growth = FALSE, exo = NULL){
  #Set any fixed parameters
  fixed = NULL
  if(det_obs == TRUE){
    if("sig_c" %in% names(par)){
      par["sig_c"] = par["sig_c"] + par["sig_e"]
    }else if("sig_d" %in% names(par) & "sig_t" %in% names(par)){
      par["sig_d"] = par["sig_d"] + par["sig_e"]/2
      par["sig_t"] = par["sig_t"] + par["sig_e"]/2
    }else if("sig_d" %in% names(par) & !"sig_t" %in% names(par)){
      par["sig_d"] = par["sig_d"] + par["sig_e"]
    }else if("sig_t" %in% names(par) & !"sig_d" %in% names(par)){
      par["sig_t"] = par["sig_t"] + par["sig_e"]
    }else if("sig_c" %in% names(par)){
      par["sig_c"] = par["sig_c"] + par["sig_e"]
    }
    par["sig_e"] = 0
    fixed = c(fixed, "sig_e")
  }
  if(det_trend == TRUE){
    if("sig_d" %in% names(par)){
      par["sig_d"] = par["sig_d"] + par["sig_t"]
    }
    par["sig_t"] = 0
    fixed = c(fixed, "sig_t")
  }
  if(det_drift == TRUE){
    if("sig_t" %in% names(par) & det_trend == FALSE){
      par["sig_t"] = par["sig_t"] + par["sig_d"]
    }
    par["sig_d"] = 0
    fixed = c(fixed, "sig_d")
  }
  if(det_cycle == TRUE){
    if("sig_e" %in% names(par) & det_obs == FALSE){
      par["sig_e"] = par["sig_e"] + par["sig_c"]  
    }else if("sig_s" %in% names(par) & det_seas == FALSE){
      par["sig_s"] = par["sig_s"] + par["sig_c"]
    }
    par["sig_c"] = 0
    fixed = c(fixed, "sig_c")
  }
  if(det_seas == TRUE){
    if("sig_e" %in% names(par) & det_obs == FALSE){
      par["sig_e"] = par["sig_e"] + par["sig_s"]
    }else if("sig_c" %in% names(par) & det_cycle == FALSE){
      par["sig_c"] = par["sig_c"] + par["sig_s"]
    }
    par[grepl("sig_s", names(par))] = 0
    fixed = c(fixed, names(par)[grepl("sig_s", names(par))])
  }
  if(saturating_growth == TRUE){
    par[names(par) == "d"] = 0
    fixed = c(fixed, "d")
  }
  if(is.null(exo)){
    X = t(matrix(0, nrow = length(y), ncol = 1))
    rownames(X) = "X"
    par = c(par, beta_X = 0)
    fixed = c(fixed, "beta_X")
  }else{
    X = t(exo)
    par = c(par, beta_ = stats::coef(stats::lm(y ~ . - 1, data = data.frame(cbind(y, exo)))))
  }
  return(list(par = par, fixed = fixed, X = X))
}