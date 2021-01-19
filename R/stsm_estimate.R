#' Trend cycle seasonal decomposition using the Kalman filter.
#'
#' Estimates a structural time series model using the Kalman filter and maximum likelihood.
#' The seasonal and cycle components are assumed to be of a trigonometric form.
#' The function checks three trend specifications to decompose a univariate time series
#' into trend, cycle, and/or seasonal components plus noise. The function automatically
#' detects the frequency and checks for a seasonal and cycle component if the user does not specify
#' the frequency or decomposition model. This can be turned off by setting freq or specifying decomp.
#' State space model for decomposition follows
#' Yt = T_t + C_t + S_t + A*X_t + e_t, e_t ~ N(0, sig_e^2)
#' Y is the data
#' T is the trend component
#' C is the cycle component
#' S is the seasonal component
#' X is the exogenous data with parameter vector B
#' e is the observation error
#' @param y Univariate time series of data values. May also be a 2 column data frame containing a date column.
#' @param exo Matrix of exogenous variables. Can be used to specify regression effects or other seasonal effects like holidays, etc.
#' @param freq Frequency of the data (1 (yearly), 4 (quarterly), 12 (monthly), 365.25/7 (weekly), 365.25 (daily)), default is NULL and will be automatically detected
#' @param seasons The seasonal periods: i.e. c(365.25, 7 if yearly and weekly seasonality). Default is NULL and will be estimated via wavelet analysis.
#' @param harmonics Vector for number of cycles per year for all of the seasons: i.e. c(unique(1:floor(12/2)), unique(1:floor(7/2))).
#' @param cycle, The period for the longer-term cycle
#' The default is NULL, which will use 1 harmonic per season (freq/seasons) if seasons is estimated via wavelet analysis for parsimony.
#' Otherwise it will use c(1, 2, 4, 6, 12, 24, 365.25/7, 365.25/3.5, 365.25/2) to capture harmonics at
#' yearly, semi-yearly, quarterly, every other month, monthly, twice per month, weekly, twice per week, and every other day to fully parameterize yearly, quarterly, monthly, and weekly
#' but will keep only harmonics > 0 and <= freq/2. This is more parsimonious than using harmonics at 1:floor(seasons/2).
#' @param decomp Decomposition model ("tend-cycle-seasonal", "trend-seasonal", "trend-cycle", "trend-noise")
#' @param trend Trend specification ("random-walk", "random-walk-drift", "double-random-walk", "random-walk2"). The default is NULL which will choose the best of all specifications based on the maximum likelihood.
#' "random-walk" is the random walk trend.
#' "random-walk-drift" is the random walk with constant drift trend.
#' "double-random-walk" is the random walk with random walk drift trend.
#' "random-walk2" is a 2nd order random walk trend as in the Hodrick-Prescott filter.
#' If trend is "random-walk", the trend model is T_t = T_{t-1} + e_t, 
#' e_t ~ N(0, sig_t^2)
#' If trend is "random-walk-drift", the trend model is T_t = T_{t-1} + M_{t-1} + e_t, 
#' e_t ~ N(0, sig_t^2) with
#' M_t = d + phi*M_{t-1} + n_t, n_t ~ N(0, sig_m^2)
#' If trend is "double-random-walk", the trend model is T_t = M_{t-1} + T_{t-1} + e_t, 
#' e_t ~ N(0, sig_t^2) with
#' M_t = M_{t-1} + n_t, n_t ~ N(0, sig_m^2)
#' If trend is "random-walk2", the trend model is T_t = 2T_{t-1} - T_{t-2} + e_t, 
#' e_t ~ N(0, sig_t^2)
#' @param multiplicative If data should be logged to create a multiplicative model.
#' If multiplicative = TRUE, then the data is logged and the original model becomes multiplicative 
#' (Y_t = T_t * C_t * S_t * BX_t * e_t)
#' @param optim_methods Vector of 1 to 3 optimization methods in order of preference ("NR", "BFGS", "CG", "BHHH", or "SANN")
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
#' If det_drift = TRUE then the error variance of the drift equation (sig_m) is set to 0 and 
#' is refereed to as a deterministic drift
#' @param maxit Maximum number of iterations for the optimization
#' @param par Initial parameters, default is NULL
#' @param sig_level Significance level to determine statistically significant seasonal frequencies
#' @param verbose Logical whether to print messages or not
#' @param unconstrained Logical whether to remove inequality constraints on the trend during estimation
#' @param saturating_growth Force the growth rate to converge to 0 in the long term 
#' @param full_seasonal_spectrum Whether to check for full spectrum of seasonality rather than weekly and yearly suggestion
#' @param seasonal_spectrum, Vector of seasonal periods to check
#' @import data.table
#' @useDynLib autostsm, .registration=TRUE
#' @return List of estimation values including a data table with coefficients, convergence code, frequency, decomposition, seasonality, cyclicality, and trend specification
#' as well as the a data table with the original data with dates. Any exogenous data given is also returned. 
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
#' stsm = stsm_estimate(NA000334Q)
#' }
#' @export
stsm_estimate = function(y, exo = NULL, freq = NULL, decomp = NULL, trend = NULL, unconstrained = FALSE, saturating_growth = FALSE,
                         multiplicative = NULL, par = NULL, seasons = NULL, harmonics = NULL, cycle = NULL,
                         det_obs = FALSE, det_trend = NULL, det_seas = FALSE, det_drift = FALSE, det_cycle = FALSE,
                         sig_level = 0.01, optim_methods = c("BFGS", "NM", "CG", "SANN"), maxit = 10000, verbose = FALSE,
                         full_seasonal_spectrum = FALSE, seasonal_spectrum = NULL){
  #exo = freq = decomp = trend = multiplicative = par = seasons = harmonics = cycle = seasonal_spectrum = det_trend = NULL
  #det_obs = det_seas = det_drift = det_cycle = unconstrained = saturating_growth = full_seasonal_spectrum = FALSE
  #sig_level = 0.01
  #optim_methods = "BFGS"
  #maxit = 10000
  #verbose = TRUE
  if(sig_level < 0.01 | sig_level > 0.1){
    stop("sig_level must be between 0.01 and 0.1.")
  }
  if(any(!optim_methods %in% c("NR", "BFGS", "BHHH", "SANN", "CG", "NM")) | length(optim_methods) < 1){
    stop("optim_methods must be a vector containing 'NR', 'BFGS', 'BHHH', 'SANN', 'CG', and/or 'NM'")
  }
  if(!is.numeric(maxit)){
    stop("maxit must be numeric and greater than 0.")
  }else if(maxit <= 0){
    stop("maxit must be numeric and greater than 0.")
  }
  if(!is.null(decomp)){
    if(!decomp %in% c('trend-noise', 'trend-cycle', 'trend-seasonal', 'trend-cycle-seasonal', 'trend-seasonal-cycle')){
      stop("decomp must be one of 'trend-noise', 'trend-cycle', 'trend-seasonal', 'trend-cycle-seasonal', or 'trend-seasonal-cycle'.")
    }
  }
  
  #Get the frequency of the data
  y = stsm_detect_frequency(y, freq)
  y = stsm_build_dates(y)
  dates = y$dates
  freq = y$freq
  y = y$data
  
  #Remove leading and trailing NAs
  range = which(!is.na(y))
  y = unname(y[range[1]:range[length(range)]])
  dates = dates[range[1]:range[length(range)]]
  if(!is.null(exo)){
    exo = exo[range[1]:range[length(range)], ]
  }
  
  #Set the decomposition
  if(verbose == TRUE & (is.null(decomp) | is.null(seasons) | is.null(trend))){
    message("Detecting the appropriate decomposition...")
  }
  
  #Detect trend type
  if(is.null(trend)){
    if(verbose == TRUE){
      message("Detecting the trend type...")
    }
    trend = stsm_detect_trend(y, freq, sig_level = sig_level)
    if(is.null(det_trend)){
      det_trend = trend$det_trend
    }
    trend = trend$trend
  }
  if(is.null(det_trend)){
    det_trend = FALSE
  }
  
  #Detect multiplicative model
  if(is.null(multiplicative)){
    if(verbose == TRUE){
      message("Detecting if the model should be multiplicative...")
    }
    multiplicative = stsm_detect_multiplicative(y, freq, sig_level)
  }
  if(multiplicative == TRUE){
    y = log(y)
  }
  
  #Set the prior
  prior = stsm_prior(y, freq)
  
  #Detect seasonality
  if(is.null(seasons) & is.null(harmonics) & ifelse(!is.null(decomp), grepl("seasonal", decomp), TRUE)){
    if(verbose == TRUE){
      message("Detecting seasonality...")
    }
    seasons = stsm_detect_seasonality(y, freq, sig_level = 0.0001, prior, full = full_seasonal_spectrum, spectrum = seasonal_spectrum)
    wavelet = TRUE
  }else if(is.null(seasons) & !is.null(harmonics)){
    seasons = freq/harmonics
    wavelet = FALSE
  }else{
    seasons = numeric(0)
    harmonics = numeric(0)
    wavelet = FALSE
  }

  #Detect cyclicality
  if(is.null(cycle) & ifelse(!is.null(decomp), grepl("cycle", decomp), TRUE) & length(y) >= 3*freq){
    if(verbose == TRUE){
      message("Detecting cyclicallity...")
    }
    cycle = stsm_detect_cycle(y, freq, sig_level = 0.0001, prior)
  }else{
    cycle = numeric(0)
  }
  
  #Assign the decomposition
  if(is.null(decomp)){
    decomp = "trend"
    if(length(seasons) > 0){
      decomp = paste0(decomp, "-seasonal")
    }
    if(length(cycle) > 0){
      decomp = paste0(decomp, "-cycle")
    }
  }
  
  #Remove seasonal and cycle from decomp if harmonics and cycle are missing
  if(grepl("seasonal", decomp) & (is.null(seasons) | length(seasons) == 0)){
    decomp = strsplit(decomp, "-")[[1]]
    decomp = paste(decomp[decomp != "seasonal"], collapse = "-")
  }
  if(grepl("cycle", decomp) & (is.null(cycle) | length(cycle) == 0)){
    decomp = strsplit(decomp, "-")[[1]]
    decomp = paste(decomp[decomp != "cycle"], collapse = "-")
  }
  
  #If no seasonality or cycle detected, include noise component only
  if(decomp == "trend"){
    decomp = "trend-noise"
  }
  
  #Standardize
  mean = mean(y, na.rm = TRUE)
  sd = stats::sd(y, na.rm = TRUE)
  y = (y - mean)/sd
  
  #Get seasonal harmonics if not already set
  if(is.null(harmonics)){
    harmonics = stsm_harmonics(freq, decomp, seasons, wavelet)
  }else{
    harmonics = sort(unique(unlist(harmonics)))
  }
  
  ##### Set up the initial values for estimation #####
  #Use naive model as the prior
  if(floor(freq) == 260){
    #Ensure the weekday and daily model have the same starting values
    y_use = stsm_detect_frequency(data.table(y = y, date = dates), 365.25)
    y_use = stsm_build_dates(y_use)
    prior = stsm_prior(y_use$data, y_use$freq, decomp, harmonics*7/5)
    rm(y_use)
  }else{
    prior = stsm_prior(y, freq, decomp, harmonics)
  }
  
  #Set the initial parameter values
  if(is.null(par)){
    par = stsm_init_pars(y, freq, trend, cycle, decomp, harmonics, prior)
  }
  
  #Set any fixed parameters
  fixed = NULL
  if(det_obs == TRUE){
    if("sig_c" %in% names(par)){
      par["sig_c"] = par["sig_c"] + par["sig_e"]
    }else if("sig_m" %in% names(par) & "sig_t" %in% names(par)){
      par["sig_m"] = par["sig_m"] + par["sig_e"]/2
      par["sig_t"] = par["sig_t"] + par["sig_e"]/2
    }else if("sig_m" %in% names(par) & !"sig_t" %in% names(par)){
      par["sig_m"] = par["sig_m"] + par["sig_e"]
    }else if("sig_t" %in% names(par) & !"sig_m" %in% names(par)){
      par["sig_t"] = par["sig_t"] + par["sig_e"]
    }else if("sig_c" %in% names(par)){
      par["sig_c"] = par["sig_c"] + par["sig_e"]
    }
    par["sig_e"] = 0
    fixed = c(fixed, "sig_e")
  }
  if(det_trend == TRUE){
    if("sig_m" %in% names(par)){
      par["sig_m"] = par["sig_m"] + par["sig_t"]
    }
    par["sig_t"] = 0
    fixed = c(fixed, "sig_t")
  }
  if(det_drift == TRUE){
    if("sig_t" %in% names(par) & det_trend == FALSE){
      par["sig_t"] = par["sig_t"] + par["sig_m"]
    }
    par["sig_m"] = 0
    fixed = c(fixed, "sig_m")
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
  
  ###### Initial values for the unobserved components #####
  init = stsm_init_vals(y = y, par = par, freq = freq, trend = trend, 
                        decomp = decomp, harmonics = harmonics, prior = prior)
  par = c(par, P0 = 1)
  #Set uncertainty to be the rescaled variance of the time series
  
  ##### Inequality constraints: ineqA %*% par + ineqB > 0 => ineqA %*% par > -ineqB #####
  constraints = stsm_constraints(prior, par, freq, unconstrained, det_trend, det_drift, det_cycle, det_seas, det_obs, saturating_growth)
  par = constraints[["par"]]
  constraints = constraints[["constraints"]]
  #Test that constraints hold for initial parameter values: constraints[[1]] %*% matrix(par, ncol = 1) > -constraints[[2]]

  #Define the objective function
  objective = function(par, y, freq, decomp, trend, init){
    yt = matrix(y, nrow = 1)
    sp = stsm_ssm(par = par, yt = yt, freq = freq, decomp = decomp, trend = trend, init = init)
    ans = kalman_filter(B0 = matrix(sp$B0, ncol = 1), P0 = sp$P0, Dt = sp$Dt, At = sp$At, Ft = sp$Ft, Ht = sp$Ht, Qt = sp$Qt, Rt = sp$Rt,
                        yt = yt, X = X, beta = sp$beta)
    return(ans$loglik)
  }
  
  #Estimate the model
  for(o in optim_methods){
    if(verbose == TRUE){
      message("Estimating with method ", o)
    }
    out = tryCatch(maxLik::maxLik(logLik = objective,
                                  start = par, method = o, fixed = fixed,
                                  finalHessian = FALSE, hess = NULL, control = list(printLevel = ifelse(verbose == TRUE, 2, 0), iterlim = maxit),
                                  constraints = constraints,
                                  y = y, freq = freq, decomp = decomp, trend = trend, init = init),
                   error = function(err){NULL})
    if(!is.null(out)){
      break
    }else{
      message("Estimation method ", o, " failed")
    }
  }
  if(is.null(out)){
    stop("Estimation failed.")
  }
  
  #Retrieve the model output
  out$estimate[grepl("sig_", names(out$estimate)) | names(out$estimate) == "P0"] = 
    abs(out$estimate[grepl("sig_", names(out$estimate)) | names(out$estimate) == "P0"])
  k = length(out$estimate) - 1
  n = length(y[!is.na(y)])
  fit = suppressWarnings(data.table(trend = trend, freq = freq, 
                               seasons = paste(seasons, collapse = ", "),
                               harmonics = ifelse(!is.null(harmonics), paste(harmonics, collapse = ", "), NA),
                               cycle = 2*pi/out$estimate[grepl("lambda", names(out$estimate))],
                               decomp = decomp, multiplicative = multiplicative, 
                               convergence = (out$code == 0), 
                               loglik = out$maximum, BIC = k*log(n) - 2*stats::logLik(out),
                               AIC = stats::AIC(out), AICc = stats::AIC(out) + (2*k^2 + 2*k)/(length(y) - k - 1),
                               coef = paste(paste(names(stats::coef(out)), unname(stats::coef(out)), sep = " = "), collapse = ", ")
  ))
  return(fit)
}

########
#Call these to build the package
#devtools::document()
#devtools::build_vignettes()
#devtools::install()
#library(autostsm)
#git config remote.origin.url git@github.com:user/autostsm.git
