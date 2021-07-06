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
#' Can set to FALSE if want no seasonality
#' @param cycle, The period for the longer-term cycle. Deafult is NULL and will be estimated via wavelet analysis.
#' Can set to FALSE if want no cycle.
#' @param decomp Decomposition model ("tend-cycle-seasonal", "trend-seasonal", "trend-cycle", "trend-noise")
#' @param trend Trend specification ("random-walk", "random-walk-drift", "double-random-walk", "random-walk2"). The default is NULL which will choose the best of all specifications based on the maximum likelihood.
#' "random-walk" is the random walk trend.
#' "random-walk-drift" is the random walk with constant drift trend.
#' "double-random-walk" is the random walk with random walk drift trend.
#' "random-walk2" is a 2nd order random walk trend as in the Hodrick-Prescott filter.
#' If trend is "random-walk", the trend model is T_t = T_{t-1} + e_t, 
#' e_t ~ N(0, sig_t^2)
#' If trend is "random-walk-drift", the trend model is T_t = T_{t-1} + D_{t-1} + e_t, 
#' e_t ~ N(0, sig_t^2) with
#' D_t = d + phi_d*D_{t-1} + n_t, n_t ~ N(0, sig_d^2)
#' If trend is "double-random-walk", the trend model is T_t = M_{t-1} + T_{t-1} + e_t, 
#' e_t ~ N(0, sig_t^2) with
#' M_t = M_{t-1} + n_t, n_t ~ N(0, sig_d^2)
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
#' If det_drift = TRUE then the error variance of the drift equation (sig_d) is set to 0 and 
#' is refereed to as a deterministic drift
#' @param maxit Maximum number of iterations for the optimization
#' @param par Initial parameters, default is NULL
#' @param sig_level Significance level to determine statistically significant seasonal frequencies
#' @param verbose Logical whether to print messages or not
#' @param unconstrained Logical whether to remove inequality constraints on the trend during estimation
#' @param saturating_growth Force the growth rate to converge to 0 in the long term 
#' @param cores Number of cores to use for seasonality and cycle detection
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
                         multiplicative = NULL, par = NULL, seasons = NULL, cycle = NULL, cores = NULL,
                         det_obs = FALSE, det_trend = NULL, det_seas = FALSE, det_drift = FALSE, det_cycle = FALSE,
                         sig_level = 0.01, optim_methods = c("BFGS", "NM", "CG", "SANN"), maxit = 10000, verbose = FALSE){
  #exo = freq = decomp = trend = multiplicative = par = seasons = cycle = det_trend = cores = NULL
  #det_obs = det_seas = det_drift = det_cycle = unconstrained = saturating_growth = FALSE
  #sig_level = 0.01
  #optim_methods = "BFGS"
  #maxit = 10000
  #verbose = TRUE
  # for(i in list.files(path = "R", pattern = ".R", full.names = T)){
  #   tryCatch(source(i), error = function(err){NULL})
  # }
  # Rcpp::sourceCpp("src/kalmanfilter.cpp")
  
  #Argument checks
  if(sig_level <= 0 | sig_level > 0.1){
    stop("sig_level must be > 0 and <= 0.1.")
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
  if(!all(sapply(c(det_obs, det_seas, det_drift, det_cycle, unconstrained, saturating_growth, verbose), is.logical))){
    stop("det_obs, det_seas, det_drift, det_cycle, unconstrained, saturating_growth must be logical.")
  }
  if(!is.null(seasons)){
    if(!(is.logical(seasons) | is.numeric(seasons))){
      stop("seasons must be NULL, a numeric vector, or logical")
    }
  }
  if(!is.null(cycle)){
    if(!(is.logical(cycle) | is.numeric(cycle))){
      stop("cycle must be NULL, a numeric vector of length 1, or logical")
    }
    if(length(cycle) > 1){
      stop("cycle can only a numeric vector of length 1")
    }
  }
  if(!is.null(freq)){
    if(!is.numeric(freq)){
      stop("freq must be numeric")
    }
  }
  if(!is.null(cores)){
    if(cores > parallel::detectCores()){
      cores = parallel::detectCores()
      message("cores was set to be more than the available number of cores on the machine. Setting 'cores' to parallel::detectCores().")
    }
  }
  stsm_check_y(y)
  stsm_check_exo(exo, y)
  
  #Get the frequency of the data
  y = stsm_detect_frequency(y, freq)
  y = stsm_build_dates(y)
  dates = y$dates
  freq = y$freq
  freq_name = y$name
  standard_freq = y$standard_freq
  y = y$data
  
  #Remove leading and trailing NAs
  range = which(!is.na(y))
  y = unname(y[range[1]:range[length(range)]])
  dates = dates[range[1]:range[length(range)]]
  exo = stsm_format_exo(exo, dates, range)
  
  #Set up parallel computing
  if(is.null(cores)){
    cores = parallel::detectCores()
  }
  cl = tryCatch(parallel::makeCluster(max(c(1, cores))),
                error = function(err){
                  message("Parallel setup failed. Using single core.")
                  return(NULL)
                })
  if(!is.null(cl)){
    doSNOW::registerDoSNOW(cl)
  }else{
    cores = 1
  }
  
  #Set the decomposition
  if(verbose == TRUE & (is.null(decomp) | is.null(seasons) | is.null(trend))){
    message("Detecting the appropriate decomposition...")
  }
  
  #Set the prior
  prior = stsm_prior(y, freq)
  
  #Detect seasonality
  if(is.null(seasons) & ifelse(!is.null(decomp), grepl("seasonal", decomp), TRUE)){
    if(verbose == TRUE){
      message("Checking for seasonality...")
    }
    seasons = stsm_detect_seasonality(y, freq, sig_level, prior, cl, cores, show_progress = verbose)
  }
  if(is.null(seasons) | all(seasons == FALSE)){
    seasons = numeric(0)
  }
  
  #Detect cycle
  if(is.null(cycle) & ifelse(!is.null(decomp), grepl("cycle", decomp), TRUE) & length(y) >= 3*freq){
    if(verbose == TRUE){
      message("Checking for a cycle...")
    }
    cycle = stsm_detect_cycle(y, freq, sig_level, prior, cl , cores, show_progress = verbose)
  }
  if(is.null(cycle) | all(cycle == FALSE)){
    cycle = numeric(0)
  }
  
  #Reset the prior
  if(length(seasons) > 0 | length(cycle) > 0){
    prior = stsm_prior(y, freq, seasons = seasons, cycle = cycle)
  }
  
  #Detect multiplicative model
  if(is.null(multiplicative)){
    if(verbose == TRUE){
      message("Checking for a multiplicative model...")
    }
    multiplicative = stsm_detect_multiplicative(y, freq, sig_level, prior)
  }
  if(multiplicative == TRUE){
    y = log(y)
    prior = stsm_prior(y, freq, seasons = seasons, cycle = cycle)
  }
  
  #Detect trend type
  if(is.null(trend)){
    if(verbose == TRUE){
      message("Selecting the trend type...")
    }
    trend = stsm_detect_trend(y, freq, sig_level = sig_level, prior = prior, seasons = seasons, cycle = cycle)
    if(is.null(det_trend)){
      det_trend = trend$det_trend
    }
    trend = trend$trend
  }
  if(is.null(det_trend)){
    det_trend = FALSE
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
  
  #Standardize and reset the prior
  mean = mean(y, na.rm = TRUE)
  sd = stats::sd(y, na.rm = TRUE)
  y = (y - mean)/sd
  prior = stsm_prior(y, freq, decomp, seasons, cycle)
  
  #Set the initial parameter values
  if(is.null(par)){
    par = stsm_init_pars(y, freq, trend, cycle, decomp, seasons, prior, sig_level)
  }
  
  #Set any fixed parameters
  fixed = stsm_fixed_pars(par, y, det_obs, det_trend, det_drift, det_cycle, 
                          det_seas, saturating_growth, exo)
  par = fixed[["par"]]
  X = fixed[["X"]]
  fixed = fixed[["fixed"]]
  
  ###### Initial values for the unobserved components #####
  ssm = stsm_ssm(par, y, decomp, trend, init = NULL, prior = prior, seasons = seasons)
  init = ssm[c("B0", "P0")]
  par = c(par, P0 = stats::var(y, na.rm = T)) #Set uncertainty to be the rescaled variance of the time series
  
  ##### Inequality constraints: ineqA %*% par + ineqB > 0 => ineqA %*% par > -ineqB #####
  constraints = stsm_constraints(prior, par, freq, unconstrained, det_trend, det_drift, det_cycle, det_seas, det_obs, saturating_growth)
  par = constraints[["par"]]
  constraints = constraints[["constraints"]]
  #Test that constraints hold for initial parameter values: constraints[[1]] %*% matrix(par, ncol = 1) > -constraints[[2]]
  
  #Define the objective function
  objective = function(par, yt, freq, decomp, trend, init){
    ssm = stsm_ssm(par, yt, decomp, trend, init)
    ans = kalman_filter(ssm, yt, X, smooth = FALSE)
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
                                  yt = matrix(y, nrow = 1), freq = freq, decomp = decomp, trend = trend, init = init),
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
  
  #Stop the cluster
  if(!is.null(cl)){
    parallel::stopCluster(cl)
  }
  
  #Retrieve the model output
  out$estimate[grepl("sig_", names(out$estimate)) | names(out$estimate) == "P0"] = 
    abs(out$estimate[grepl("sig_", names(out$estimate)) | names(out$estimate) == "P0"])
  k = length(out$estimate) - 1
  n = length(y[!is.na(y)])
  fit = suppressWarnings(data.table(trend = trend, freq = freq, freq_name = freq_name, standard_freq = standard_freq,
                                    seasons = paste(seasons, collapse = ", "),
                                    cycle = 2*pi/out$estimate[grepl("lambda", names(out$estimate))],
                                    decomp = decomp, multiplicative = multiplicative, 
                                    convergence = (out$code == 0), 
                                    loglik = out$maximum, BIC = k*log(n) - 2*stats::logLik(out),
                                    AIC = stats::AIC(out), AICc = stats::AIC(out) + (2*k^2 + 2*k)/(n - k - 1),
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
