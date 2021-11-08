#' Kalman Filter and Forecast
#'
#' Kalman filter and forecast an estimated model from stsm_estimate output
#' @param model Structural time series model estimated using stsm_estimate.
#' @param n.ahead Number of periods to forecast
#' @param ci Confidence interval, value between 0 and 1 exclusive.
#' @param y Univariate time series of data values. May also be a 2 column data frame containing a date column.
#' @param freq Frequency of the data (1 (yearly), 4 (quarterly), 12 (monthly), 365.25/7 (weekly), 365.25 (daily)), default is NULL and will be automatically detected
#' @param exo Matrix of exogenous variables used for the historical data. Can be used to specify regression effects or other seasonal effects like holidays, etc.
#' @param exo.fc Matrix of exogenous variables used for the forecast
#' @param plot, Logical, whether to plot everything
#' @param plot.decomp Logical, whether to plot the filtered historical data
#' @param plot.fc Logical, whether to plot the forecast
#' @param n.hist Number of historical periods to include in the forecast plot. If plot = TRUE and n.hist = NULL, defaults to 3 years.
#' @param smooth Whether or not to use the Kalman smoother
#' @param dampen_cycle Whether to remove oscillating cycle dynamics and smooth the cycle forecast into the trend using a sigmoid function that maintains the rate of convergence
#' @param envelope_ci Whether to create a envelope for the confidence interval to smooth out seasonal fluctuations
#' @import data.table ggplot2
#' @useDynLib autostsm, .registration=TRUE
#' @return data table (or list of data tables) containing the filtered and/or smoothed series.
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
#' fc = stsm_forecast(stsm, y = NA000334Q, n.ahead = floor(stsm$freq)*3, plot = TRUE)
#' }
#' @export
stsm_forecast = function(model, y, n.ahead = 0, freq = NULL, exo = NULL, exo.fc = NULL, ci = 0.8, 
                         plot = FALSE, plot.decomp = FALSE, plot.fc = FALSE, n.hist = NULL, smooth = TRUE, 
                         dampen_cycle = FALSE, envelope_ci = FALSE){
  #model = stsm
  #n.ahead = floor(model$freq*3)
  #ci = 0.8
  #freq = exo = exo.fc = n.hist = NULL
  #plot = plot.decomp = plot.fc = smooth = TRUE
  #dampen_cycle = FALSE
  # for(i in list.files(path = "R", pattern = ".R", full.names = T)){
  #   tryCatch(source(i), error = function(err){NULL})
  # }
  # Rcpp::sourceCpp("src/kalmanfilter.cpp")

  #Bind data.table variables to the global environment
  fev = extrema = fitted = variable = value = cycle = seasonal = observed = remainder = 
    forecast = group = Tt_0 = Ct_0 = Dt_0 = interpolated = NULL
  
  #Argument checks
  if(!all(sapply(c(plot, plot.decomp, plot.fc, smooth, dampen_cycle), is.logical))){
    stop("plot, plot.decomp, plot.fc, smooth, dampen_cycle must be logical")
  }
  if(n.ahead < 0 | (n.ahead %% 1) != 0){
    stop("n.ahead must be an integer > 0")
  }
  if(!is.null(n.hist)){
    if(n.hist <= 0){
      stop("n.hist must be > 0")
    }
  }
  if(ci <= 0 | ci >= 1){
    stop("ci must between 0 and 1. Suggested values are 0.8, 0.9, 0.95, 0.99.")
  }
  stsm_check_y(y)
  stsm_check_exo(exo, y)
  stsm_check_exo_fc(exo.fc, n.ahead)
  
  #Assign plot values
  if(plot == TRUE){
    plot.decomp = plot.fc = TRUE
  }
  if(plot.fc == TRUE & n.ahead == 0){
    plot.fc = FALSE
  }
  
  #Get the frequency of the data
  y = stsm_detect_frequency(y, freq)
  y = stsm_build_dates(y)
  dates = y$dates
  if(model$freq != y$freq){
    warning("Frequency of data is different than the model. Defaulting to the model frequency")
    freq = model$freq
    standard_freq = model$standard_freq
  }else{
    freq = y$freq
    standard_freq = y$standard_freq
  }
  y = y$data
  
  #Remove leading and trailing NAs
  range = which(!is.na(y))
  y = unname(y[range[1]:range[length(range)]])
  dates = dates[range[1]:range[length(range)]]
  exo = stsm_format_exo(exo, dates, range)
  
  #Interpolate dates
  if(!is.null(model$interpolate)){
    y = stsm_dates_to_interpolate(y = y, dates = dates, exo = exo, interpolate = model$interpolate)
    dates = y$dates
    exo = y$exo
    y = y$y
  }
  
  #Set the historical exogenous variables
  if(is.null(exo)){
    X = t(matrix(0, nrow = length(y), ncol = 1))
    X[is.na(X)] = 0
    rownames(X) = "X"
  }else{
    X = t(exo)
  }
  
  #Set the forecast exogenous variables
  if(!is.null(exo.fc)){
    if(is.null(exo)){
      warning("exo.fc is supplied but exo is not. Not using exo.")
    }else{
      exo.fc = as.data.table(exo.fc)
      exo.fc = exo.fc[, names(which(unlist(exo.fc[, lapply(.SD, is.numeric)]))), with = FALSE]
      X = cbind(X, t(exo.fc))
      if(ncol(X) != length(y) + n.ahead){
        stop("nrow of exo + exo.fc does not equal length of y + n.ahead.")
      }
    }
  }else{
    if(n.ahead > 0){
      X = cbind(X, matrix(0, ncol = n.ahead))
    }
  }
  
  #Set n.hist to the minimum of the length of the data or 3 times the frequency of the data if it is not specified
  #or longer than th length of the data and plotting is set to TRUE
  if(plot.fc == TRUE){
    if(is.null(n.hist) | ifelse(!is.null(n.hist), n.hist > length(y), FALSE)){
      n.hist = min(c(length(y), floor(freq*3)))
    }
  }
  
  #Apply multiplicative model
  if(model$multiplicative == TRUE){
    y = log(y)
  }
  
  #Standardize
  mean = mean(y, na.rm = TRUE)
  sd = stats::sd(y, na.rm = TRUE)
  y = (y - mean)/sd
  
  #Get the seasons
  if(!is.na(model$seasons)){
    seasons = as.numeric(strsplit(model$seasons, ", ")[[1]])
  }else{
    seasons = c()
  }
  
  #Get the coefficients
  par = eval(parse(text = paste0("c(", model$coef, ")")))
  
  #Filter and smooth the data
  ssm = stsm_ssm(yt = y, model = model)
  if(!is.null(model$interpolate)){
    int_per = ssm[["int_per"]]
    n.ahead = n.ahead*int_per
    if(is.null(exo.fc)){
      if(n.ahead > 0){
        X = cbind(X, matrix(0, ncol = n.ahead - (ncol(X) - length(y))))
      }
    }
  }
  msg = utils::capture.output(B_tt <- kalman_filter(ssm, matrix(y, nrow = 1), X, smooth)$B_tt, type = "message")
  rownames(B_tt) = rownames(ssm$Fm)

  #Get the unobserved series
  series = data.table(t(B_tt))
  series$fitted = c(matrix(ssm$Am, nrow = 1, ncol = nrow(series)) + 
                      ssm$Hm %*% t(as.matrix(series[, rownames(B_tt), with = FALSE])) + 
                      matrix(par[grepl("beta_", names(par))], nrow = 1) %*% X[, 1:length(y)])

  #Forecast
  if(n.ahead > 0){
    Fm_pow = diag(nrow(ssm$Fm))
    for(j in 1:n.ahead){
      Fm_pow = Fm_pow %*% ssm$Fm
      series = rbind(series, cbind(t(ssm$Dm + ssm$Fm %*% t(as.matrix(series[.N, rownames(B_tt), with = FALSE]))),
                                   c(ssm$Hm %*% (Fm_pow %*% ssm$Qm %*% t(Fm_pow)) %*% t(ssm$Hm))),
                     use.names = TRUE, fill = TRUE)
    }
    colnames(series)[ncol(series)] = "fev"
    series[!is.na(fev), "fev" := cumsum(fev) + c(ssm$Rm)]
    
    if(dampen_cycle == TRUE & is.complex(eigen(ssm$Fm[grepl("Ct|et", colnames(ssm$Fm)), grepl("Ct|et", rownames(ssm$Fm))])$values) & n.ahead > 1){
      #Smooth the cycle forecast into the trend using a sigmoid function so as not to get oscillating cycle forecasts
      Fm = ssm$Fm[grepl("Ct|et", colnames(ssm$Fm)), grepl("Ct|et", rownames(ssm$Fm))]
      B = matrix(unlist(series[length(y), colnames(series)[grepl("Ct|et", colnames(ssm$Fm))], with = FALSE]), ncol = 1)
      for(j in 1:n.ahead){
        B = cbind(B, Fm %*% B[, ncol(B)])
      }
      B = data.table(c = c(B[1, 2:ncol(B)]))
      B[, "extrema" := shift(abs(c), type = "lag", n = 1) < abs(c) & shift(abs(c), type = "lead", n = 1) < abs(c)]
      B[, "switch" := sign(c) != shift(sign(c), type = "lag", n = 1)]
      if(nrow(B[switch == TRUE | extrema == TRUE, ]) > 0){
        first_sign_switch = B[switch == TRUE, which = TRUE][1]
        first_sign_siwtch = ifelse(length(first_sign_switch) == 0, Inf, first_sign_switch)
        first_extrema = B[extrema == TRUE, which = TRUE][1]
        first_extrema = ifelse(length(first_extrema) == 0, Inf, first_extrema)
        wh1 = ifelse(first_sign_switch <= first_extrema, 1, first_extrema)
        B = B[wh1:.N, ]
        wh2 = B[switch == TRUE, which = TRUE][1]
        wh2 = ifelse(length(wh2) == 0, nrow(B), wh2)
        sigmoid = function(p, t){
          p[1] + (0 - p[1])/(1 + p[2]*exp(-p[3]*(t - p[4])))^(1/p[5])
        }
        t = 1:wh2
        obj = function(p, t){
          f1 = (B$c[1] - sigmoid(p, t[1]))^2
          f2 = sum((B$c[1:wh2] - sigmoid(p, t[1:wh2]))^2)
          return(f1 + f2)
        }
        out = stats::optim(c(B[1, 1], 1, 1, stats::median(t[1:wh2]), 1), fn = obj, method = "BFGS", t = t)
        z = sigmoid(out$par, 1:(n.ahead - wh1 + 1))
        series[(length(y) + wh1):.N, "Ct_0" := z]
      }
    }
    
    series[!is.na(fev), "fitted" := c(matrix(ssm$Am, nrow = 1, ncol = nrow(series[!is.na(fev), ])) +
                                        ssm$Hm %*% t(as.matrix(series[!is.na(fev), rownames(B_tt), with = FALSE])) + 
                                        matrix(par[grepl("beta_", names(par))], nrow = 1) %*% X[, (length(y) + 1):ncol(X)])]
    
    series[!is.na(fev), paste0((1/2 - ci/2)*100, "%") := stats::qnorm(1/2 - ci/2, fitted, sqrt(fev))]
    series[!is.na(fev), paste0((1/2 + ci/2)*100, "%") := stats::qnorm(1/2 + ci/2, fitted, sqrt(fev))]
    series[length(y), c(paste0((1/2 - ci/2)*100, "%"), paste0((1/2 + ci/2)*100, "%")) := fitted]
    
    #Create envelope around seasonal fluctuations in the confidence intervals
    if(envelope_ci == TRUE){
      if(!is.null(seasons) & all(!is.na(seasons)) & length(seasons) > 0){
        series[, "rn" := 1:.N]
        series = melt(series, id.vars = "rn")
        series[variable %in% paste0((1/2 + c(-ci, ci)/2)*100, "%") & !is.na(value),
               "smooth" := stats::predict(stats::smooth.spline(value))$y,
               by = "variable"]
        series[variable %in% paste0((1/2 + c(-ci, ci)/2)*100, "%") & !is.na(value),
                "shift" := frollapply(abs(value - smooth), align = "center", n = round(min(as.numeric(strsplit(model$seasons, ",")[[1]]))),
                                      FUN = max),
                by = "variable"]
        series[variable %in% paste0((1/2 + c(-ci, ci)/2)*100, "%") & !is.na(value),
                "shift" := nafill(shift, type = "locf"),
                by = "variable"]
        series = series[.N:1, ]
        series[variable %in% paste0((1/2 + c(-ci, ci)/2)*100, "%") & !is.na(value),
                "shift" := nafill(shift, type = "locf"),
                by = "variable"]
        series = series[.N:1, ]
        series[variable == paste0((1/2 + ci/2)*100, "%"), "value" := smooth + shift]
        series[variable == paste0((1/2 - ci/2)*100, "%"), "value" := smooth - shift]
        series = dcast(series, "rn ~ variable", value.var = "value")
        series[, "rn" := NULL]
      }
    }
    series[, "fev" := NULL]
    rm(Fm_pow)
  }
  
  #Build the forecast dates
  `%m+%` = lubridate::`%m+%`
  if(n.ahead > 0){
    y.fc = c(ssm$Am[1, ] + ssm$Hm %*% t(as.matrix(series[(length(y) + 1):.N, rownames(B_tt), with = FALSE])))
    if(standard_freq == TRUE){
      if(floor(freq) == floor(60*60*24*365.25)){
        #Secondly frequency
        dates.fc = dates[length(dates)] %m+% lubridate::seconds(1:n.ahead)
      }else if(floor(freq) == floor(60*60*24*365.25*5/7)){
        #Secondly frequency, weekdays only
        dates.fc = dates[length(dates)] %m+% lubridate::seconds(1:n.ahead*2)
        dates.fc = dates.fc[which(!weekdays(dates.fc) %in% c("Saturday", "Sunday"))][1:n.ahead]
      }else if(floor(freq) == floor(60*24*365.25)){
        #Minutely frequency
        dates.fc = dates[length(dates)] %m+% lubridate::minutes(1:n.ahead)
      }else if(floor(freq) == floor(60*24*365.25*5/7)){
        #Minutely frequency, weekday only
        dates.fc = dates[length(dates)] %m+% lubridate::minutes(1:n.ahead)
        dates.fc = dates.fc[which(!weekdays(dates.fc) %in% c("Saturday", "Sunday"))][1:n.ahead]
      }else if(floor(freq) == floor(24*365.25)){
        #Hourly frequency
        dates.fc = dates[length(dates)] %m+% lubridate::hours(1:n.ahead)
      }else if(floor(freq) == floor(24*365.25*5/7)){
        #Hourly frequency, weekday only
        dates.fc = dates[length(dates)] %m+% lubridate::hours(1:n.ahead)
        dates.fc = dates.fc[which(!weekdays(dates.fc) %in% c("Saturday", "Sunday"))][1:n.ahead]
      }else if(floor(freq) == 365){
        #Daily frequency
        dates.fc = dates[length(dates)] %m+% lubridate::days(1:n.ahead)
      }else if(floor(freq) == floor(365.25*5/7)){
        #Daily frequency, weekday only
        dates.fc = dates[length(dates)] %m+% lubridate::days(1:ceiling(n.ahead*2))
        dates.fc = dates.fc[which(!weekdays(dates.fc) %in% c("Saturday", "Sunday"))][1:n.ahead]
      }else if(floor(freq) == 52){
        #Weekly frequency
        if(is.null(model$interpolate)){
          dates.fc = dates[length(dates)] %m+% lubridate::weeks(1:n.ahead)
        }else if(model$interpolate == "daily"){
          dates.fc = dates[length(dates)] %m+% lubridate::days(1:n.ahead)
        }
      }else if(floor(freq) == 12){
        #Monthly frequency
        if(is.null(model$interpolate)){
          dates.fc = dates[length(dates)] %m+% months(1:n.ahead)
        }else if(model$interpolate == "weekly"){
          dates.fc = dates[length(dates)] %m+% lubridate::weeks(1:n.ahead)
        }else if(model$interpolate == "daily"){
          dates.fc = dates[length(dates)] %m+% lubridate::days(1:n.ahead)
        }
      }else if(floor(freq) == 4){
        #Quarterly frequency
        if(is.null(model$interpolate)){
          dates.fc = dates[length(dates)] %m+% months((1:n.ahead)*3)
        }else if(model$interpolate == "monthly"){
          dates.fc = dates[length(dates)] %m+% months(1:n.ahead)
        }else if(model$interpolate == "weekly"){
          dates.fc = dates[length(dates)] %m+% lubridate::weeks(1:n.ahead)
        }else if(model$interpolate == "daily"){
          dates.fc = dates[length(dates)] %m+% lubridate::days(1:n.ahead)
        }
      }else if(floor(freq) == 1){
        #Yearly frequency
        dates.fc = dates[length(dates)] %m+% lubridate::years(1:n.ahead)
      }
    }else{
      #Non-standard frequencies
      dates.dt = data.table(date = dates)[, "diff" := difftime(date, shift(date, type = "lag", n = 1), units = "days")]
      interval = as.numeric(mean(dates.dt$diff, na.rm = T))
      dates.fc = dates[length(dates)] + as.difftime(seq(interval, as.numeric(interval*n.ahead), interval), units = "days")
    }
  }else{
    y.fc = NULL
    dates.fc = NULL
  }
  
  #Set the trend
  if(!"Tt_0" %in% colnames(series)){
    series[, "trend" := 0]
  }else{
    if(is.null(model$interpolate)){
      colnames(series)[colnames(series) == "Tt_0"] = "trend"
    }else if(model$interpolate_method == "eop"){
      series[, "trend" := shift(Tt_0, n = 1, type = "lag")]
    }else if(model$interpolate_method == "avg"){
      series[, "trend" := shift(frollmean(Tt_0, n = int_per, align = "right"), type = "lag", n = 1)]
    }else if(model$interpolate_method == "sum"){
      series[, "trend" := shift(frollsum(Tt_0, n = int_per, align = "right"), type = "lag", n = 1)]
    }
  }
  
  #Set the drift
  if(!"Dt_0" %in% colnames(series)){
    series[, "drift" := 0]
  }else{
    if(is.null(model$interpolate)){
      colnames(series)[colnames(series) == "Dt_0"] = "drift"
    }else if(model$interpolate_method == "eop"){
      series[, "drift" := shift(Dt_0, n = 1, type = "lag")]
    }else if(model$interpolate_method == "avg"){
      series[, "drift" := shift(frollmean(Dt_0, n = int_per, align = "right"), type = "lag", n = 1)]
    }else if(model$interpolate_method == "sum"){
      series[, "drift" := shift(frollsum(Dt_0, n = int_per, align = "right"), type = "lag", n = 1)]
    }
  }
  
  #Set the cycle
  if(!"Ct_0" %in% colnames(series)){
    series[, "cycle" := 0]
  }else{
    if(is.null(model$interpolate)){
      colnames(series)[colnames(series) == "Ct_0"] = "cycle"
    }else if(model$interpolate_method == "eop"){
      series[, "cycle" := shift(Ct_0, n = 1, type = "lag")]
    }else if(model$interpolate_method == "avg"){
      series[, "cycle" := shift(frollmean(Ct_0, n = int_per, align = "right"), type = "lag", n = 1)]
    }else if(model$interpolate_method == "sum"){
      series[, "cycle" := shift(frollsum(Ct_0, n = int_per, align = "right"), type = "lag", n = 1)]
    } 
  }
  
  #Set the interpolation
  if("It_0" %in% colnames(series)){
    colnames(series)[colnames(series) == "It_0"] = "interpolated"
  }
  
  #Set the seasonal
  if(!any(grepl("St", colnames(series)))){
    series[, "seasonal" := 0]
  }else{
    colnames(series)[grepl("St", colnames(series)) & !grepl("Sts", colnames(series))] = 
      gsub("\\_0", "", gsub("St", "seasonal", colnames(series)[grepl("St", colnames(series)) & !grepl("Sts", colnames(series))]))
    if(sum(grepl("seasonal", colnames(series))) > 1){
      series[, "seasonal" := rowSums(series[, grepl("seasonal", colnames(series)), with = FALSE])]
    }else{
      colnames(series)[grepl("seasonal", colnames(series))] = "seasonal"
    }
    if(!is.null(model$interpolate)){
      if(model$interpolate_method == "eop"){
        series[, "seasonal" := shift(seasonal, n = 1, type = "lag")]
      }else if(model$interpolate_method == "avg"){
        series[, "seasonal" := shift(frollmean(seasonal, n = int_per, align = "right"), type = "lag", n = 1)]
      }else if(model$interpolate_method == "sum"){
        series[, "seasonal" := shift(frollsum(seasonal, n = int_per, align = "right"), type = "lag", n = 1)]
      } 
    }
  }
  
  #Set the exogenous
  if(!is.null(exo)){
    XB = do.call("cbind", lapply(1:nrow(X), function(x){
      name = names(par)[grepl("beta_", names(par))][x]
      par[name]*X[x, ]
    }))
    colnames(XB) = rownames(X)
    series = cbind(series, XB)
    rm(XB)
  }
  
  #Get the remainder
  series[, "remainder" := c(y, y.fc) - fitted]
  
  #Combine the filtered series
  final = data.table(date = c(dates, dates.fc), observed = c(y, y.fc), 
                     series[, colnames(series) %in% c("trend", "drift", "cycle", "remainder", "fitted", "interpolated", rownames(X)) | grepl("seasonal|\\%", colnames(series)), with = FALSE])
  rm(series, B_tt)
  
  #Calculate adjusted series
  final[, "seasonal_adjusted" := observed - seasonal]
  final[, "cycle_adjusted" := observed - cycle]
  final[, "seasonal_cycle_adjusted" := observed - seasonal - cycle]
  
  #Calculate noise adjusted series
  final[, "seasonal_noise_adjusted" := observed - seasonal - remainder]
  final[, "cycle_noise_adjusted" := observed - cycle - remainder]
  final[, "seasonal_cycle_noise_adjusted" := observed - cycle - seasonal - remainder]
  final[, "noise_adjusted" := observed - remainder]
  
  cols = colnames(final)[colnames(final) %in% c("observed", "trend", "fitted", paste0((1/2 + c(-ci, ci)/2)*100, "%")) | grepl("_adjusted", colnames(final))]
  if(!is.null(model$interpolate)){
    if(model$interpolate_method %in% c("avg", "eop")){
      cols = c(cols, "interpolated")
    }
  }
  final[, c(cols) := lapply(.SD, function(x){x*sd + mean}), .SDcols = c(cols)]
  
  cols = colnames(final)[!colnames(final) %in% c("date", cols, rownames(X))]
  final[,  c(cols) := lapply(.SD, function(x){x*sd}), .SDcols = c(cols)]
  if(!is.null(model$interpolate)){
    if(model$interpolate_method == "sum"){
      final[, "interpolated" := interpolated + mean/int_per]
    }
  }
  
  if(model$multiplicative == TRUE){
    final[, c(colnames(final)[!colnames(final) %in% c("date")]) := lapply(.SD, exp),
          .SDcols = c(colnames(final)[!colnames(final) %in% c("date")])]
  }
  final[date %in% dates, "forecast" := FALSE]
  final[is.na(forecast), "forecast" := TRUE]
  final[forecast == TRUE, "observed" := NA]
  
  if(plot.decomp == TRUE){
    g = list()
    g[["trend"]] = ggplot(data.table::melt(final[forecast == FALSE, ], id.vars = "date", measure.vars = c("observed", "trend"))[!is.na(value), ]) +
      labs(title = "Observed vs Trend", x = "", y = "") +
      geom_line(aes(x = date, y = value, group = variable, color = variable)) + theme_minimal() +
      scale_color_viridis_d() +
      guides(color = guide_legend(title = NULL), size = guide_legend(title = NULL)) +
      theme(legend.position = "bottom")
    
    g[["cycle"]] = ggplot(data.table::melt(final[forecast == FALSE, ], id.vars = "date", measure.vars = c("cycle"))) +
      labs(title = "Cycle", x = "", y = "") +
      geom_line(aes(x = date, y = value, group = variable, color = variable)) +
      geom_hline(aes(yintercept = ifelse(model$multiplicative == TRUE, 1, 0)), color = "black") + 
      scale_color_viridis_d() +
      theme_minimal() + guides(color = guide_legend(title = NULL)) +
      theme(legend.position = "bottom")
    
    if(grepl("seasonal", model$decomp)){
      seas_plot = data.table::melt(final[forecast == FALSE, ], id.vars = "date", measure.vars = colnames(final)[colnames(final) %in% c("seasonal", paste0("seasonal", seasons))])
      seas_plot[, "variable" := factor(variable, levels = c("seasonal", rev(sort(paste0("seasonal", seasons)))), ordered = TRUE)]
      g[["seasonal"]] = ggplot(seas_plot) +
        labs(title = "Seasonal", x = "", y = "") +
        geom_line(aes(x = date, y = value, group = variable, color = variable)) +
        geom_hline(aes(yintercept = ifelse(model$multiplicative == TRUE, 1, 0)), color = "black") + 
        scale_color_viridis_d() +
        theme_minimal() + guides(color = guide_legend(title = NULL)) +
        theme(legend.position = "bottom")
      rm(seas_plot)
    }
    
    if(length(unique(final$remainder)) > 1){
      g[["remainder"]] = ggplot(data.table::melt(final[forecast == FALSE, ], id.vars = "date", measure.vars = c("remainder"))) +
        labs(title = "Remainder", x = "", y = "") +
        geom_line(aes(x = date, y = value, group = variable, color = variable)) +
        geom_hline(aes(yintercept = ifelse(model$multiplicative == TRUE, 1, 0)), color = "black") + 
        scale_color_viridis_d() +
        theme_minimal() + guides(color = guide_legend(title = NULL)) +
        theme(legend.position = "bottom")
    }
    
    layout_matrix = 1:(2*ceiling(length(g)/2))
    if(length(layout_matrix) > length(g)){
      layout_matrix = c(layout_matrix[1], 1:(length(layout_matrix) - 1))
    }
    layout_matrix = matrix(layout_matrix, ncol = 2, byrow = TRUE)
    suppressWarnings(gridExtra::grid.arrange(grobs = g, layout_matrix = layout_matrix))
  }
  
  if(plot.fc == TRUE){
    fc_plot = data.table::melt(final[(.N - (n.hist + n.ahead) + 1):.N, ], id.vars = c("date", "forecast"), measure.vars = c("observed", "trend", paste0((1/2 + c(-ci, ci)/2)*100, "%")))
    fc_plot[variable != "trend" & forecast == FALSE, "group" := "observed"]
    fc_plot[variable != "trend" & forecast == TRUE, "group" := "forecast"]
    fc_plot[variable == "trend", "group" := "trend"]
    fc_plot[variable %in% paste0((1/2 + c(-ci, ci)/2)*100, "%"), "group" := "conf. int."]
    fc_plot[, "group" := factor(group, levels = c("observed", "forecast", "trend", "conf. int."), ordered = TRUE)]
    
    suppressWarnings(plot(ggplot(fc_plot) +
                            labs(title = "Observed vs Forecast", subtitle = paste0(ci*100, "% Confidence Interval"), x = "", y = "") +
                            geom_ribbon(data = dcast(fc_plot, "date ~ variable", value.var = "value"),
                                        aes(x = date,
                                            ymin = eval(parse(text = paste0("`", (1/2 - ci/2)*100, "%", "`"))),
                                            ymax = eval(parse(text = paste0("`", (1/2 + ci/2)*100, "%", "`")))), alpha = 0.5) +
                            geom_line(aes(x = date, y = value, group = variable, color = group)) + theme_minimal() +
                            scale_color_viridis_d() +
                            guides(color = guide_legend(title = NULL)) +
                            theme(legend.position = "bottom")))
  }
  colnames(final)[grepl("%", colnames(final))] = paste0(colnames(final)[grepl("%", colnames(final))], "_fc")
  return(final)
}

#' Kalman Filter
#'
#' Kalman filter an estimated model from stsm_estimate output. This is a wrapper to stsm_forecast with n.ahead = 0.
#' @param model Structural time series model estimated using stsm_estimate.
#' @param ci Confidence interval, value between 0 and 1 exclusive.
#' @param y Univariate time series of data values. May also be a 2 column data frame containing a date column.
#' @param freq Frequency of the data (1 (yearly), 4 (quarterly), 12 (monthly), 365.25/7 (weekly), 365.25 (daily)), default is NULL and will be automatically detected
#' @param exo Matrix of exogenous variables used for the historical data. Can be used to specify regression effects or other seasonal effects like holidays, etc.
#' @param plot, Logical, whether to plot everything
#' @param plot.decomp Logical, whether to plot the filtered historical data
#' @param smooth Whether or not to use the Kalman smoother
#' @import data.table ggplot2
#' @useDynLib autostsm, .registration=TRUE
#' @return data table (or list of data tables) containing the filtered and/or smoothed series.
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
#' fc = stsm_filter(stsm, y = NA000334Q, plot = TRUE)
#' }
#' @export
stsm_filter = function(model, y, freq = NULL, exo = NULL, ci = 0.8, 
                       plot = FALSE, plot.decomp = FALSE, smooth = TRUE){
  stsm = stsm_forecast(model, y, n.ahead = 0, freq = NULL, exo = NULL, exo.fc = NULL, ci = 0.8, 
                plot = plot, plot.decomp = plot.decomp, plot.fc = FALSE, n.hist = NULL, smooth = TRUE, dampen_cycle = FALSE)
  return(stsm)
}
