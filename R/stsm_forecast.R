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
#' @param n.ahead the number of periods to forecast
#' @param dampen_cycle Whether to remove oscillating cycle dynamics and smooth the cycle forecast into the trend using a sigmoid function that maintains the rate of convergence
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
                         plot = FALSE, plot.decomp = FALSE, plot.fc = FALSE, n.hist = NULL, smooth = TRUE, dampen_cycle = FALSE){
  #model = stsm
  #n.ahead = floor(model$freq*3)
  #ci = 0.8
  #freq = exo = exo.fc = n.hist = NULL
  #plot = plot.decomp = plot.fc = smooth = TRUE
  #dampen_cycle = FALSE
  if(!is.logical(plot)){
    stop("plot.decomp must be TRUE or FALSE")
  }
  if(!is.logical(plot.decomp)){
    stop("plot.decomp must be TRUE or FALSE")
  }
  if(!is.logical(plot.fc)){
    stop("plot.fc must be TRUE or FALSE")
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
  }else{
    freq = y$freq
  }
  y = y$data
  
  #Remove leading and trailing NAs
  range = which(!is.na(y))
  y = unname(y[range[1]:range[length(range)]])
  dates = dates[range[1]:range[length(range)]]
  if(!is.null(exo)){
    exo = exo[range[1]:range[length(range)], ]
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
      if(ncol(exo.fc) != n.ahead){
        stop("nrow of exo.fc does not equal n.ahead.")
      }
      X = cbind(X, exo.fc)
      if(ncol(X) != length(y) + n.ahead){
        stop("nrow of exo + exo.fc does not equal length of y + n.ahead.")
      }
    }
  }else{
    if(n.ahead > 0){
      X = cbind(X, matrix(0, ncol = n.ahead))
    }
  }
  
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
  
  #Get model specifications
  decomp = model$decomp
  trend = model$trend
  if(!is.na(model$harmonics)){
    harmonics = as.numeric(strsplit(model$harmonics, ", ")[[1]])
  }else{
    harmonics = c()
  }
  
  #Get the coefficients
  par = eval(parse(text = paste0("c(", model$coef, ")")))
  init = stsm_init_vals(y, par, freq, trend, decomp, harmonics)
  
  #Filter and smooth the data
  sp = stsm_ssm(par = par, yt = y, freq = freq, decomp = decomp, trend = trend, init = init)
  B_tt = kalman_filter(B0 = matrix(sp$B0, ncol = 1), P0 = sp$P0, Dt = sp$Dt, At = sp$At, Ft = sp$Ft, Ht = sp$Ht, Qt = sp$Qt, Rt = sp$Rt,
                       yt = matrix(y, nrow = 1), X = X, beta = sp$beta)[c("B_tl", "B_tt", "P_tt", "P_tl")]
  if(smooth == TRUE){
    B_tt = kalman_smoother(B_tl = B_tt$B_tl, B_tt = B_tt$B_tt, P_tl = B_tt$P_tl, P_tt = B_tt$P_tt, Ft = sp$Ft)$B_tt
  }else{
    B_tt = B_tt$B_tt
  }
  rownames(B_tt) = rownames(sp$Ft)

  #Get the unobserved series
  series = data.table(t(B_tt))
  series$fitted = c(matrix(sp$At, nrow = 1, ncol = nrow(series)) + 
                      sp$Ht %*% t(as.matrix(series[, rownames(B_tt), with = FALSE])) + 
                      matrix(par[grepl("beta_", names(par))], nrow = 1) %*% X[, 1:length(y)])

  #Forecast
  fev = extrema = fitted = variable = value = cycle = seasonal = observed = remainder = 
    forecast = group = NULL
  if(n.ahead > 0){
    Fm_pow = diag(nrow(sp$Ft))
    for(j in 1:n.ahead){
      Fm_pow = Fm_pow %*% sp$Ft
      series = rbind(series, cbind(t(sp$Dt + sp$Ft %*% t(as.matrix(series[.N, rownames(B_tt), with = FALSE]))),
                                   c(sp$Ht %*% (Fm_pow %*% sp$Qt %*% t(Fm_pow)) %*% t(sp$Ht))),
                     use.names = TRUE, fill = TRUE)
    }
    colnames(series)[ncol(series)] = "fev"
    series[!is.na(fev), "fev" := cumsum(fev) + c(sp$Rt)]
    
    if(dampen_cycle == TRUE & is.complex(eigen(sp$Ft[grepl("Ct|et", colnames(sp$Ft)), grepl("Ct|et", rownames(sp$Ft))])$values) & n.ahead > 1){
      #Smooth the cycle forecast into the trend using a sigmoid function so as not to get oscillating cycle forecasts
      Fm = sp$Ft[grepl("Ct|et", colnames(sp$Ft)), grepl("Ct|et", rownames(sp$Ft))]
      B = matrix(unlist(series[length(y), colnames(series)[grepl("Ct|et", colnames(sp$Ft))], with = FALSE]), ncol = 1)
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
        series[(length(y) + wh1):.N, "Ct" := z]
      }
    }
    
    series[!is.na(fev), "fitted" := c(matrix(sp$At, nrow = 1, ncol = nrow(series[!is.na(fev), ])) +  + 
                                        sp$Ht %*% t(as.matrix(series[!is.na(fev), rownames(B_tt), with = FALSE])) + 
                                        matrix(par[grepl("beta_", names(par))], nrow = 1) %*% X[, (length(y) + 1):ncol(X)])]
    
    series[!is.na(fev), paste0((1/2 - ci/2)*100, "%") := stats::qnorm(1/2 - ci/2, fitted, sqrt(fev))]
    series[!is.na(fev), paste0((1/2 + ci/2)*100, "%") := stats::qnorm(1/2 + ci/2, fitted, sqrt(fev))]
    series[length(y), c(paste0((1/2 - ci/2)*100, "%"), paste0((1/2 + ci/2)*100, "%")) := fitted]
    
    #Create envelope around seasonal fluctuations in the confidence intervals
    if(!is.null(harmonics) & all(!is.na(harmonics)) & length(harmonics) > 0){
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
    series[, "fev" := NULL]
    rm(Fm_pow)
  }
  
  #Set the forecast dates
  `%m+%` = lubridate::`%m+%`
  if(n.ahead > 0){
    y.fc = c(sp$At[1, ] + sp$Ht %*% t(as.matrix(series[(length(y) + 1):.N, rownames(B_tt), with = FALSE])))
    if(floor(freq) == 365){
      dates.fc = dates[length(dates)] %m+% lubridate::days(1:n.ahead)
    }else if(floor(freq) == 260){
      dates.fc = dates[length(dates)] %m+% lubridate::days(1:ceiling(n.ahead*2))
      dates.fc = dates.fc[which(!weekdays(dates.fc) %in% c("Saturday", "Sunday"))][1:n.ahead]
    }else if(floor(freq) == 52){
      dates.fc = dates[length(dates)] %m+% lubridate::weeks(1:n.ahead)
    }else if(floor(freq) == 12){
      dates.fc = dates[length(dates)] %m+% months(1:n.ahead)
    }else if(floor(freq) == 4){
      dates.fc = dates[length(dates)] %m+% months((1:n.ahead)*3)
    }else if(floor(freq) == 1){
      dates.fc = dates[length(dates)] %m+% lubridate::years(1:n.ahead)
    }
  }else{
    y.fc = NULL
    dates.fc = NULL
  }
  
  #Set the trend
  if(!"Tt0" %in% colnames(series)){
    series[, "trend" := 0]
  }else{
    colnames(series)[colnames(series) == "Tt0"] = "trend"
  }
  
  #Set the drift
  if(!"Mt0" %in% colnames(series)){
    series[, "drift" := 0]
  }else{
    colnames(series)[colnames(series) == "Mt0"] = "drift"
  }
  
  #Set the cycle
  if(!"Ct" %in% colnames(series)){
    series[, "Ct" := 0]
  }else{
    colnames(series)[colnames(series) == "Ct"] = "cycle"
  }
  
  #Set the seasonal
  if(!any(grepl("St", colnames(series)))){
    series[, "seasonal" := 0]
  }else{
    colnames(series)[grepl("St", colnames(series)) & !grepl("Sts", colnames(series))] = 
      gsub("St", "seasonal", colnames(series)[grepl("St", colnames(series)) & !grepl("Sts", colnames(series))])
    if(sum(grepl("seasonal", colnames(series))) > 1){
      series[, "seasonal" := rowSums(series[, grepl("seasonal", colnames(series)), with = FALSE])]
    }else{
      colnames(series)[grepl("seasonal", colnames(series))] = "seasonal"
    }
  }
  
  #Set the exogenous
  if(!is.null(exo)){
    XB = do.call("cbind", lapply(names(par)[grepl("beta_", names(par))], function(x){
      par[x]*X[gsub("beta_", "", x), ]
    }))
    colnames(XB) = rownames(X)
    series = cbind(series, XB)
    rm(XB)
  }
  
  #Get the remainder
  series[, "remainder" := c(y, y.fc) - trend - cycle - seasonal]
  
  #Combine the filtered series
  final = data.table(date = c(dates, dates.fc), observed = c(y, y.fc), 
                     series[, colnames(series) %in% c("trend", "drift", "cycle", "remainder", "fitted", rownames(X)) | grepl("seasonal|\\%", colnames(series)), with = FALSE])
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
  final[, c(cols) := lapply(.SD, function(x){x*sd + mean}), .SDcols = c(cols)]
  
  cols = colnames(final)[!colnames(final) %in% c("date", cols)]
  final[,  c(cols) := lapply(.SD, function(x){x*sd}), .SDcols = c(cols)]
  
  if(model$multiplicative == TRUE){
    final[, c(colnames(final)[!colnames(final) %in% c("date")]) := lapply(.SD, exp),
          .SDcols = c(colnames(final)[!colnames(final) %in% c("date")])]
  }
  final[date %in% dates, "forecast" := FALSE]
  final[is.na(forecast), "forecast" := TRUE]
  
  if(plot.decomp == TRUE){
    g = list()
    g[["trend"]] = ggplot(data.table::melt(final[forecast == FALSE, ], id.vars = "date", measure.vars = c("observed", "trend"))) +
      labs(title = "Observed vs Trend", x = "", y = "") +
      geom_line(aes(x = date, y = value, group = variable, color = variable)) + theme_minimal() +
      scale_color_viridis_d() +
      guides(color = guide_legend(title = NULL), size = guide_legend(title = NULL)) +
      theme(legend.position = "bottom")
    
    g[["cycle"]] = ggplot(data.table::melt(final[forecast == FALSE, ], id.vars = "date", measure.vars = c("cycle"))) +
      labs(title = "Cycle", x = "", y = "") +
      geom_line(aes(x = date, y = value, group = variable, color = variable)) +
      scale_color_viridis_d() +
      theme_minimal() + guides(color = guide_legend(title = NULL)) +
      theme(legend.position = "bottom")
    
    if(grepl("seasonal", decomp)){
      seas_plot = data.table::melt(final[forecast == FALSE, ], id.vars = "date", measure.vars = colnames(final)[colnames(final) %in% c("seasonal", paste0("seasonal", round(harmonics)))])
      seas_plot[, "variable" := factor(variable, levels = c("seasonal", rev(sort(paste0("seasonal", round(harmonics))))), ordered = TRUE)]
      g[["seasonal"]] = ggplot(seas_plot) +
        labs(title = "Seasonal", x = "", y = "") +
        geom_line(aes(x = date, y = value, group = variable, color = variable)) +
        scale_color_viridis_d() +
        theme_minimal() + guides(color = guide_legend(title = NULL)) +
        theme(legend.position = "bottom")
      rm(seas_plot)
    }
    
    if(length(unique(final$remainder)) > 1){
      g[["remainder"]] = ggplot(data.table::melt(final[forecast == FALSE, ], id.vars = "date", measure.vars = c("remainder"))) +
        labs(title = "Remainder", x = "", y = "") +
        geom_line(aes(x = date, y = value, group = variable, color = variable)) +
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