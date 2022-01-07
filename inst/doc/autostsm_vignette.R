## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  echo = TRUE, include = TRUE, message = FALSE, warning = FALSE, eval = FALSE
)

## -----------------------------------------------------------------------------
#  library(ggplot2)
#  library(gridExtra)
#  library(autostsm)
#  library(lubridate)
#  
#  set.seed(1024)
#  
#  #Daily data
#  freq = 365.25
#  
#  #Build the trend and drift
#  t = c()
#  m = c()
#  t[1] = 100
#  m[1] = 1
#  sig_e = 0.1
#  sig_t = 1
#  sig_m = 0.1
#  sig_s = 0.01
#  for(i in 2:3000){
#    m[i] = 0.05 + 0.75*m[i-1] + rnorm(1, 0, sig_m)
#    t[i] = t[i-1] + m[i-1] + rnorm(1, 0, sig_t)
#  }
#  
#  #Build the seasonality including yearly and weekly
#  s365 = sin(2*pi/freq*(1:length(t))) + rnorm(length(t), 0, sig_s)
#  s365 = (s365 - min(s365))/diff(range(s365))*(1.125 - 0.865) + 0.865
#  s7 = sin(2*pi/7*(1:length(t))) + rnorm(length(t), 0, sig_s)
#  s7 = (s7 - min(s7))/diff(range(s7))*(1.125 - 0.865) + 0.865
#  s = s365 + s7
#  s = (s - min(s))/diff(range(s))*(1.25 - 0.75) + 0.75
#  
#  #Build the cyclicality using every 3 years periodicity
#  c = sin(2*pi*0.33/freq*(1:length(t))) + rnorm(length(t), 0, sig_s)
#  c = (c - min(c))/diff(range(c))*(1.25 - 0.75) + 0.75
#  
#  #Build the data using a multiplicative model
#  ts = data.table(date = as.Date("2016-01-01") + 1:length(t),
#                  y = t*c*s*exp(rnorm(length(t), 0, sig_e)),
#                  trend = t, seasonal = s, seasonal7 = s7,
#                  seasonal365 = s365, cycle = c)
#  
#  #Create some missing values
#  ts[sample(1:nrow(ts), round(0.05*nrow(ts))), "y" := NA]
#  
#  #View the data
#  g1 = ggplot(melt(ts, id.vars = "date", measure.vars = c("y", "trend"))) +
#    labs(title = "Observed vs Trend") +
#    geom_line(aes(x = date, y = value, group = variable, color = variable)) +
#    scale_color_viridis_d() +
#    theme_minimal() + guides(color = guide_legend(title = NULL)) +
#    theme(legend.position = "bottom")
#  g2 = ggplot(melt(ts, id.vars = "date", measure.vars = c("cycle"))) +
#    labs(title = "Cycle") +
#    geom_line(aes(x = date, y = value, group = variable, color = variable)) +
#    scale_color_viridis_d() +
#    theme_minimal() + guides(color = guide_legend(title = NULL)) +
#    theme(legend.position = "bottom")
#  g3 = ggplot(melt(ts, id.vars = "date", measure.vars = colnames(ts)[grepl("seasonal", colnames(ts))])) +
#    labs(title = "Seasonal") +
#    geom_line(aes(x = date, y = value, group = variable, color = variable)) +
#    scale_color_viridis_d() +
#    theme_minimal() + guides(color = guide_legend(title = NULL)) +
#    theme(legend.position = "bottom")
#  grid.arrange(g1, g2, g3, layout_matrix = rbind(c(1, 1), c(2, 3)))
#  
#  #Estimate the model
#  stsm = stsm_estimate(ts[, c("date", "y"), with = FALSE], verbose = TRUE)
#  
#  #Forecast and plot the results
#  stsm_fc = stsm_forecast(stsm, y = ts[, c("date", "y"), with = FALSE], n.ahead = floor(stsm$freq)*3, plot = TRUE)
#  
#  #Detect anomalies
#  stsm_fc = merge(stsm_fc,
#                  stsm_detect_anomalies(stsm, y = ts[, c("date", "y"), with = FALSE], plot = TRUE),
#                  by = "date", all = TRUE)
#  
#  #Detect structural breaks
#  stsm_fc = merge(stsm_fc,
#                  stsm_detect_breaks(stsm, y = ts[, c("date", "y"), with = FALSE], plot = TRUE, show_progress = TRUE),
#                  by = "date", all = TRUE)

## -----------------------------------------------------------------------------
#  library(autostsm)
#  
#  ##### Unemployment rate examples #####
#  #Not seasonally adjusted
#  data("UNRATENSA", package = "autostsm") #From FRED
#  UNRATENSA = data.table(UNRATENSA, keep.rownames = TRUE)
#  colnames(UNRATENSA) = c("date", "y")
#  UNRATENSA[, "date" := as.Date(date)]
#  UNRATENSA[, "y" := as.numeric(y)]
#  stsm = stsm_estimate(UNRATENSA, verbose = TRUE)
#  stsm_fc = stsm_forecast(stsm, y = UNRATENSA, n.ahead = floor(stsm$freq)*3, plot = TRUE)
#  stsm_fc = merge(stsm_fc,
#                  stsm_detect_anomalies(stsm, y = UNRATENSA, plot = TRUE),
#                  by = "date", all = TRUE)
#  stsm_fc = merge(stsm_fc,
#                  stsm_detect_breaks(stsm, y = UNRATENSA, plot = TRUE, show_progress = TRUE),
#                  by = "date", all = TRUE)
#  
#  #Seasonally adjusted
#  data("UNRATE", package = "autostsm") #From FRED
#  UNRATE = data.table(UNRATE, keep.rownames = TRUE)
#  colnames(UNRATE) = c("date", "y")
#  UNRATE[, "date" := as.Date(date)]
#  UNRATE[, "y" := as.numeric(y)]
#  stsm = stsm_estimate(UNRATE, verbose = TRUE)
#  stsm_fc = stsm_forecast(stsm, y = UNRATE, n.ahead = floor(stsm$freq)*3, plot = TRUE)
#  stsm_fc = merge(stsm_fc,
#                  stsm_detect_anomalies(stsm, y = UNRATE, plot = TRUE),
#                  by = "date", all = TRUE)
#  stsm_fc = merge(stsm_fc,
#                  stsm_detect_breaks(stsm, y = UNRATE, plot = TRUE, show_progress = TRUE),
#                  by = "date", all = TRUE)
#  
#  ##### GDP examples #####
#  #Not seasonally adjusted
#  data("NA000334Q", package = "autostsm") #From FRED
#  NA000334Q = data.table(NA000334Q, keep.rownames = TRUE)
#  colnames(NA000334Q) = c("date", "y")
#  NA000334Q[, "date" := as.Date(date)]
#  NA000334Q[, "y" := as.numeric(y)]
#  NA000334Q = NA000334Q[date >= "1990-01-01", ]
#  stsm = stsm_estimate(NA000334Q, verbose = TRUE)
#  stsm_fc = stsm_forecast(stsm, y = NA000334Q, n.ahead = floor(stsm$freq)*3, plot = TRUE)
#  stsm_fc = merge(stsm_fc,
#                  stsm_detect_anomalies(stsm, y = NA000334Q, plot = TRUE),
#                  by = "date", all = TRUE)
#  stsm_fc = merge(stsm_fc,
#                  stsm_detect_breaks(stsm, y = NA000334Q, plot = TRUE, show_progress = TRUE),
#                  by = "date", all = TRUE)
#  
#  #Seasonally adjusted
#  data("GDP", package = "autostsm") #From FRED
#  GDP = data.table(GDP, keep.rownames = TRUE)
#  colnames(GDP) = c("date", "y")
#  GDP[, "date" := as.Date(date)]
#  GDP[, "y" := as.numeric(y)]
#  GDP = GDP[date >= "1990-01-01", ]
#  stsm = stsm_estimate(GDP, verbose = TRUE)
#  stsm_fc = stsm_forecast(stsm, y = GDP, n.ahead = floor(stsm$freq)*3, plot = TRUE)
#  stsm_fc = merge(stsm_fc,
#                  stsm_detect_anomalies(stsm, y = GDP, plot = TRUE),
#                  by = "date", all = TRUE)
#  stsm_fc = merge(stsm_fc,
#                  stsm_detect_breaks(stsm, y = GDP, plot = TRUE, show_progress = TRUE),
#                  by = "date", all = TRUE)
#  
#  ##### S&P 500 example #####
#  data("SP500", package = "autostsm") #From FRED
#  SP500 = data.table(SP500, keep.rownames = TRUE)
#  colnames(SP500) = c("date", "y")
#  SP500[, "date" := as.Date(date)]
#  SP500[, "y" := as.numeric(y)]
#  stsm = stsm_estimate(SP500, verbose = TRUE)
#  stsm_fc = stsm_forecast(stsm, y = SP500, n.ahead = floor(stsm$freq)*3, plot = TRUE)
#  stsm_fc = merge(stsm_fc,
#                  stsm_detect_anomalies(stsm, y = SP500, plot = TRUE),
#                  by = "date", all = TRUE)
#  stsm_fc = merge(stsm_fc,
#                  stsm_detect_breaks(stsm, y = SP500, plot = TRUE, show_progress = TRUE),
#                  by = "date", all = TRUE)
#  
#  ##### 5 Year Treasury Yield example #####
#  data("DGS5", package = "autostsm") #From FRED
#  DGS5 = data.table(DGS5, keep.rownames = TRUE)
#  colnames(DGS5) = c("date", "y")
#  DGS5[, "date" := as.Date(date)]
#  DGS5[, "y" := as.numeric(y)]
#  stsm = stsm_estimate(DGS5, verbose = TRUE)
#  stsm_fc = stsm_forecast(stsm, y = DGS5, n.ahead = floor(stsm$freq)*3, plot = TRUE)
#  stsm_fc = merge(stsm_fc,
#                  stsm_detect_anomalies(stsm, y = DGS5, plot = TRUE),
#                  by = "date", all = TRUE)
#  stsm_fc = merge(stsm_fc,
#                  stsm_detect_breaks(stsm, y = DGS5, plot = TRUE, show_progress = TRUE),
#                  by = "date", all = TRUE)

## -----------------------------------------------------------------------------
#  library(ggplot2)
#  library(autostsm)
#  library(lubridate)
#  
#  #Rebuild the data
#  set.seed(1024)
#  
#  #Daily data
#  freq = 365.25
#  
#  #Build the trend and drift
#  t = c()
#  m = c()
#  t[1] = 100
#  m[1] = 1
#  sig_e = 0.1
#  sig_t = 1
#  sig_m = 0.1
#  sig_s = 0.01
#  for(i in 2:3000){
#    m[i] = 0.05 + 0.75*m[i-1] + rnorm(1, 0, sig_m)
#    t[i] = t[i-1] + m[i-1] + rnorm(1, 0, sig_t)
#  }
#  
#  #Build the seasonality including yearly and weekly
#  s365 = sin(2*pi/freq*(1:length(t))) + rnorm(length(t), 0, sig_s)
#  s365 = (s365 - min(s365))/diff(range(s365))*(1.125 - 0.865) + 0.865
#  s7 = sin(2*pi/7*(1:length(t))) + rnorm(length(t), 0, sig_s)
#  s7 = (s7 - min(s7))/diff(range(s7))*(1.125 - 0.865) + 0.865
#  s = s365 + s7
#  s = (s - min(s))/diff(range(s))*(1.25 - 0.75) + 0.75
#  
#  #Build the cyclicality using every 3 years periodicity
#  c = sin(2*pi*0.33/freq*(1:length(t))) + rnorm(length(t), 0, sig_s)
#  c = (c - min(c))/diff(range(c))*(1.25 - 0.75) + 0.75
#  
#  #Build the data using a multiplicative model
#  ts = data.table(date = as.Date("2016-01-01") + 1:length(t),
#                  y = t*c*s*exp(rnorm(length(t), 0, sig_e)))
#  
#  #Convert the data to monthly by using the end of period
#  ts[, "mnth" := floor_date(date, "month")]
#  ts[, "rn" := .N:1, by = "mnth"]
#  ts = ts[rn == 1, c("mnth", "y"), with = FALSE]
#  colnames(ts)[colnames(ts) == "mnth"] = "date"
#  
#  #Get the quarter and set the date to be the end of the quarter
#  ts[, "qtr" := ceiling_date(date, "quarter") %m-% months(1)]
#  ts[, "y_avg" := frollmean(y, n = 3, align = "right")]
#  ts[, "y_sum" := frollsum(y, n = 3, align = "right")]
#  
#  #We already know the data has a multiplicative structure with yearly seasonality and
#  #the frequency below will be set to quarterly,
#  #so we will set seasons = 4 and multiplicative to TRUE to help the model with identification
#  #since the data set goes from 3000 daily observations to 99 monthly observations and then to
#  #33 quarterly observations
#  
#  #End of period
#  ts[, "rn" := .N:1, by = "qtr"]
#  ts.eop = ts[rn == 1, ][, "rn" := NULL]
#  stsm = stsm_estimate(y = ts.eop[, c("qtr", "y"), with = FALSE], seasons = 4, multiplicative = TRUE,
#                       interpolate = "monthly", interpolate_method = "eop", verbose = TRUE)
#  stsm_fc = stsm_filter(stsm, y = ts.eop[, c("qtr", "y"), with = FALSE], plot = TRUE)
#  stsm_fc = merge(stsm_fc[, c("date", "observed", "interpolated")],
#                  ts[, c("date", "y")], by = "date", all.x = TRUE, all.y = TRUE)
#  eop.err = cbind(method = "eop", stsm_fc[is.na(observed), .(mape = mean(abs(y - interpolated)/(1 + y)*100, na.rm = TRUE))])
#  ggplot(melt(stsm_fc, id.vars = "date", measure.vars = c("y", "interpolated"))[!is.na(value), ]) +
#    geom_line(aes(x = date, y = value, group = variable, color = variable)) +
#    scale_color_viridis_d() +
#    theme_minimal() + guides(color = guide_legend(title = NULL)) +
#    theme(legend.position = "bottom")
#  
#  #Period average
#  ts.avg = ts[, .(y = mean(y, na.rm = TRUE)), by = "qtr"]
#  stsm = stsm_estimate(y = ts.avg[, c("qtr", "y"), with = FALSE], seasons = 4, multiplicative = TRUE,
#                       interpolate = "monthly", interpolate_method = "avg", verbose = TRUE)
#  stsm_fc = stsm_filter(stsm, y = ts.avg[, c("qtr", "y"), with = FALSE], plot = TRUE)
#  stsm_fc = merge(stsm_fc[, c("date", "observed", "interpolated")],
#                  ts[, c("date", "y", "y_avg")], by = "date", all.x = TRUE, all.y = TRUE)
#  if(stsm$multiplicative == TRUE){
#    stsm_fc[, "interpolated_rollup" := exp(frollmean(log(interpolated), n = 3, align = "right"))]
#  }else{
#    stsm_fc[, "interpolated_rollup" := frollmean(interpolated, n = 3, align = "right")]
#  }
#  avg.err = cbind(method = "avg",
#                  stsm_fc[, .(mape = mean(abs(y - interpolated)/(1 + y)*100, na.rm = TRUE))],
#                  stsm_fc[is.na(observed), .(mape_rollup = mean(abs(y_avg - interpolated_rollup)/(1 + y_avg)*100, na.rm = TRUE))])
#  ggplot(melt(stsm_fc, id.vars = "date", measure.vars = c("y", "interpolated"))[!is.na(value), ]) +
#    geom_line(aes(x = date, y = value, group = variable, color = variable)) +
#    scale_color_viridis_d() +
#    theme_minimal() + guides(color = guide_legend(title = NULL)) +
#    theme(legend.position = "bottom")
#  
#  #Period sum
#  ts.sum = ts[, .(y = sum(y, na.rm = TRUE)), by = "qtr"]
#  stsm = stsm_estimate(y = ts.sum[, c("qtr", "y"), with = FALSE], seasons = 4, multiplicative = FALSE,
#                       interpolate = "monthly", interpolate_method = "sum", verbose = TRUE)
#  stsm_fc = stsm_filter(stsm, y = ts.sum[, c("qtr", "y"), with = FALSE], plot = TRUE)
#  stsm_fc = merge(stsm_fc[, c("date", "observed", "interpolated")],
#                  ts[, c("date", "y", "y_sum")], by = "date", all.x = TRUE, all.y = TRUE)
#  if(stsm$multiplicative == TRUE){
#    stsm_fc[, "interpolated_rollup" := exp(frollsum(log(interpolated), n = 3, align = "right"))]
#  }else{
#    stsm_fc[, "interpolated_rollup" := frollsum(interpolated, n = 3, align = "right")]
#  }
#  sum.err = cbind(method = "sum",
#                  stsm_fc[, .(mape = mean(abs(y - interpolated)/(1 + y)*100, na.rm = TRUE))],
#                  stsm_fc[is.na(observed), .(mape_rollup = mean(abs(y_sum - interpolated_rollup)/(1 + y_sum)*100, na.rm = TRUE))])
#  ggplot(melt(stsm_fc, id.vars = "date", measure.vars = c("y", "interpolated"))[!is.na(value), ]) +
#    geom_line(aes(x = date, y = value, group = variable, color = variable)) +
#    scale_color_viridis_d() +
#    theme_minimal() + guides(color = guide_legend(title = NULL)) +
#    theme(legend.position = "bottom")
#  
#  #Errors depend on the volatility of the series
#  #the end of period method is more susceptible to the volatility of the series than either
#  #the average or sum methods
#  rbind(eop.err, avg.err, sum.err, use.names = TRUE, fill = TRUE)

## -----------------------------------------------------------------------------
#  library(ggplot2)
#  library(autostsm)
#  library(lubridate)
#  
#  ##### Unemployment rate examples #####
#  #Not seasonally adjusted
#  data("UNRATENSA", package = "autostsm") #From FRED
#  UNRATENSA = data.table(UNRATENSA, keep.rownames = TRUE)
#  colnames(UNRATENSA) = c("date", "y")
#  UNRATENSA[, "date" := as.Date(date)]
#  UNRATENSA[, "y" := as.numeric(y)]
#  UNRATENSA[, "qtr" := ceiling_date(date, "quarter") %m-% months(1)]
#  UNRATENSA[, "y_avg" := frollmean(y, n = 3, align = "right")]
#  
#  #Cycle taken from previous example
#  stsm = stsm_estimate(y = UNRATENSA[, .(y = mean(y, na.rm = TRUE)), by = "qtr"],
#                       seasons = 4, cycle = 156/3,
#                       interpolate = "monthly", interpolate_method = "avg", verbose = TRUE)
#  stsm_fc = stsm_filter(stsm, y = UNRATENSA[, .(y = mean(y, na.rm = TRUE)), by = "qtr"], plot = TRUE)
#  stsm_fc = merge(stsm_fc[, c("date", "observed", "interpolated")],
#                  UNRATENSA[, c("date", "y", "y_avg")], by = "date", all.x = TRUE, all.y = TRUE)
#  if(stsm$multiplicative == TRUE){
#    stsm_fc[, "interpolated_rollup" := exp(frollmean(log(interpolated), n = 3, align = "right"))]
#  }else{
#    stsm_fc[, "interpolated_rollup" := frollmean(interpolated, n = 3, align = "right")]
#  }
#  cbind(stsm_fc[, .(mape = mean(abs(y - interpolated)/(1 + y)*100, na.rm = TRUE))],
#        stsm_fc[is.na(observed), .(mape_rollup = mean(abs(y_avg - interpolated_rollup)/(1 + y_avg)*100, na.rm = TRUE))])
#  ggplot(melt(stsm_fc, id.vars = "date", measure.vars = c("y", "interpolated"))[!is.na(value), ]) +
#    geom_line(aes(x = date, y = value, group = variable, color = variable)) +
#    scale_color_viridis_d() +
#    theme_minimal() + guides(color = guide_legend(title = NULL)) +
#    theme(legend.position = "bottom")
#  
#  #Seasonally adjusted
#  data("UNRATE", package = "autostsm") #From FRED
#  UNRATE = data.table(UNRATE, keep.rownames = TRUE)
#  colnames(UNRATE) = c("date", "y")
#  UNRATE[, "date" := as.Date(date)]
#  UNRATE[, "y" := as.numeric(y)]
#  UNRATE[, "qtr" := ceiling_date(date, "quarter") %m-% months(1)]
#  UNRATE[, "y_avg" := frollmean(y, n = 3, align = "right")]
#  
#  #Cycle taken from previous example
#  stsm = stsm_estimate(y = UNRATE[, .(y = mean(y, na.rm = TRUE)), by = "qtr"], cycle = 156/3,
#                       interpolate = "monthly", interpolate_method = "avg", verbose = TRUE)
#  stsm_fc = stsm_filter(stsm, y = UNRATE[, .(y = mean(y, na.rm = TRUE)), by = "qtr"], plot = TRUE)
#  stsm_fc = merge(stsm_fc[, c("date", "observed", "interpolated")],
#                  UNRATE[, c("date", "y", "y_avg")], by = "date", all.x = TRUE, all.y = TRUE)
#  if(stsm$multiplicative == TRUE){
#    stsm_fc[, "interpolated_rollup" := exp(frollmean(log(interpolated), n = 3, align = "right"))]
#  }else{
#    stsm_fc[, "interpolated_rollup" := frollmean(interpolated, n = 3, align = "right")]
#  }
#  cbind(stsm_fc[, .(mape = mean(abs(y - interpolated)/(1 + y)*100, na.rm = TRUE))],
#        stsm_fc[is.na(observed), .(mape_rollup = mean(abs(y_avg - interpolated_rollup)/(1 + y_avg)*100, na.rm = TRUE))])
#  ggplot(melt(stsm_fc, id.vars = "date", measure.vars = c("y", "interpolated"))[!is.na(value), ]) +
#    geom_line(aes(x = date, y = value, group = variable, color = variable)) +
#    scale_color_viridis_d() +
#    theme_minimal() + guides(color = guide_legend(title = NULL)) +
#    theme(legend.position = "bottom")
#  
#  #5 Year Treasury Yield
#  data("DGS5", package = "autostsm") #From FRED
#  DGS5 = data.table(DGS5, keep.rownames = TRUE)
#  colnames(DGS5) = c("date", "y")
#  DGS5[, "date" := as.Date(date)]
#  DGS5[, "y" := as.numeric(y)]
#  DGS5[, "qtr" := ceiling_date(date, "quarter") %m-% months(1)]
#  DGS5[, "y_avg" := frollmean(y, n = 3, align = "right")]
#  
#  stsm = stsm_estimate(y = DGS5[, .(y = mean(y, na.rm = TRUE)), by = "qtr"], seasons = FALSE,
#                       interpolate = "monthly", interpolate_method = "avg", verbose = TRUE)
#  stsm_fc = stsm_filter(stsm, y = DGS5[, .(y = mean(y, na.rm = TRUE)), by = "qtr"], plot = TRUE)
#  stsm_fc = merge(stsm_fc[, c("date", "observed", "interpolated")],
#                  DGS5[, c("date", "y", "y_avg")], by = "date", all.x = TRUE, all.y = TRUE)
#  if(stsm$multiplicative == TRUE){
#    stsm_fc[, "interpolated_rollup" := exp(frollmean(log(interpolated), n = 3, align = "right"))]
#  }else{
#    stsm_fc[, "interpolated_rollup" := frollmean(interpolated, n = 3, align = "right")]
#  }
#  cbind(stsm_fc[, .(mape = mean(abs(y - interpolated)/(1 + y)*100, na.rm = TRUE))],
#        stsm_fc[is.na(observed), .(mape_rollup = mean(abs(y_avg - interpolated_rollup)/(1 + y_avg)*100, na.rm = TRUE))])
#  ggplot(melt(stsm_fc, id.vars = "date", measure.vars = c("y", "interpolated"))[!is.na(value), ]) +
#    geom_line(aes(x = date, y = value, group = variable, color = variable)) +
#    scale_color_viridis_d() +
#    theme_minimal() + guides(color = guide_legend(title = NULL)) +
#    theme(legend.position = "bottom")

