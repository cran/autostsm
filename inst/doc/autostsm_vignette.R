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

