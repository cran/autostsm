#' Detect Structural Breaks
#'
#' Detect structural breaks using the estimated structural time series model
#' @param model Structural time series model estimated using stsm_estimate.
#' @param components Vector of components to test for structural breaks
#' @param y Univariate time series of data values. May also be a 2 column data frame containing a date column.
#' @param freq Frequency of the data (1 (yearly), 4 (quarterly), 12 (monthly), 365.25/7 (weekly), 365.25 (daily)), default is NULL and will be automatically detected
#' @param exo Matrix of exogenous variables used for the historical data. Can be used to specify regression effects or other seasonal effects like holidays, etc.
#' @param plot Whether to plot everything
#' @param sig_level Significance level to determine statistically significant anomalies
#' @param ci Confidence interval, value between 0 and 1 exclusive.
#' @param smooth Whether or not to use the Kalman smoother
#' @param cores Number of cores to use for break detection
#' @return data table (or list of data tables) containing the dates of detected anomalies from the filtered and/or smoothed series
#' @import ggplot2 data.table foreach
#' @useDynLib autostsm, .registration=TRUE
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
#' breaks = stsm_detect_breaks(model = stsm, y = NA000334Q, plot = TRUE, cores = 2)
#' }
#' @export
stsm_detect_breaks = function(model, y, components = c("trend", "cycle", "seasonal"), freq = NULL, exo = NULL, 
                              sig_level = 0.01, ci = 0.8, smooth = TRUE, plot = FALSE, cores = NULL){
  #model = stsm
  #sig_level = 0.01
  #ci = 0.8
  #plot = smooth = TRUE
  #freq = exo = cores = NULL
  #components = c("trend", "cycle", "seasonal")
  if(!is.logical(plot)){
    stop("plot must be TRUE or FALSE")
  }
  if(sig_level < 0 | sig_level > 0.1){
    stop("sig_level must be positive and <= 0.1")
  }
  if(ci <= 0 | ci >= 1){
    stop("ci must between 0 and 1. Suggested values are 0.8, 0.9, 0.95, 0.99.")
  }
  if(any(!components %in% c("trend", "cycle", "seasonal"))){
    stop("components must be 'trend', 'cycle', or 'seasonal'")
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
  
  #Setup parallel computing
  cl = parallel::makeCluster(max(c(1, ifelse(is.null(cores), parallel::detectCores() - 1, 2))))
  doSNOW::registerDoSNOW(cl)
  
  #Get the unobserved series
  series = data.table(t(B_tt))
  
  #Detect structural breaks
  bp_dt = data.table()
  comps = c("Tt0", "Ct", "St")[c(c("Tt0", "Ct") %in% colnames(series), grepl("seasonal", decomp))]
  error = break_type = date.lag = date.lead = value = variable = lower = upper = Tt0_interval_lower = 
    Tt0_interval_upper = Ct_interval_lower = Ct_interval_upper = St_interval_lower = St_interval_upper = . = NULL
  for(j in comps[comps %in% c("Tt0", "Ct", "St")[c("trend", "cycle", "seasonal") %in% components]]){
    if(j == "Tt0"){
      dat = series[, c(j), with = FALSE][[1]]
      #Test if trend is trending up/down to test for breaks in growth, if not then breaks in level
      ol = forecast::tsoutliers(stats::ts(dat, frequency = freq))
      dat2 = copy(dat)
      dat2[ol$index] = ol$replacements
      if(tsutils::coxstuart(dat2, type = "trend")$p.value <= sig_level){
        nd = ifelse(model$trend %in% c("random-walk", "random-walk-drift"), 1, 
                    ifelse(model$trend %in% c("double-random-walk", "random-walk2"), 2, 0))
        dat_trans = diff(dat, differences = nd)
      }else{
        dat_trans = dat
      }
      h = round(freq) + 2
      bp = suppressWarnings(strucchange::breakpoints(dat_trans ~ 1, h = ifelse(h > 0.33*length(dat), round(0.33*length(dat)), h), 
                                                     alpha = sig_level, hpc = "foreach"))
    }else if(j == "St"){
      dat = series[, grepl("St", colnames(series)) & !grepl("Sts", colnames(series)), with = FALSE]
      dat_trans = dat = rowSums(dat)
      h = round(max(freq/harmonics)) + 2
      lm.data = data.table(y = dat_trans)[, "t" := 1:.N]
      for(k in floor(harmonics)){
        lm.data[, paste0("sin", k) := sin(2*pi*k/freq*t)]
        lm.data[, paste0("cos", k) := cos(2*pi*k/freq*t)]
      }
      lm.data[, "t" := NULL]
      lm_s = stats::lm(y ~ ., data = lm.data)
      lm.data = cbind(lm.data, do.call("cbind", lapply(floor(harmonics), function(x){
        matrix(stats::coef(lm_s)["(Intercept)"]/length(harmonics) + 
                 c(as.matrix(lm.data[, paste0(c("sin", "cos"), x), with = FALSE]) %*% 
                     matrix(stats::coef(lm_s)[paste0(c("sin", "cos"), x)], ncol = 1)), 
               ncol = 1)
      })))
      colnames(lm.data)[(ncol(lm.data) - length(harmonics) + 1):ncol(lm.data)] = paste0("St", floor(harmonics))
      lm.data[, colnames(lm.data)[grepl("sin|cos", colnames(lm.data))] := NULL]
      bp = suppressWarnings(strucchange::breakpoints(y ~ ., data = lm.data,
                                                     h = ifelse(h > 0.33*length(dat), round(0.33*length(dat)), h), 
                                                     alpha = sig_level, hpc = "foreach"))
    }else if(j == "Ct"){
      dat_trans = dat = series[, c(j), with = FALSE][[1]]
      h = ifelse(!is.na(model$cycle), round(model$cycle) + 2, round(freq) + 2)
      arima = tryCatch(forecast::Arima(y = dat_trans, order = c(2, 0, 1)), 
                       error = function(err){NULL})
      lm.data = data.table(y = dat_trans)[, "lag1" := shift(y, type = "lag", n = 1)][, "lag2" := shift(y, type = "lag", n = 2)]
      if(!is.null(arima)){
        lm.data[, "error" := arima$residuals][, "error" := shift(error, type = "lag", n = 1)]
      }else{
        lm = stats::lm(y ~ lag1 + lag2, data = lm.data)
        lm.data[(.N - length(lm$residuals) + 1):.N, "error" := lm$residuals][, "error" := shift(error, type = "lag", n = 1)]
      }
      bp = suppressWarnings(strucchange::breakpoints(y ~ lag1 + lag2 + error, data = lm.data,
                                                     h = ifelse(h > 0.33*length(dat), round(0.33*length(dat)), h), 
                                                     alpha = sig_level, hpc = "foreach"))
    }
    
    if(!all(is.na(bp$breakpoints))){
      if(j == "Tt0"){
        lm_dat = data.table(y = dat[(length(dat) - length(dat_trans) + 1):length(dat)], 
                            segment = strucchange::breakfactor(bp, breaks = length(bp$breakpoints)))[, "t" := 1:.N, by = "segment"]
        lm = stats::lm(y ~ t*segment, lm_dat)
      }else if(j == "St"){
        lm_dat = data.table(breakpoints = strucchange::breakfactor(bp, breaks = length(bp$breakpoints)))
        lm_dat[, "y" := dat[abs(length(dat) - length(lm_dat$breakpoints) + 1):length(dat)]]
        lm = stats::lm(y ~ breakpoints, lm_dat)
      }else if(j == "Ct"){
        lm_dat = data.table(breakpoints = strucchange::breakfactor(bp, breaks = length(bp$breakpoints)))
        lm_dat[, "y" := dat_trans[abs(length(dat_trans) - length(lm_dat$breakpoints) + 1):length(dat_trans)]]
        lm = stats::lm(y ~ breakpoints, lm_dat)
      }
      fit = data.table(stats::predict(lm, level = ci, interval = "confidence"))
      fit = rbind(matrix(NA, nrow = nrow(series) - nrow(fit), ncol = ncol(fit)), fit, use.names = FALSE)
      colnames(fit) = c(paste0(j, "_interval_mean"), paste0(j, "_interval_lower"), paste0(j, "_interval_upper"))
      series = cbind(series, fit)
      bp = tryCatch(suppressWarnings(data.table(stats::confint(bp, level = ci)$confint)), 
                    error = function(err){temp = data.table(x1 = NA, 
                                                     breakpoints = bp$breakpoints, 
                                                     x2 = NA)
                    colnames(temp)[c(1, 3)] = paste(c(1/2*(1 - ci), 1/2*(1 + ci))*100, "%")
                    return(temp)})
      colnames(bp) = gsub(" ", "", colnames(bp))
      bp[eval(parse(text = paste0("`", colnames(bp)[1], "`"))) < 0, colnames(bp)[1] := 0]
      bp[eval(parse(text = paste0("`", colnames(bp)[3], "`"))) > nrow(series), colnames(bp)[3] := nrow(series)]
      bp = bp[, lapply(.SD, function(x){dates[x]})]
      bp[, "break_type" := ifelse(j == "Tt0", "trend", 
                                  ifelse(j == "St", "seasonal", 
                                         ifelse(j == "Ct", "cycle")))]
      bp_dt = rbind(bp_dt, bp, use.names = TRUE, fill = TRUE)
    }
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
    colnames(series)[grepl("St", colnames(series)) & !grepl("Sts|_interval_", colnames(series))] = 
      gsub("St", "seasonal", colnames(series)[grepl("St", colnames(series)) & !grepl("Sts|_interval_", colnames(series))])
    series[, "seasonal" := rowSums(series[, grepl("seasonal", colnames(series)), with = FALSE])]
  }
  
  #Combine the filtered series
  final = data.table(date = dates, observed = y, 
                     series[, colnames(series) %in% c("trend", "cycle") | 
                              grepl("seasonal|\\%", colnames(series)) |
                              grepl("_interval_", colnames(series)), with = FALSE]
  )
  rm(series, B_tt)
  
  if(nrow(bp_dt) > 0){
    final = merge(final, bp_dt, by.x = "date", by.y = "breakpoints", all = TRUE)
    final[, "date.lag" := shift(date, type = "lag", n = 1)]
    final[, "date.lead" := shift(date, type = "lead", n = 1)]
    suppressWarnings(final[is.na(eval(parse(text = paste0("`", (1/2 - 1/2*ci)*100, "%`")))) & !is.na(break_type), 
                           paste0((1/2 - 1/2*ci)*100, "%") := date.lag])
    suppressWarnings(final[is.na(eval(parse(text = paste0("`", (1/2 + 1/2*ci)*100, "%`")))) & !is.na(break_type), 
                            paste0((1/2 + 1/2*ci)*100, "%") := date.lead])
    suppressWarnings(final[!is.na(break_type) & eval(parse(text = paste0("`", (1/2 - 1/2*ci)*100, "%`"))) < min(date), 
          paste0("`", (1/2 - 1/2*ci)*100, "%`") := min(date)])
    suppressWarnings(final[!is.na(break_type) & eval(parse(text = paste0("`", (1/2 - 1/2*ci)*100, "%`"))) > max(date), 
          paste0("`", (1/2 - 1/2*ci)*100, "%`") := max(date)])
    suppressWarnings(final[!is.na(break_type) & eval(parse(text = paste0("`", (1/2 + 1/2*ci)*100, "%`"))) < min(date), 
          paste0("`", (1/2 - 1/2*ci)*100, "%`") := min(date)])
    suppressWarnings(final[!is.na(break_type) & eval(parse(text = paste0("`", (1/2 + 1/2*ci)*100, "%`"))) > max(date), 
          paste0("`", (1/2 - 1/2*ci)*100, "%`") := max(date)])
    final[, c("date.lag", "date.lead") := NULL]
  }else{
    final[, "break_type" := NA]
  }
  
  parallel::stopCluster(cl)
  
  cols = colnames(final)[colnames(final) %in% c("observed", "trend") | grepl("Tt0_interval_", colnames(final))]
  final[, c(cols) := lapply(.SD, function(x){x*sd + mean}), .SDcols = c(cols)]
  
  cols = colnames(final)[!colnames(final) %in% c("date", "break_type", cols) & !grepl("\\%", colnames(final))]
  final[,  c(cols) := lapply(.SD, function(x){x*sd}), .SDcols = c(cols)]
  
  if(model$multiplicative == TRUE){
    final[, c(colnames(final)[!colnames(final) %in% c("date", "break_type") & !grepl("\\%", colnames(final))]) := lapply(.SD, exp),
          .SDcols = c(colnames(final)[!colnames(final) %in% c("date", "break_type") & !grepl("\\%", colnames(final))])]
  }
  
  if(plot == TRUE){
    if(nrow(final[break_type == "trend", ]) > 0){
      bp = final[break_type == "trend", c(paste0((1/2 - ci/2)*100, "%"), "date", paste0((1/2 + ci/2)*100, "%")), with = FALSE]
      colnames(bp) = c("lower", "date", "upper")
      suppressWarnings(plot(ggplot(data.table::melt(final, id.vars = "date", measure.vars = c("observed", "trend"))) +
                              labs(title = "Trend Break Detection", x = "", y = "") +
                              geom_line(aes(x = date, y = value, group = variable, color = variable)) + theme_minimal() +
                              scale_color_viridis_d() +
                              guides(color = guide_legend(title = NULL), size = guide_legend(title = NULL)) +
                              geom_vline(data = bp, aes(xintercept = date), color = "black") +
                              geom_rect(data = bp,
                                        aes(xmin = lower, xmax = upper, ymin = -Inf, ymax = Inf, fill = "break"), alpha = 0.5) +
                              geom_ribbon(data = final[, c("date", "Tt0_interval_lower", "Tt0_interval_upper"), with = FALSE],
                                          aes(x = date, ymin = Tt0_interval_lower, ymax = Tt0_interval_upper), alpha = 0.5) +
                              geom_line(data = data.table::melt(final[, c("date", "Tt0_interval_mean"), with = FALSE], id.vars = "date"), 
                                        aes(x = date, y = value), color = "red") + 
                              theme(legend.position = "bottom") + 
                              guides(fill = guide_legend(title = NULL))))
    }
    
    if(nrow(final[break_type == "cycle", ]) > 0){
      bp = final[break_type == "cycle", c(paste0((1/2 - ci/2)*100, "%"), "date", paste0((1/2 + ci/2)*100, "%")), with = FALSE]
      colnames(bp) = c("lower", "date", "upper")
      suppressWarnings(plot(ggplot(data.table::melt(final, id.vars = "date", measure.vars = c("cycle"))) +
                              labs(title = "Cycle Break Detection", x = "", y = "") +
                              geom_line(aes(x = date, y = value - mean(final$cycle, na.rm = TRUE), group = variable, color = variable)) +
                              scale_color_viridis_d() +
                              theme_minimal() + guides(color = guide_legend(title = NULL)) +
                              geom_vline(data = bp, aes(xintercept = date), color = "black") +
                              geom_rect(data = bp,
                                        aes(xmin = lower, xmax = upper, ymin = -Inf, ymax = Inf, fill = "break"), alpha = 0.5) +
                              geom_ribbon(data = final[, c("date", "Ct_interval_lower", "Ct_interval_upper"), with = FALSE],
                                          aes(x = date,
                                              ymin = Ct_interval_lower - mean(final$cycle, na.rm = TRUE),
                                              ymax = Ct_interval_upper - mean(final$cycle, na.rm = TRUE)), alpha = 0.5) +
                              geom_line(data = data.table::melt(final[, c("date", "Ct_interval_mean"), with = FALSE], id.vars = "date"), 
                                        aes(x = date, y = value - mean(final$cycle, na.rm = TRUE)), color = "red") + 
                              theme(legend.position = "bottom") + 
                              guides(fill = guide_legend(title = NULL))))
    }
    
    if(any(grepl("seasonal", unique(final$break_type)))){
      bp = final[break_type == "seasonal", c(paste0((1/2 - ci/2)*100, "%"), "date", paste0((1/2 + ci/2)*100, "%")), with = FALSE]
      colnames(bp) = c("lower", "date", "upper")
      suppressWarnings(plot(ggplot(data.table::melt(final, id.vars = "date", measure.vars = "seasonal")) +
                              labs(title = "Seasonal Break Detection", x = "", y = "") +
                              geom_line(aes(x = date, y = value - mean(final$seasonal, na.rm = TRUE), group = variable, color = variable)) +
                              scale_color_viridis_d() +
                              theme_minimal() + guides(color = guide_legend(title = NULL)) +
                              geom_vline(data = bp, aes(xintercept = date), color = "black") +
                              geom_rect(data = bp,
                                        aes(xmin = lower, xmax = upper, ymin = -Inf, ymax = Inf, fill = "break"), alpha = 0.5) +
                              geom_ribbon(data = final[, c("date", "St_interval_lower", "St_interval_upper"), with = FALSE],
                                          aes(x = date,
                                              ymin = eval(parse(text = "St_interval_lower")) - mean(final$seasonal, na.rm = TRUE),
                                              ymax = eval(parse(text = "St_interval_upper")) - mean(final$seasonal, na.rm = TRUE)), alpha = 0.5) +
                              geom_line(data = data.table::melt(final[, c("date", "St_interval_mean"), with = FALSE], id.vars = "date"), 
                                        aes(x = date, y = value - mean(final$seasonal, na.rm = TRUE)), color = "red") + 
                              theme(legend.position = "bottom") + 
                              guides(fill = guide_legend(title = NULL))))
    }
  }
  final = final[, colnames(final) %in% c("date", "break_type", paste0((1/2 - ci/2)*100, "%"), paste0((1/2 + ci/2)*100, "%")), with = FALSE]
  if(any(grepl("\\%", colnames(final)))){
    colnames(final)[grepl("\\%", colnames(final))] = paste0(colnames(final)[grepl("\\%", colnames(final))], "_break")
  }
  return(final)
}
