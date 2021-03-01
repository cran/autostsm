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
#' @param show_progress Whether to show progress bar
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
                              sig_level = 0.01, ci = 0.8, smooth = TRUE, plot = FALSE, cores = NULL, show_progress = FALSE){
  #model = stsm
  #sig_level = 0.01
  #ci = 0.8
  #plot = smooth = TRUE
  #freq = exo = cores = NULL
  #components = c("trend", "cycle", "seasonal")
  #stsm_build_dates = autostsm:::stsm_build_dates
  #stsm_init_vals = autostsm:::stsm_init_vals
  #stsm_check_y = autostsm:::stsm_check_y
  #stsm_check_exo = autostsm:::stsm_check_exo
  #stsm_format_exo = autostsm:::stsm_format_exo
  #Rcpp::sourceCpp("src/kalmanfilter.cpp")
  
  #Bind data.table variables to the global environment
  break_type = date.lag = date.lead = variable = j = x = value = 
    breakpoints = Estimate = `Std. Error` = lower = upper = 
    Tt0_interval_lower = Tt0_interval_upper = . = test_val = 
    V2 = lwr = upr = group = lwr_sig = upr_sig = break_type2 = 
    interval_end = start = end = segment = Period = season = text = 
    Ct_amp_interval_lower = Ct_amp_interval_upper = NULL
  
  #Argument checks
  if(!is.logical(plot)){
    stop("plot must be TRUE or FALSE")
  }
  if(sig_level < 0 | sig_level > 0.1){
    stop("sig_level must be > 0 and <= 0.1")
  }
  if(ci <= 0 | ci >= 1){
    stop("ci must between 0 and 1. Suggested values are 0.8, 0.9, 0.95, 0.99.")
  }
  if(any(!components %in% c("trend", "cycle", "seasonal"))){
    stop("components must be 'trend', 'cycle', or 'seasonal'")
  }
  if(!all(sapply(c(smooth, plot, show_progress), is.logical))){
    stop("smooth, plot, show_progress must be logical")
  }
  if(!all(components %in% c("trend", "cycle", "seasonal"))){
    stop("components must be trend, cycle, seasonal.")
  }
  if(!is.null(freq)){
    if(!is.numeric(freq)){
      stop("freq must be numeric")
    }
  }
  stsm_check_y(y)
  stsm_check_exo(exo, y)
  
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
  exo = stsm_format_exo(exo, dates, range)
  
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
  if(!is.na(model$seasons)){
    seasons = as.numeric(strsplit(model$seasons, ", ")[[1]])
  }else{
    seasons = c()
  }
  if(!is.na(model$cycle)){
    cycle = model$cycle
  }else{
    cycle = c()
  }
  
  #Get the coefficients
  par = eval(parse(text = paste0("c(", model$coef, ")")))
  init = stsm_init_vals(y, par, freq, trend, decomp, seasons, NULL, cycle)
  
  #Filter and smooth the data
  sp = stsm_ssm(par, y, decomp, trend, init)
  B_tt = kalman_filter(matrix(sp$B0, ncol = 1), sp$P0, sp$Dt, sp$At, sp$Ft, sp$Ht, sp$Qt, sp$Rt,
                       matrix(y, nrow = 1), X, sp$beta)[c("B_tl", "B_tt", "P_tt", "P_tl")]
  if(smooth == TRUE){
    B_tt = kalman_smoother(B_tt$B_tl, B_tt$B_tt, B_tt$P_tl, B_tt$P_tt, sp$Ft)$B_tt
  }else{
    B_tt = B_tt$B_tt
  }
  rownames(B_tt) = rownames(sp$Ft)
  
  #Setup parallel computing
  cl = parallel::makeCluster(max(c(1, ifelse(is.null(cores), parallel::detectCores() - 1, 2))))
  doSNOW::registerDoSNOW(cl)
  `%fun%` = foreach::`%dopar%`
  
  #Get the unobserved series
  series = data.table(t(B_tt))
  
  #Detect structural breaks
  comps = c("Tt0", "Ct", "St")[c(TRUE, grepl("cycle", decomp), grepl("seasonal", decomp))]
  
  #Detect structural breaks
  combine = function(x, y){
    return(list(bp_dt = rbind(x[["bp_dt"]], y[["bp_dt"]], use.names = TRUE, fill = TRUE), 
                series_fit = cbind(x[["series_fit"]], y[["series_fit"]])))
  }
  iter = comps[comps %in% c("Tt0", "Ct", "St")[c("trend", "cycle", "seasonal") %in% components]]
  if(show_progress == TRUE){
      pb = progress::progress_bar$new(
        format = " [:bar] :percent complete |:elapsed elapsed |:eta remaining ",
        total = length(iter),  
        clear = FALSE, show_after = 0
    )
    invisible(pb$tick(0))
    progress = function(n){pb$tick()}
  }else{
    progress = NULL
  }
  out = foreach::foreach(j = iter, .combine = combine, 
                         .packages = c("data.table", "forecast", "strucchange", "stats", "tsutils", "foreach", "doSNOW", "parallel"),
                         .options.snow = list(progress = progress)) %fun% {
    bp_dt = data.table()
    series_fit = data.table()
    if(j == "Tt0"){
      #Trend break test
      dat = series[, c(j), with = FALSE][[1]]
      #Test if trend is trending up/down to test for breaks in growth, if not then breaks in level
      ol = forecast::tsoutliers(stats::ts(dat, frequency = freq))
      dat2 = copy(dat)
      dat2[ol$index] = ol$replacements
      if(tsutils::coxstuart(dat2, type = "trend")$p.value <= sig_level){
        dat_trans = diff(dat)
      }else{
        dat_trans = dat
      }
      h = round(freq) + 2
      bp = suppressWarnings(strucchange::breakpoints(dat_trans ~ 1, hpc = "foreach"))
      
      if(!all(is.na(bp$breakpoints))){
        #Trend break test results
        lm_dat = data.table(y = dat[(length(dat) - length(dat_trans) + 1):length(dat)], 
                            segment = strucchange::breakfactor(bp, breaks = length(bp$breakpoints)))[, "t" := 1:.N]
        lm = stats::lm(y ~ t*segment, lm_dat)
        lm2 = stats::lm(y ~ t, lm_dat)
        fit = data.table(stats::predict(lm, level = 1 - sig_level, interval = "confidence"), 
                         stats::predict(lm2, interval = "none"))
        if(nrow(fit[!(V2 >= lwr & V2 <= upr), ]) > 0){
          fit = data.table(stats::predict(lm, level = ci, interval = "confidence"))
          fit = rbind(matrix(NA, nrow = nrow(series) - nrow(fit), ncol = ncol(fit)), fit, use.names = FALSE)
          colnames(fit) = c(paste0(j, "_interval_mean"), paste0(j, "_interval_lower"), paste0(j, "_interval_upper"))
          series_fit = cbind(series_fit, fit)
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
          bp[, "break_type" := "trend"]
          bp_dt = rbind(bp_dt, bp, use.names = TRUE, fill = TRUE)
        }
      }
    }else if(j %in% c("Ct", "St")){
      #Seasonal/cycle amplitude break test
      dat = series[, grepl(j, colnames(series)) & !grepl(paste0(j, "s"), colnames(series)), with = FALSE]
      iter = colnames(dat)
      if(j == "St"){
        h = round(max(seasons)) + 2
      }else if(j == "Ct"){
        h = ifelse(!is.na(model$cycle), round(model$cycle) + 2, round(freq) + 2)
      }
      
      for(k in iter){
        dat_trans = abs(dat[, k, with = FALSE][[1]])
        bp1 = suppressWarnings(strucchange::breakpoints(dat_trans ~ 1,
                                                       h = min(c(h, round(0.33*length(dat_trans)))), 
                                                       hpc = "foreach"))
        
        if(!all(is.na(bp1$breakpoints))){
          #Seasonal/cycle break test results
          lm_dat = data.table(breakpoints = strucchange::breakfactor(bp1, breaks = length(bp1$breakpoints)))
          lm_dat[, "y" := dat_trans[abs(length(dat_trans) - length(lm_dat$breakpoints) + 1):length(dat_trans)]]
          lm = stats::lm(y ~ breakpoints, lm_dat)
          fit = data.table(stats::predict(lm, level = 1 - sig_level, interval = "confidence"))
          if(nrow(fit[!(mean(abs(dat_trans), na.rm = TRUE) >= lwr & mean(abs(dat_trans), na.rm = TRUE) <= upr), ]) > 0){
            fit = pi/2*data.table(stats::predict(lm, level = ci, interval = "confidence"))
            fit = rbind(matrix(NA, nrow = nrow(series) - nrow(fit), ncol = ncol(fit)), fit, use.names = FALSE)
            colnames(fit) = c(paste0(k, "_amp_interval_mean"), paste0(k, "_amp_interval_lower"), paste0(k, "_amp_interval_upper"))
            series_fit = cbind(series_fit, fit)
            bp = tryCatch(suppressWarnings(data.table(stats::confint(bp1, level = ci)$confint)), 
                          error = function(err){temp = data.table(x1 = NA, 
                                                                  breakpoints = bp1$breakpoints, 
                                                                  x2 = NA)
                          colnames(temp)[c(1, 3)] = paste(c(1/2*(1 - ci), 1/2*(1 + ci))*100, "%")
                          return(temp)})
            colnames(bp) = gsub(" ", "", colnames(bp))
            bp[eval(parse(text = paste0("`", colnames(bp)[1], "`"))) < 0, colnames(bp)[1] := 0]
            bp[eval(parse(text = paste0("`", colnames(bp)[3], "`"))) > nrow(series), colnames(bp)[3] := nrow(series)]
            bp = bp[, lapply(.SD, function(x){dates[x]})]
            bp[, "break_type" := ifelse(j == "St", "seasonal", 
                                        ifelse(j == "Ct", "cycle"))]
            bp[, "break_type2" := "amp"]
            bp[, "season" := ifelse(j == "St", gsub("St", "", k), "")]
            bp_dt = rbind(bp_dt, bp, use.names = TRUE, fill = TRUE)
          }
        }
      }
      
      if(j == "Ct"){ #Only use for cycle not seasonality
        #Seasonal/Cycle periodicity break test
        dat = series[, grepl(j, colnames(series)) & !grepl("_", colnames(series)), with = FALSE]
        dat[, "y" := rowSums(dat[, grepl(j, colnames(dat)) & !grepl(paste0(j, "s"), colnames(dat)), with = F])]
        dat[, "ys" := rowSums(dat[, grepl(paste0(j, "s"), colnames(dat)), with = F])]
        dat_trans = dat$y
        dat[, paste0(colnames(dat), ".l1") := lapply(.SD, function(x){shift(x, type = "lag", n = 1)})]
        formula = formula(paste("y ~", paste(colnames(dat)[!grepl("y", colnames(dat)) & grepl("\\.l1", colnames(dat))], collapse = " + "), " - 1"))
        bp2 = suppressWarnings(strucchange::breakpoints(formula, data = dat,
                                                       h = min(c(h, round(0.33*length(dat_trans)))),
                                                       hpc = "foreach"))
        
        if(!all(is.na(bp2$breakpoints))){
          #Seasonal/cycle periodicity break test results
          lm_dat = data.table(breakpoints = strucchange::breakfactor(bp2, breaks = length(bp2$breakpoints)))
          lm_dat = cbind(lm_dat, dat[abs(length(dat_trans) - length(lm_dat$breakpoints) + 1):length(dat_trans), ])
          fit = do.call("rbind", lapply(unique(lm_dat$breakpoints), function(x){
            lm = stats::lm(formula, data = lm_dat[breakpoints == x, ])
            coef = data.table(summary(lm)$coefficients, keep.rownames = T)
            coef[, "upr_sig" := Estimate + stats::qnorm(1 - sig_level/2)*`Std. Error`]
            coef[, "lwr_sig" := Estimate + stats::qnorm(sig_level/2)*`Std. Error`]
            coef[, "upr" := Estimate + stats::qnorm(1/2 + ci/2)*`Std. Error`]
            coef[, "lwr" := Estimate + stats::qnorm(1/2 - ci/2)*`Std. Error`]
            lm2 = stats::lm(formula, data = lm_dat)
            coef2 = data.table(summary(lm2)$coefficients)
            #rho = beta1/cos(lambda)
            #lambda = atan(beta2/beta1)
            fit = abs(sapply(seq(1, nrow(coef), 2), function(y){
              2*pi/atan(coef[y + 1, ]$Estimate/coef[y, ]$Estimate)
            }))
            fit2 = abs(sapply(seq(1, nrow(coef2), 2), function(y){
              2*pi/atan(coef2[y + 1, ]$Estimate/coef2[y, ]$Estimate)
            }))
            confint = abs(unlist(lapply(c("upr", "lwr"), function(y){
              unlist(lapply(c("upr", "lwr"), function(z){
                sapply(seq(1, nrow(coef), 2), function(w){
                  2*pi/atan(coef[w + 1, c(y), with = FALSE][[1]]/coef[w, c(z), with = FALSE][[1]])
                })
              }))
            })))
            sigint = abs(unlist(lapply(c("upr_sig", "lwr_sig"), function(y){
              unlist(lapply(c("upr_sig", "lwr_sig"), function(z){
                sapply(seq(1, nrow(coef), 2), function(w){
                  2*pi/atan(coef[w + 1, c(y), with = FALSE][[1]]/coef[w, c(z), with = FALSE][[1]])
                })
              }))
            })))
            if(j == "St"){
              names(fit) = seasons
              names(confint) = rep(seasons, length(confint)/2)
              names(sigint) = rep(seasons, length(sigint)/2)
              group = names(fit)
            }else if(j == "Ct"){
              names(fit) = ""
              names(confint) = rep("", length(confint))
              names(sigint) = rep("", length(sigint))
              group = model$cycle
            }
            return(data.table(breakpoints = x, group = group, fit = round(fit), 
                              lwr = floor(sapply(unique(names(confint)), function(y){min(confint[names(confint) == y])})), 
                              upr = ceiling(sapply(unique(names(confint)), function(y){max(confint[names(confint) == y])})), 
                              lwr_sig = floor(sapply(unique(names(sigint)), function(y){min(sigint[names(sigint) == y])})), 
                              upr_sig = ceiling(sapply(unique(names(sigint)), function(y){max(sigint[names(sigint) == y])})), 
                              test_val = round(fit2)))
          }))
          if(nrow(fit[!(test_val >= lwr_sig & test_val <= upr_sig), ]) > 0){
            fit[, c("lwr_sig", "upr_sig", "test_val") := NULL]
            fit = do.call("cbind", lapply(unique(fit$group), function(x){
              temp = data.table::merge.data.table(fit[group == x, ], lm_dat[, "breakpoints"], all.x = TRUE)[, "breakpoints" := NULL]
              temp = rbind(matrix(NA, nrow = nrow(series) - nrow(temp), ncol = ncol(temp)), temp, use.names = FALSE)
              colnames(temp) = c("group", paste0(j, ifelse(j == "St", x, ""), "_period_interval_mean"), 
                                 paste0(j, ifelse(j == "St", x, ""), "_period_interval_lower"), 
                                 paste0(j, ifelse(j == "St", x, ""), "_period_interval_upper"))
              return(temp[, 2:ncol(temp)])
            }))
            series_fit = cbind(series_fit, fit)
            bp = tryCatch(suppressWarnings(data.table(stats::confint(bp2, level = ci)$confint)), 
                          error = function(err){temp = data.table(x1 = NA, 
                                                                  breakpoints = bp2$breakpoints, 
                                                                  x2 = NA)
                          colnames(temp)[c(1, 3)] = paste(c(1/2*(1 - ci), 1/2*(1 + ci))*100, "%")
                          return(temp)})
            colnames(bp) = gsub(" ", "", colnames(bp))
            bp[eval(parse(text = paste0("`", colnames(bp)[1], "`"))) < 0, colnames(bp)[1] := 0]
            bp[eval(parse(text = paste0("`", colnames(bp)[3], "`"))) > nrow(series), colnames(bp)[3] := nrow(series)]
            bp = bp[, lapply(.SD, function(x){dates[x]})]
            bp[, "break_type" := ifelse(j == "St", "seasonal", 
                                        ifelse(j == "Ct", "cycle"))]
            bp[, "break_type2" := "period"]
            bp_dt = rbind(bp_dt, bp, use.names = TRUE, fill = TRUE)
          }
        }
      }
    }
    suppressWarnings(rm(bp, fit, bp, bp1, bp2, dat, dat_trans, lm, lm_dat, formula, h, k, lm2, ol, dat2))
    gc()
    return(list(bp_dt = bp_dt, series_fit = series_fit))
  }
  bp_dt = out[["bp_dt"]]
  series = cbind(series, out[["series_fit"]])
  rm(out)
  
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
          paste0((1/2 - 1/2*ci)*100, "%") := min(date)])
    suppressWarnings(final[!is.na(break_type) & eval(parse(text = paste0("`", (1/2 - 1/2*ci)*100, "%`"))) > max(date), 
          paste0((1/2 - 1/2*ci)*100, "%") := max(date)])
    suppressWarnings(final[!is.na(break_type) & eval(parse(text = paste0("`", (1/2 + 1/2*ci)*100, "%`"))) < min(date), 
          paste0((1/2 - 1/2*ci)*100, "%") := min(date)])
    suppressWarnings(final[!is.na(break_type) & eval(parse(text = paste0("`", (1/2 + 1/2*ci)*100, "%`"))) > max(date), 
          paste0((1/2 - 1/2*ci)*100, "%") := max(date)])
    final[, c("date.lag", "date.lead") := NULL]
  }else{
    final[, "break_type" := NA]
  }
  
  parallel::stopCluster(cl)
  
  cols = colnames(final)[colnames(final) %in% c("observed", "trend") | grepl("Tt0_interval_", colnames(final))]
  final[, c(cols) := lapply(.SD, function(x){x*sd + mean}), .SDcols = c(cols)]
  
  cols = colnames(final)[!colnames(final) %in% c("date", "break_type", "break_type2", "season", cols) & !grepl("\\%", colnames(final)) & !grepl("period", colnames(final))]
  final[,  c(cols) := lapply(.SD, function(x){x*sd}), .SDcols = c(cols)]
  
  if(model$multiplicative == TRUE){
    final[, c(colnames(final)[!colnames(final) %in% c("date", "break_type", "break_type2", "season") & !grepl("\\%", colnames(final)) & !grepl("period", colnames(final))]) := lapply(.SD, exp),
          .SDcols = c(colnames(final)[!colnames(final) %in% c("date", "break_type", "break_type2", "season") & !grepl("\\%", colnames(final)) & !grepl("period", colnames(final))])]
  }
  
  if(plot == TRUE){
    if(nrow(final[break_type == "trend", ]) > 0){
      bp = final[break_type == "trend", c(paste0((1/2 - ci/2)*100, "%"), "date", paste0((1/2 + ci/2)*100, "%")), with = FALSE]
      colnames(bp) = c("lower", "date", "upper")
      suppressWarnings(plot(ggplot(data.table::melt(final, id.vars = "date", measure.vars = c("observed", "trend"))) +
                              labs(title = "Trend Break Detection", subtitle = paste0(round(ci*100), "% Confidence Interval"), x = "", y = "") +
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
      bp = final[break_type == "cycle", c(paste0((1/2 - ci/2)*100, "%"), "date", paste0((1/2 + ci/2)*100, "%"), "break_type2"), with = FALSE]
      colnames(bp) = c("lower", "date", "upper", "break_type")
      
      bp2 = final[break_type == "cycle" & break_type2 == "period", c(paste0((1/2 - ci/2)*100, "%"), "date", paste0((1/2 + ci/2)*100, "%")), with = FALSE]
      if(nrow(bp2) > 0){
        bp2 = rbind(data.table(date = min(final$date)), bp2, data.table(date = max(final$date)), use.names = TRUE, fill = TRUE)
        bp2[, "interval_start" := date][, "interval_end" := shift(date, type = "lead", n = 1)]
        bp2 = bp2[!is.na(interval_end), ]
        bp2 = cbind(bp2, do.call("rbind", lapply(1:nrow(bp2), function(x){
          temp = final[date >= bp2[x, ]$interval_start & date <= bp2[x, ]$interval_end, c(colnames(final)[grepl("period", colnames(final))]), with = FALSE]
          colMeans(temp[, grepl("Ct_", colnames(temp)), with = FALSE], na.rm = TRUE)
        })))
        periods = tryCatch(bp2[, 4:ncol(bp2)], error = function(err){NULL})
        if(!is.null(periods)){
          periods[, colnames(periods)[3:ncol(periods)] := lapply(.SD, round), 
                .SDcols = colnames(periods)[3:ncol(periods)]]
          colnames(periods) = gsub("Ct|interval|period|_", "", colnames(periods))
          colnames(periods)[grepl("mean", colnames(periods))] = gsub("mean", " Period", colnames(periods)[grepl("mean", colnames(periods))])
          colnames(periods)[grepl("lower", colnames(periods))] = gsub("lower", paste0(" ", (1/2 - ci/2)*100, "%"), colnames(periods)[grepl("lower", colnames(periods))])
          colnames(periods)[grepl("upper", colnames(periods))] = gsub("upper", paste0(" ", (1/2 + ci/2)*100, "%"), colnames(periods)[grepl("upper", colnames(periods))])
          periods[, "x" := mean(c(start, end)), by = "start"]
          periods[, "y" := mean(final$cycle, na.rm = T)]
          periods[, "segment" := 1:.N]
          colnames(periods) = trimws(colnames(periods))
          periods[, "text" := paste0("Segment ", segment, "\nPeriod: ", Period, " [",
                                     eval(parse(text = paste0("`", (1/2 - ci/2)*100, "%`"))), ", ",
                                     eval(parse(text = paste0("`", (1/2 + ci/2)*100, "%`"))), "]")]
          periods = periods[, .(text = paste(text, collapse = "\n")), by = c("start", "end", "x", "y")]
        }
      }else{
        periods = NULL
      }
        
      adj = mean(final$seasonal, na.rm = TRUE)
      subtitle = paste0("Amplitude", ifelse(!is.null(periods),  " & Periodicity", ""), ", ", round(ci*100), "% Confidence Interval")
      g1 = ggplot(data.table::melt(final, id.vars = "date", measure.vars = c("cycle"))) +
          labs(title = "Cycle Break Detection", subtitle = subtitle, x = "", y = "") +
          geom_line(aes(x = date, y = value, group = variable, color = variable)) +
          scale_color_viridis_d() +
          theme_minimal() + guides(color = guide_legend(title = NULL)) +
          geom_vline(data = bp, aes(xintercept = date), color = "black") +
          geom_hline(aes(yintercept = ifelse(model$multiplicative == TRUE, 1, 0)), color = "black") + 
          geom_rect(data = bp,
                    aes(xmin = lower, xmax = upper, ymin = -Inf, ymax = Inf, fill = break_type), alpha = 0.5) + 
          theme(legend.position = "bottom") + 
          guides(fill = guide_legend(title = NULL))
        
      if(any(grepl("Ct_amp", colnames(final)))){
        g1 = g1 + geom_ribbon(data = final[, c("date", "Ct_amp_interval_lower", "Ct_amp_interval_upper"), with = FALSE],
                      aes(x = date,
                          ymin = Ct_amp_interval_lower,
                          ymax = Ct_amp_interval_upper), alpha = 0.5) +
          geom_ribbon(data = final[, c("date", "Ct_amp_interval_lower", "Ct_amp_interval_upper"), with = FALSE],
                      aes(x = date,
                          ymin = -(Ct_amp_interval_lower - adj) + adj,
                          ymax = -(Ct_amp_interval_upper - adj) + adj), alpha = 0.5) + 
          geom_line(data = data.table::melt(final[, c("date", "Ct_amp_interval_mean"), with = FALSE], id.vars = "date"), 
                    aes(x = date, y = value), color = "red") + 
          geom_line(data = data.table::melt(final[, c("date", "Ct_amp_interval_mean"), with = FALSE], id.vars = "date"), 
                    aes(x = date, y = -(value - adj) + adj), color = "red")
      }

      if(!is.null(periods)){
        g1 = g1 + ggrepel::geom_label_repel(data = periods, aes(x = x, y = y, label = text), alpha = 0.8)
      }
      suppressWarnings(plot(g1))
    }

    if(any(grepl("seasonal", unique(final$break_type)))){
      g = list()
      for(k in unique(bp_dt[break_type == "seasonal" & !is.na(season), ]$season)){
        bp = final[break_type == "seasonal" & (season == k | is.na(season)), c(paste0((1/2 - ci/2)*100, "%"), "date", paste0((1/2 + ci/2)*100, "%"), "break_type2"), with = FALSE]
        colnames(bp) = c("lower", "date", "upper", "break_type")
        
        bp2 = final[break_type == "seasonal" & break_type2 == "period", c(paste0((1/2 - ci/2)*100, "%"), "date", paste0((1/2 + ci/2)*100, "%")), with = FALSE]
        if(nrow(bp2) > 0){
          bp2 = rbind(data.table(date = min(final$date)), bp2, data.table(date = max(final$date)), use.names = TRUE, fill = TRUE)
          bp2[, "interval_start" := date][, "interval_end" := shift(date, type = "lead", n = 1)]
          bp2 = bp2[!is.na(interval_end), ]
          bp2 = cbind(bp2, do.call("rbind", lapply(1:nrow(bp2), function(x){
            temp = final[date >= bp2[x, ]$interval_start & date >= bp2[x, ]$interval_end, c(colnames(final)[grepl("period", colnames(final))]), with = FALSE]
            colMeans(temp[, grepl(paste(paste0("St", seasons), collapse = "|"), colnames(temp)), with = FALSE], na.rm = TRUE)
          })))
          periods = tryCatch(bp2[, 4:ncol(bp2)], error = function(err){NULL})
          if(!is.null(periods)){
            periods[, colnames(periods)[3:ncol(periods)] := lapply(.SD, round), 
                  .SDcols = colnames(periods)[3:ncol(periods)]]
            colnames(periods) = gsub("St|interval|period|_", "", colnames(periods))
            colnames(periods)[grepl("mean", colnames(periods))] = gsub("mean", " Period", colnames(periods)[grepl("mean", colnames(periods))])
            colnames(periods)[grepl("lower", colnames(periods))] = gsub("lower", paste0(" ", (1/2 - ci/2)*100, "%"), colnames(periods)[grepl("lower", colnames(periods))])
            colnames(periods)[grepl("upper", colnames(periods))] = gsub("upper", paste0(" ", (1/2 + ci/2)*100, "%"), colnames(periods)[grepl("upper", colnames(periods))])
            periods = data.table::melt(periods, id.vars = c("start", "end"))
            periods$season = sapply(periods$variable, function(x){strsplit(as.character(x), " ")[[1]][1]})
            periods$variable = sapply(periods$variable, function(x){strsplit(as.character(x), " ")[[1]][2]})
            periods = dcast(periods, "start + end + season ~ variable", value.var = "value")
            periods[, "x" := mean(c(start, end)), by = "start"]
            periods[, "y" := mean(final[, paste0("seasonal", k), with = FALSE][[1]], na.rm = T)]
            colnames(periods) = trimws(colnames(periods))
            periods = periods[season == k, ]
            periods[, "text" := paste0(Period, " [",
                                       eval(parse(text = paste0("`", (1/2 - ci/2)*100, "%`"))), ", ",
                                       eval(parse(text = paste0("`", (1/2 + ci/2)*100, "%`"))), "]")]
            periods = periods[, .(text = paste(text, collapse = "\n")), by = c("start", "end", "x", "y")]
            periods[, "segment" := 1:.N]
            periods[, "text" := paste0("Segment ", segment, "\n", text)]
          }
        }else{
          periods = NULL
        }
        
        adj = mean(final[, paste0("seasonal", k), with = FALSE][[1]], na.rm = TRUE)
        toplot = data.table::melt(final, id.vars = "date", measure.vars = c(paste0("seasonal", k)))
        subtitle = paste0("Amplitude", ifelse(!is.null(periods),  " & Periodicity", ""), ", ", round(ci*100), "% Confidence Interval")
        g1 = ggplot(toplot) +
          labs(title = paste("Seasonal", k, "Break Detection"), subtitle = subtitle, x = "", y = "") +
          geom_line(aes(x = date, y = value, group = variable, color = variable)) +
          scale_color_viridis_d() +
          theme_minimal() + guides(color = guide_legend(title = NULL)) +
          geom_vline(data = bp, aes(xintercept = date), color = "black") +
          geom_hline(aes(yintercept = ifelse(model$multiplicative == TRUE, 1, 0)), color = "black") + 
          geom_rect(data = bp,
                    aes(xmin = lower, xmax = upper, ymin = -Inf, ymax = Inf, fill = break_type), alpha = 0.5) + 
          theme(legend.position = "bottom") + 
          guides(fill = guide_legend(title = NULL))
        if(any(grepl(paste0("St", k, "_amp"), colnames(final)))){
          toplot2 = final[, c("date", paste0("St", k, c("_amp_interval_lower", "_amp_interval_upper", "_amp_interval_mean"))), with = FALSE]
          colnames(toplot2)[grepl("lower", colnames(toplot2))] = "lower"
          colnames(toplot2)[grepl("upper", colnames(toplot2))] = "upper"
          colnames(toplot2)[grepl("mean", colnames(toplot2))] = "mean"
          g1 = g1 + geom_ribbon(data = toplot2,
                      aes(x = date, ymin = lower, ymax = upper), alpha = 0.5) +
          geom_ribbon(data = toplot2,
                      aes(x = date, ymin = -(lower - adj) + adj, ymax = -(upper - adj) + adj), alpha = 0.5) +
          geom_line(data = data.table::melt(toplot2[, c("date", "mean"), with = FALSE], id.vars = "date"), 
                    aes(x = date, y = value), color = "red") + 
          geom_line(data = data.table::melt(toplot2[, c("date", "mean"), with = FALSE], id.vars = "date"), 
                    aes(x = date, y = -(value - adj) + adj), color = "red")
        }
     
        if(!is.null(periods)){
          g1 = g1 + ggrepel::geom_label_repel(data = periods, aes(x = x, y = y, label = text), alpha = 0.8)
        }
        g[[k]] = g1
      }
      
      layout_matrix = 1:(2*ceiling(length(g)/2))
      if(length(layout_matrix) > length(g)){
        layout_matrix = c(layout_matrix[1], 1:(length(layout_matrix) - 1))
      }
      layout_matrix = t(matrix(layout_matrix, ncol = 2, byrow = TRUE))
      suppressWarnings(gridExtra::grid.arrange(grobs = g, layout_matrix = layout_matrix))
    }
  }
  final = final[, colnames(final)[grepl("\\%|interval|break|date", colnames(final))], with = FALSE]
  if(any(grepl("\\%", colnames(final)))){
    colnames(final)[grepl("\\%", colnames(final))] = paste0(colnames(final)[grepl("\\%", colnames(final))], "_break")
  }
  colnames(final)[grepl("interval_lower", colnames(final))] = gsub("_lower", paste0(c(1/2 - ci/2)*100, "%"), colnames(final)[grepl("interval_lower", colnames(final))])
  colnames(final)[grepl("interval_upper", colnames(final))] = gsub("_upper", paste0(c(1/2 + ci/2)*100, "%"), colnames(final)[grepl("interval_upper", colnames(final))])
  return(final)
}
