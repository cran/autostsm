#' Detect Anomalies
#'
#' Detect anomalies using the estimated structural time series model
#' @param model Structural time series model estimated using stsm_estimate.
#' @param y Univariate time series of data values. May also be a 2 column data frame containing a date column.
#' @param freq Frequency of the data (1 (yearly), 4 (quarterly), 12 (monthly), 365.25/7 (weekly), 365.25 (daily)), default is NULL and will be automatically detected
#' @param exo_obs Matrix of exogenous variables to be used in the observation equation. 
#' @param exo_state Matrix of exogenous variables to be used in the state matrix. 
#' @param plot Whether to plot everything
#' @param smooth Whether or not to use the Kalman smoother
#' @param sig_level Significance level to determine statistically significant anomalies
#' @return data table (or list of data tables) containing the dates of detected anomalies from the filtered and/or smoothed series
#' @import ggplot2 data.table kalmanfilter
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
#' anomalies = stsm_detect_anomalies(model = stsm, y = NA000334Q, plot = TRUE)
#' }
#' @export
stsm_detect_anomalies = function(model, y = NULL, freq = NULL, exo_obs = NULL, exo_state = NULL,
                                 sig_level = 0.01, smooth = TRUE, plot = FALSE){
  #model = stsm
  #sig_level = 0.01
  #plot = smooth = TRUE
  #freq = exo = NULL
  # for(i in list.files(path = "R", pattern = ".R", full.names = T)){
  #   tryCatch(source(i), error = function(err){NULL})
  # }
  
  #Bind data.table variables to the global environment
  anomaly = predicted = value = variable = observed = NULL
  
  #Argument tests
  if(!is.logical(plot)){
    stop("plot must be TRUE or FALSE")
  }
  if(sig_level < 0 | sig_level > 0.1){
    stop("sig_level must be positive and <= 0.1")
  }
  if(!is.null(freq)){
    if(!is.numeric(freq)){
      stop("freq must be numeric")
    }
  }
  if(!all(sapply(c(smooth, plot), is.logical))){
    stop("smooth, plot must be logical")
  }
  stsm_check_y(y)
  stsm_check_exo(exo_obs, y)
  stsm_check_exo(exo_state, y)
  
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
  exo = stsm_format_exo(exo_obs, exo_state, dates, range)
  
  #Set the historical exogenous variables
  if(is.null(exo_obs)){
    Xo = t(matrix(0, nrow = length(y), ncol = 1))
    Xo[is.na(Xo)] = 0
    rownames(Xo) = "Xo"
  }else{
    Xo = t(exo[, grepl("obs\\.", colnames(exo)), with = FALSE])
  }
  if(is.null(exo_state)){
    Xs = t(matrix(0, nrow = length(y), ncol = 1))
    Xs[is.na(Xs)] = 0
    rownames(Xs) = "Xs"
  }else{
    Xs = t(exo[, grepl("state\\.", colnames(exo)), with = FALSE])
  }
  
  #Apply multiplicative model
  if(model$multiplicative == TRUE){
    y = log(y)
  }
  
  #Standardize
  mean = mean(y, na.rm = TRUE)
  sd = stats::sd(y, na.rm = TRUE)
  y = (y - mean)/sd
  
  #Filter and smooth the data
  ssm = stsm_ssm(yt = y, model = model)
  msg = utils::capture.output(B_tt <- kalman_filter(ssm, matrix(y, nrow = 1), Xo, Xs, smooth = smooth)$B_tt, type = "message")
  rownames(B_tt) = rownames(ssm[["Fm"]])
  
  #Get the unobserved series
  series = data.table(t(B_tt))
  fev = (ssm[["Hm"]] %*% (ssm[["Fm"]] %*% ssm[["Qm"]] %*% t(ssm[["Fm"]])) %*% t(ssm[["Hm"]])) + ssm[["Rm"]]
    
  #Lag the series
  series.l = copy(series)
  series.l[, colnames(series.l) := lapply(.SD, shift, type = "lag", n = 1), .SDcols = colnames(series.l)]
  
  #Get the model errors
  pred_uc = (matrix(ssm[["Dm"]], nrow = nrow(ssm[["Dm"]]), ncol = nrow(series)) + 
               ssm[["Fm"]] %*% t(as.matrix(series.l[, rownames(B_tt), with = FALSE])) + 
               ssm[["betaS"]] %*% Xs)
  pred = t(matrix(ssm[["Am"]], nrow = 1, ncol = ncol(pred_uc)) + ssm[["Hm"]] %*% pred_uc + 
             ssm[["betaO"]] %*% Xo)
  errors = data.table(t(t(as.matrix(series[, rownames(B_tt), with = FALSE])) - pred_uc))
  errors[1:length(y), "pred" := y - pred]
  
  #Test for data anomalies
  anomalies = data.table(date = dates[which(abs(errors$pred) > stats::qnorm(1 - sig_level, 0, sqrt(fev)))], anomaly = TRUE)
  final = merge(data.table(date = dates), anomalies, by = "date", all = TRUE)[is.na(anomaly), "anomaly" := FALSE]
  final = final[!is.na(date), ]
  final[, "predicted" := pred]
  final[, paste0((1 - sig_level)*100, "%") := predicted + stats::qnorm(1 - sig_level, 0, sqrt(fev))]
  final[, paste0((sig_level)*100, "%") := predicted - stats::qnorm(1 - sig_level, 0, sqrt(fev))]
  final[, "observed" := y]
  
  final[, c("observed", "predicted", paste0((1-sig_level)*100, "%"),  paste0((sig_level)*100, "%")) := lapply(.SD, function(x){
    x*sd + mean
  }), .SDcols = c("observed", "predicted", paste0((1-sig_level)*100, "%"),  paste0((sig_level)*100, "%"))]
  
  if(model$multiplicative == TRUE){
    final[, c("observed", "predicted", paste0((1-sig_level)*100, "%"),  paste0((sig_level)*100, "%")) := lapply(.SD, exp), 
          .SDcols = c("observed", "predicted", paste0((1-sig_level)*100, "%"),  paste0((sig_level)*100, "%"))]
  }
  
  if(plot == TRUE){
    if(nrow(final[anomaly == TRUE, ]) > 0){
      suppressWarnings(plot(ggplot(data.table::melt(final, id.vars = c("date"), measure.vars = c("observed"))) +
                              labs(title = "Anomaly Detection", subtitle = paste0((1 - sig_level)*100, "% Threshold Level"), x = "", y = "") +
                              geom_line(aes(x = date, y = value, group = variable, color = variable)) + theme_minimal() +
                              geom_line(data = data.table::melt(final[, c("date", paste0((sig_level)*100, "%")), with = FALSE], id.vars = "date"), 
                                        aes(x = date, y = value, color = "threshold"), alpha = 0.5, linetype = "dotdash") + 
                              geom_line(data = data.table::melt(final[, c("date", paste0((1-sig_level)*100, "%")), with = FALSE], id.vars = "date"), 
                                        aes(x = date, y = value, color = "threshold"), alpha = 0.5, linetype = "dotdash") + 
                              geom_point(data = final[anomaly == TRUE, ], aes(x = date, y = observed, size = "anomaly"), 
                                         color = "red") + 
                              scale_color_viridis_d() +
                              scale_linetype_manual(values = c("solid", "dashed")) + 
                              guides(color = guide_legend(title = NULL), size = guide_legend(title = NULL)) +
                              theme(legend.position = "bottom")))
    }
  }
  final = final[, colnames(final) %in% c("date", "anomaly", paste0((1-sig_level)*100, "%"),  paste0((sig_level)*100, "%")), with = FALSE]
  if(any(grepl("\\%", colnames(final)))){
    colnames(final)[grepl("\\%", colnames(final))] = paste0(colnames(final)[grepl("\\%", colnames(final))], "_anomaly")
  }
  return(final)
}