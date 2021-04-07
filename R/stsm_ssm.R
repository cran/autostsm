#' State space model
#' 
#' Creates a state space model in list form
#' yt = H*B + e_t
#' B = F*B_{t-1} + u_t
#'
#' @param par Vector of named parameter values, includes the harmonics
#' @param yt Univariate time series of data values
#' @param decomp Decomposition model ("tend-cycle-seasonal", "trend-seasonal", "trend-cycle", "trend-noise")
#' @param trend Trend specification ("random-walk", "random-walk-drift", "double-random-walk", "random-walk2"). The default is NULL which will choose the best of all specifications based on the maximum likielhood.
#' "random-walk" is the random walk trend.
#' "random-walk-drift" is the random walk with constant drift trend.
#' "double-random-walk" is the random walk with random walk drift trend.
#' "random-walk2" is a 2nd order random walk trend as in the Hodrick-Prescott filter.
#' @param init Initial state values for the Kalman filter
#' @param model a stsm_estimate model object
#' @import data.table
#' @return List of space space matrices
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
#' ssm = stsm_ssm(model = stsm)
#' }
#' @export
stsm_ssm = function(par = NULL, yt = NULL, decomp = NULL,
                    trend = NULL, init = NULL, model = NULL){
  if(!is.null(model)){
    par = eval(parse(text = paste0("c(", model$coef, ")")))
    trend = model$trend
    decomp = model$decomp
  }
  
  if(!is.null(yt)){
    yt = matrix(yt[!is.na(yt)], nrow = 1)
  }else{
    yt = matrix(rep(0, 2), ncol = 1)
  }
  
  #Define the transition and observation equation matrices based on the trend specification
  if(trend  == "random-walk"){
    #T_t = T_{t-1} + e_t, e_t ~ N(0, sig_t^2)
    #Transition matrix
    Fm = matrix(c(1), nrow = 1, ncol = 1)
    colnames(Fm) = "Tt_1"
    rownames(Fm) = "Tt_0"
    #Observation matrix
    Hm = matrix(c(1), nrow = 1)
    colnames(Hm) = rownames(Fm)
  }else if(trend %in% c("random-walk-drift", "double-random-walk")){
    #T_t = M_{t-1} + T_{t-1} + e_t, e_t ~ N(0, sig_t^2)
    #D_t = D_{t-1} + n_t, n_t ~ N(0, sig_d^2)
    #Transition matrix
    Fm = rbind(c(1, 1),
               c(0, 1))
    colnames(Fm) = c("Tt_1", "Dt_1")
    rownames(Fm) = c("Tt_0", "Dt_0")
    if(trend == "random-walk-drift"){
      Fm["Dt_0", "Dt_1"] = par["phi_d"]
    }
    #Observation matrix
    Hm = matrix(c(1, 0), nrow = 1)
    colnames(Hm) = rownames(Fm)
  }else if(trend %in% c("random-walk2")){
    #T_t = 2T_{t-1} - T_{t-2} + e_t, e_t ~ N(0, sig_t^2)
    #Transition matrix
    Fm = rbind(c(2, -1),
               c(1, 0))
    colnames(Fm) = c("Tt_1", "Tt_2")
    rownames(Fm) = c("Tt_0", "Tt_1")
    #Observation matrix
    Hm = matrix(c(1, 0), nrow = 1)
    colnames(Hm) = rownames(Fm)
  }
  
  #Get the parameters on exogenous data
  if(any(grepl("beta_", names(par)))){
    beta = matrix(par[grepl("beta_", names(par))], nrow = 1, ncol = length(par[grepl("beta_", names(par))]))
    colnames(beta) = gsub("beta_\\.", "", names(par[grepl("beta_", names(par))]))
  }else{
    beta = matrix(0, nrow = 1, ncol = 1)
  }
  
  #Define the transition error covariance matrix
  Qm = matrix(0, nrow = nrow(Fm), ncol = ncol(Fm))
  colnames(Qm) = rownames(Qm) = rownames(Fm)
  Qm[rownames(Qm) == "Tt_0", colnames(Qm) == "Tt_0"] = par["sig_t"]^2
  Qm[rownames(Qm) == "Dt_0", colnames(Qm) == "Dt_0"] = par["sig_d"]^2
  
  #Define the cycle component
  if(grepl("cycle", decomp)){
    Cm = par["phi_c"]*rbind(c(cos(par["lambda"]), sin(par["lambda"])),
                          c(-sin(par["lambda"]), cos(par["lambda"])))
    colnames(Cm) = c("Ct_1", "Cts_1")
    rownames(Cm) = c("Ct_0", "Cts_0")
    Fm = matrix(Matrix::bdiag(Fm, Cm), ncol = ncol(Fm) + ncol(Cm),
                dimnames = list(c(rownames(Fm), rownames(Cm)), c(colnames(Fm), colnames(Cm))))
    Hm = cbind(Hm, matrix(c(1, 0), nrow = 1))
    colnames(Hm) = rownames(Fm)
    Qm = as.matrix(Matrix::bdiag(Qm, diag(2)*par["sig_c"]^2))
    colnames(Qm) = rownames(Qm) = rownames(Fm)
  }else if(sum(grepl("phi_c\\.\\d+|theta_c\\.\\d+", names(par))) > 0){
    #Define the ARMA part of the cycle/noise component
    par["phi_c.2"] = -par["phi_c.1"]^2/4
    len = sum(grepl("phi_c\\.\\d+|theta_c\\.\\d+", names(par))) 
    Cm = rbind(par[grepl("phi_c|theta_c", names(par))], 
               matrix(0, nrow = len - 1, ncol = len))
    colnames_c = rownames_c = c()
    if(sum(grepl("phi_c", names(par))) > 0){
      colnames_c = c(colnames_c, paste0("Ct_", gsub("phi_c\\.", "", names(par)[grepl("phi_c\\.\\d+", names(par))])))
      rownames_c = c(rownames_c, paste0("Ct_", as.numeric(gsub("phi_c\\.", "", names(par)[grepl("phi_c\\.\\d+", names(par))])) - 1))
    }
    if(sum(grepl("theta_c", names(par))) > 0){
      colnames_c = c(colnames_c, paste0("et_", gsub("theta_c\\.", "", names(par)[grepl("theta_c\\.\\d+", names(par))])))
      rownames_c = c(rownames_c, paste0("et_", as.numeric(gsub("theta_c\\.", "", names(par)[grepl("theta_c\\.\\d+", names(par))])) - 1))
    }
    colnames(Cm) = colnames_c
    rownames(Cm) = rownames_c
    
    Cm[rownames(Cm) %in% colnames(Cm), colnames(Cm) %in% rownames(Cm)] = diag(sum(rownames(Cm) %in% colnames(Cm)))
    Fm = matrix(Matrix::bdiag(Fm, Cm), ncol = ncol(Fm) + ncol(Cm),
                dimnames = list(c(rownames(Fm), rownames(Cm)), c(colnames(Fm), colnames(Cm))))
    Hm = tryCatch(cbind(Hm, matrix(as.numeric(rownames(Cm) == "Ct_0"), nrow = 1)), 
                  error = function(err){matrix(as.numeric(rownames(Cm) == "Ct_0"), nrow = 1)})
    colnames(Hm) = rownames(Fm)
    Qm = as.matrix(Matrix::bdiag(Qm, diag(as.numeric(rownames(Cm) %in% c("Ct_0", "et_0"))*par["sig_c"]^2)))
    colnames(Qm) = rownames(Qm) = rownames(Fm)
    
    # Cm = rbind(par[c("phi_c.1", "phi_c.2", "theta_c.1")])
    # Cm = matrix(rbind(Cm,
    #            cbind(diag(2)*c(1, 0), matrix(0, nrow = 2, ncol = 1))),
    #            ncol = ncol(Cm), dimnames = list(NULL, colnames(Cm)))
    # rownames(Cm) = c("Ct_0", "Ct_1", "et_0")
    # colnames(Cm) = c("Ct_1", "Ct_2", "et_1")
    # Fm = matrix(Matrix::bdiag(Fm, Cm), ncol = ncol(Fm) + ncol(Cm),
    #             dimnames = list(c(rownames(Fm), rownames(Cm)), c(colnames(Fm), colnames(Cm))))
    # Hm = cbind(Hm, matrix(as.numeric(rownames(Cm) == "Ct_0"), nrow = 1))
    # colnames(Hm) = rownames(Fm)
    # Qm = as.matrix(Matrix::bdiag(Qm, diag(as.numeric(rownames(Cm) %in% c("Ct_0", "et_0"))*par["sig_c"]^2)))
    # colnames(Qm) = rownames(Qm) = rownames(Fm)
  }
  
  #Define the seasonal component
  if(grepl("seasonal", decomp)) {
    jiter = as.numeric(unique(gsub("[[:alpha:]]|_", "", names(par)[grepl("sig_s", names(par))])))
    for(j in jiter) {
      colnames = colnames(Fm)
      rownames = rownames(Fm)
      Sm = rbind(c(cos(2*pi*1/j), sin(2*pi*1/j)),
                 c(-sin(2*pi*1/j), cos(2*pi*1/j)))
      Fm = as.matrix(Matrix::bdiag(Fm, Sm))
      colnames(Fm) = c(colnames, paste0("St", j, "_1"), paste0("Sts", j, "_1"))
      rownames(Fm) = c(rownames, paste0("St", j, "_0"), paste0("Sts", j, "_0"))
    }
    Hm = cbind(Hm, matrix(rep(c(1, 0), length(jiter)), nrow = 1))
    colnames(Hm) = rownames(Fm)
    Qm = as.matrix(Matrix::bdiag(Qm, diag(unlist(lapply(jiter, function(j){rep(par[paste0("sig_s", j)]^2, 2)})))))
    colnames(Qm) = rownames(Qm) = rownames(Fm)
  }
  
  #Transition equation intercept matrix
  Dm = matrix(0, nrow = nrow(Fm), dimnames = list(rownames(Fm), NULL))
  if(trend == "random-walk-drift"){
    Dm[rownames(Dm) == "Dt_0", ] = par["d"]
  }
  
  #Observation equation intercept matrix
  Am = matrix(0, nrow = 1, ncol = 1)
  
  #Observation equation error covariance matrix
  Rm = matrix(par["sig_e"]^2, nrow = 1, ncol = 1)
  
  #Initial guess for unobserved vector
  if(any(rownames(Fm) %in% names(par))){
    B0 = par[rownames(Fm)]
  }else if(!is.null(init)) {
    B0 = init[["B0"]]
  }else{
    B0 = rep(0, nrow(Fm))
    names(B0) = rownames(Fm)
    B0[grepl("Tt", names(B0))] = yt[1, 1]
  }
  
  #Initial guess for variance of the unobserved vector
  if("P0" %in% names(par)){
    P0 = diag(par["P0"]^2, nrow = nrow(Fm), ncol = nrow(Fm))
    rownames(P0) = colnames(P0) = rownames(Fm)
  }else if(!is.null(init)){
    P0 = init[["P0"]]
  }else{
    P0 = diag(ifelse(ncol(yt) > 2, stats::var(c(yt), na.rm = TRUE), 100), nrow = nrow(Fm), ncol = nrow(Fm))
    rownames(P0) = colnames(P0) = rownames(Fm)
  }
  return(list(B0 = B0, P0 = P0, Am = Am, Dm = Dm, Hm = Hm, Fm = Fm, Rm = Rm, Qm = Qm, beta = beta))
}
