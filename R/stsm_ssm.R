#' State space model
#' 
#' Creates a state space model in list form
#' yt = H*B + e_t
#' B = F*B_{t-1} + u_t
#'
#' @param par Vector of named parameter values, includes the harmonics
#' @param yt Univariate time series of data values
#' @param freq Frequency of the data (1 (yearly), 4 (quarterly), 12 (monthly), 365.25/7 (weekly), 365.25 (daily))
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
stsm_ssm = function(par = NULL, yt = NULL, freq = NULL, decomp = NULL,
                    trend = NULL, init = NULL, model = NULL){
  if(!is.null(model)){
    par = eval(parse(text = paste0("c(", model$coef, ")")))
    freq = model$freq
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
    colnames(Fm) = "Tt1"
    rownames(Fm) = "Tt0"
    #Observation matrix
    Hm = matrix(c(1), nrow = 1)
    colnames(Hm) = rownames(Fm)
  }else if(trend %in% c("random-walk-drift", "double-random-walk")){
    #T_t = M_{t-1} + T_{t-1} + e_t, e_t ~ N(0, sig_t^2)
    #M_t = M_{t-1} + n_t, n_t ~ N(0, sig_m^2)
    #Transition matrix
    Fm = rbind(c(1, 1),
               c(0, 1))
    colnames(Fm) = c("Tt1", "Mt1")
    rownames(Fm) = c("Tt0", "Mt0")
    if(trend == "random-walk-drift"){
      Fm["Mt0", "Mt1"] = par["phi"]
    }
    #Observation matrix
    Hm = matrix(c(1, 0), nrow = 1)
    colnames(Hm) = rownames(Fm)
  }else if(trend %in% c("random-walk2")){
    #T_t = 2T_{t-1} - T_{t-2} + e_t, e_t ~ N(0, sig_t^2)
    #Transition matrix
    Fm = rbind(c(2, -1),
               c(1, 0))
    colnames(Fm) = c("Tt1", "Tt2")
    rownames(Fm) = c("Tt0", "Tt1")
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
  Qm[rownames(Qm) == "Tt0", colnames(Qm) == "Tt0"] = par["sig_t"]^2
  Qm[rownames(Qm) == "Mt0", colnames(Qm) == "Mt0"] = par["sig_m"]^2
  
  #Define the cycle component
  if(grepl("cycle", decomp)){
    Cm = par["rho"]*rbind(c(cos(par["lambda"]), sin(par["lambda"])),
                          c(-sin(par["lambda"]), cos(par["lambda"])))
    colnames(Cm) = c("Ctl", "Ctsl")
    rownames(Cm) = c("Ct", "Cts")
    Fm = matrix(Matrix::bdiag(Fm, Cm), ncol = ncol(Fm) + ncol(Cm),
                dimnames = list(c(rownames(Fm), rownames(Cm)), c(colnames(Fm), colnames(Cm))))
    Hm = cbind(Hm, matrix(c(1, 0), nrow = 1))
    colnames(Hm) = rownames(Fm)
    Qm = as.matrix(Matrix::bdiag(Qm, diag(2)*par["sig_c"]^2))
    colnames(Qm) = rownames(Qm) = rownames(Fm)
  }else if("ar1" %in% names(par) | "ma1" %in% names(par)){
    #Define the ARMA part of the noise component
    par["ar2"] = -par["ar1"]^2/4
    Cm = rbind(par[c("ar1", "ar2", "ma1")])
    Cm = matrix(rbind(Cm,
               cbind(diag(2)*c(1, 0), matrix(0, nrow = 2, ncol = 1))),
               ncol = ncol(Cm), dimnames = list(NULL, colnames(Cm)))
    rownames(Cm) = c("Ct", "Ctl", "et")
    colnames(Cm) = c("Ctl", "Ctl2", "etl")
    Fm = matrix(Matrix::bdiag(Fm, Cm), ncol = ncol(Fm) + ncol(Cm),
                dimnames = list(c(rownames(Fm), rownames(Cm)), c(colnames(Fm), colnames(Cm))))
    Hm = cbind(Hm, matrix(as.numeric(rownames(Cm) == "Ct"), nrow = 1))
    colnames(Hm) = rownames(Fm)
    Qm = as.matrix(Matrix::bdiag(Qm, diag(as.numeric(rownames(Cm) %in% c("Ct", "et"))*par["sig_c"]^2)))
    colnames(Qm) = rownames(Qm) = rownames(Fm)
  }
  
  #Define the seasonal component
  if(grepl("seasonal", decomp)) {
    jiter = as.numeric(unique(gsub("[[:alpha:]]|_", "", names(par)[grepl("sig_s", names(par))])))
    for(j in jiter) {
      colnames = colnames(Fm)
      rownames = rownames(Fm)
      Sm = rbind(c(cos(2*pi*j/freq), sin(2*pi*j/freq)),
                 c(-sin(2*pi*j/freq), cos(2*pi*j/freq)))
      Fm = as.matrix(Matrix::bdiag(Fm, Sm))
      colnames(Fm) = c(colnames, paste0("Stl", round(j)), paste0("Stls", round(j)))
      rownames(Fm) = c(rownames, paste0("St", round(j)), paste0("Sts", round(j)))
    }
    Hm = cbind(Hm, matrix(rep(c(1, 0), length(jiter)), nrow = 1))
    colnames(Hm) = rownames(Fm)
    Qm = as.matrix(Matrix::bdiag(Qm, diag(unlist(lapply(jiter, function(j){rep(par[paste0("sig_s", j)]^2, 2)})))))
    colnames(Qm) = rownames(Qm) = rownames(Fm)
  }
  
  #Transition equation intercept matrix
  Dm = matrix(0, nrow = nrow(Fm), dimnames = list(rownames(Fm), NULL))
  if(trend == "random-walk-drift"){
    Dm[rownames(Dm) == "Mt0", ] = par["d"]
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
    P0 = diag(ifelse(ncol(yt) > 2, stats::var(yt, na.rm = TRUE), 100), nrow = nrow(Fm), ncol = nrow(Fm))
    rownames(P0) = colnames(P0) = rownames(Fm)
  }
  return(list(B0 = B0, P0 = P0, At = Am, Dt = Dm, Ht = Hm, Ft = Fm, Rt = Rm, Qt = Qm, beta = beta))
}


# hp_filter = function(y, freq, trend){
#   hp_objective = function(par, y, freq, trend, lambda){
#     if(trend == "random-walk"){
#       par["sig_e"] = sqrt(lambda*par["sig_t"]^2)
#     }else if(trend == "random-walk2"){
#       par["sig_e"] = sqrt(lambda*par["sig_t"]^2)
#     }else if(trend == "random-walk-drift"){
#       par["sig_e"] = sqrt(lambda*(par["sig_t"]^2 + par["sig_m"]^2/(1 - par["phi"]^2)))
#     }else if(trend == "double-random-walk"){
#       par["sig_e"] = sqrt(lambda*(par["sig_t"]^2 + par["sig_m"]^2))
#     }
#     ssm = stsm_ssm(par, y, freq = freq, trend = trend, decomp = "trend-noise")
#     return(kalman_filter(B0 = matrix(ssm$B0, nrow = nrow(ssm$Ft), ncol = 1), P0 = diag(nrow(ssm$Qt))*100, 
#                          Dt = ssm$Dt, At = ssm$At, Ft = ssm$Ft, Ht = ssm$Ht, 
#                          Qt = ssm$Qt, Rt = ssm$Rt, yt = matrix(y, nrow = 1), 
#                          X = matrix(0, nrow = 1, ncol = length(y)), 
#                          beta = matrix(0, nrow = 1, ncol = 1))$loglik)
#   }
#   
#   #Trend and drift parameters
#   lambda = 1600*(freq/4)^4
#   if(trend == "random-walk"){
#     par = c(sig_t = sqrt(var(diff(y), na.rm = TRUE)/(1 + 2*lambda)))
#   }else if(trend == "random-walk2"){
#     par = c(sig_t = sqrt(var(diff(diff(y)), na.rm = TRUE)/(1 + 6*lambda)))
#   }else if(trend == "random-walk-drift"){
#     par = c(phi = 0.5, d = mean(diff(y), na.rm = TRUE)*(1 - 0.5))
#     par["sig_m"] = par["sig_t"] = sqrt(var(diff(y), na.rm = TRUE)/(1 + 2*lambda)*(1 - par["phi"]^2)/(2 - par["phi"]^2))
#   }else if(trend == "double-random-walk"){
#     par["sig_m"] = par["sig_t"] = sqrt(var(diff(diff(y)), na.rm = TRUE)/8)
#   }
#   
#   ineqA = matrix(0, ncol = length(par), nrow = sum(grepl("sig_", names(par))), dimnames = list(NULL, names(par)))
#   ineqA[, grepl("sig_", names(par))] = diag(sum(grepl("sig_", names(par))))
#   ineqB = matrix(0, ncol = 1, nrow = nrow(ineqA))
#   constraints = list(ineqA = ineqA, ineqB = ineqB)
#   if("d" %in% names(par) & saturating_growth == FALSE){
#     #The drift constant must be the sign of the prior
#     ineqA = matrix(0, nrow = 1, ncol = length(par), 
#                    dimnames = list(NULL, names(par)))
#     ineqA[, "d"] = sign(par["d"]) #constrain do to be sign of the prior
#     ineqB = matrix(0, nrow = nrow(ineqA))
#     constraints[["ineqA"]] = rbind(constraints[["ineqA"]], ineqA)
#     constraints[["ineqB"]] = rbind(constraints[["ineqB"]], ineqB)
#   }
#   if("phi" %in% names(par)){
#     #The drift must be stationary
#     ineqA = matrix(0, nrow = 2, ncol = length(par), 
#                    dimnames = list(NULL, names(par)))
#     ineqA[, "phi"] = c(1, -1)
#     ineqB = matrix(c(0, 1), ncol = 1)
#     constraints[["ineqA"]] = rbind(constraints[["ineqA"]], ineqA)
#     constraints[["ineqB"]] = rbind(constraints[["ineqB"]], ineqB)
#   }
#   
#   hp_out = maxLik::maxLik(logLik = hp_objective,
#                           start = par, method = "BFGS", constraints = constraints,
#                           finalHessian = FALSE, hess = NULL, control = list(printLevel = 0, iterlim = 10000), 
#                           y = y, freq = freq, trend = trend, lambda = lambda)
#   if(trend == "random-walk"){
#     hp_out$estimate["sig_e"] = sqrt(lambda*hp_out$estimate["sig_t"]^2)
#   }else if(trend == "random-walk2"){
#     hp_out$estimate["sig_e"] = sqrt(lambda*hp_out$estimate["sig_t"]^2)
#   }else if(trend == "random-walk-drift"){
#     hp_out$estimate["sig_e"] = sqrt(lambda*(hp_out$estimate["sig_t"]^2 + hp_out$estimate["sig_m"]^2/(1 - hp_out$estimate["phi"]^2)))
#   }else if(trend == "double-random-walk"){
#     hp_out$estimate["sig_e"] = sqrt(lambda*(hp_out$estimate["sig_t"]^2 + hp_out$estimate["sig_m"]^2))
#   }
#   
#   ssm = stsm_ssm(hp_out$estimate, y, freq = freq, trend = trend, decomp = "trend-noise")
#   filter = kalman_filter(B0 = matrix(ssm$B0, nrow = 2, ncol = 1), P0 = diag(nrow(ssm$Qt))*100, 
#                          Dt = ssm$Dt, At = ssm$At, Ft = ssm$Ft, Ht = ssm$Ht, 
#                          Qt = ssm$Qt, Rt = ssm$Rt, yt = matrix(y, nrow = 1), 
#                          X = matrix(0, nrow = 1, ncol = length(y)), 
#                          beta = matrix(0, nrow = 1, ncol = 1))
#   smooth = kalman_smoother(B_tl = filter$B_tl, B_tt = filter$B_tt, 
#                            P_tl = filter$P_tl, P_tt = filter$P_tt, 
#                            Ft = ssm$Ft)
#   # ts.plot(ts(y), ts(filter$B_tt[1, ]), col = c("black", "red"))
#   # ts.plot(ts(y), ts(smooth$B_tt[1, ]), col = c("black", "red"))
#   return(c(smooth$B_tt[1, ]))
# }