#' Set the inequality constraints for estimation
#'
#' Inequality constraints: ineqA %*% par + ineqB > 0 => ineqA %*% par > -ineqB
#' @param prior A data table created by stsm_prior
#' @param par parameter values for the state space model
#' @param freq Frequency of the data
#' @param unconstrained Whether to remove inequality constraints on the trend during estimation
#' @param det_obs Set the observation equation error variance to 0 (deterministic observation equation)
#' @param det_trend Set the trend error variance to 0 (deterministic trend)
#' @param det_seas Set the seasonality error variances to 0 (deterministic seasonality)
#' @param det_cycle Set the cycle error variance to 0 (deterministic cycle)
#' @param det_drift Set the drift error variance to 0 (deterministic drift)
#' @param saturating_growth Force the growth rate to converge to 0 in the long term 
#' @import data.table
#' @return list containing the initial values for the Kalman filter
stsm_constraints = function(prior, par, freq, unconstrained, det_trend, det_drift, det_cycle, det_seas, det_obs, saturating_growth){
  #All sigma parameters must be positive unless they are deterministic
  ineqA = matrix(0, nrow = length(par[grepl("sig_", names(par))]) - 
                   sum(c(det_trend, det_drift, det_cycle, det_seas*sum(grepl("sig_s", names(par))), det_obs)),
                 ncol = length(par),
                 dimnames = list(NULL, names(par)))
  nr = 1
  for(j in colnames(ineqA)[grepl("sig_", colnames(ineqA))]){
    if(j == "sig_t"){
      det_test = (det_trend == FALSE)
    }else if(j == "sig_d"){
      det_test = (det_drift == FALSE)
    }else if(j == "sig_c"){
      det_test = (det_cycle == FALSE)
    }else if(grepl("sig_s", j)){
      det_test = (det_seas == FALSE)
    }else if(j == "sig_e"){
      det_test = (det_obs == FALSE)
    }
    if(det_test == TRUE){
      ineqA[nr, j] = 1
      nr = nr + 1
    }
  }
  ineqB = matrix(0, nrow = nrow(ineqA), ncol = 1)
  constraints = list(ineqA = ineqA, ineqB = ineqB)
  
  if("d" %in% names(par) & saturating_growth == FALSE){
    #The drift constant must be the sign of the prior
    ineqA = matrix(0, nrow = 1, ncol = length(par), 
                   dimnames = list(NULL, names(par)))
    ineqA[, "d"] = sign(par["d"]) #constrain do to be sign of the prior
    ineqB = matrix(0, nrow = nrow(ineqA))
    constraints[["ineqA"]] = rbind(constraints[["ineqA"]], ineqA)
    constraints[["ineqB"]] = rbind(constraints[["ineqB"]], ineqB)
  }
  if("phi_d" %in% names(par)){
    #The drift must be stationary
    ineqA = matrix(0, nrow = 2, ncol = length(par), 
                   dimnames = list(NULL, names(par)))
    ineqA[, "phi_d"] = c(1, -1)
    ineqB = matrix(c(0, 1), ncol = 1)
    constraints[["ineqA"]] = rbind(constraints[["ineqA"]], ineqA)
    constraints[["ineqB"]] = rbind(constraints[["ineqB"]], ineqB)
    if(abs(par["phi_d"]) >= 1){
      par["phi_d"] = 0.95
    }
  }
  if("lambda" %in% names(par)){
    #The periodicity of the cycle parameter must be between 0 and 2*pi/(2.5*freq)
    ineqA = matrix(0, nrow = 2, ncol = length(par),
                   dimnames = list(NULL, names(par)))
    ineqA[, "lambda"] = c(1, -1)
    ineqB = matrix(c(0, (2*pi)/(2.5*freq)), ncol = 1)
    constraints[["ineqA"]] = rbind(constraints[["ineqA"]], ineqA)
    constraints[["ineqB"]] = rbind(constraints[["ineqB"]], ineqB)
    if(par["lambda"] >= (2*pi)/(2.5*freq)){
      par["lambda"] = pi/(2.5*freq)
    }
  }
  if("phi_c" %in% names(par)){
    #Damped cycle must be stationary
    ineqA = matrix(0, nrow = 2, ncol = length(par),
                   dimnames = list(NULL, names(par)))
    ineqA[, "phi_c"] = c(1, -1)
    ineqB = matrix(c(0, 1), ncol = 1)
    constraints[["ineqA"]] = rbind(constraints[["ineqA"]], ineqA)
    constraints[["ineqB"]] = rbind(constraints[["ineqB"]], ineqB)
    if(abs(par["phi_c"]) >= 1){
      par["phi_c"] = 0.95
    }
  }
  if("phi_c.1" %in% names(par)){
    #AR component must be stationary
    ineqA = matrix(0, nrow = 2, ncol = length(par),
                   dimnames = list(NULL, names(par)))
    ineqA[1, grepl("phi_c\\.\\d+", colnames(ineqA))] = 1
    ineqA[2, grepl("phi_c\\.\\d+", colnames(ineqA))] = -1
    ineqB = matrix(1, nrow = 2, ncol = 1)
    constraints[["ineqA"]] = rbind(constraints[["ineqA"]], ineqA)
    constraints[["ineqB"]] = rbind(constraints[["ineqB"]], ineqB)
    
    test = (ineqA[, grepl("phi_c\\.\\d+", colnames(ineqA))] %*% matrix(par[grepl("phi_c\\.\\d+", names(par))], ncol = 1) <= -ineqB)
    if(any(test == TRUE)){
      vals = par[grepl("phi_c\\.\\d+", names(par))]
      adj = (sum(vals) - 0.99)/length(vals)
      par[grepl("phi_c\\.\\d+", names(par))] = sign(vals)*(abs(vals) - adj)
    }
  }
  if("theta_c.1" %in% names(par)){
    #MA component must be stationary
    ineqA = matrix(0, nrow = 2, ncol = length(par),
                   dimnames = list(NULL, names(par)))
    ineqA[1, grepl("theta_c\\.\\d+", colnames(ineqA))] = 1
    ineqA[2, grepl("theta_c\\.\\d+", colnames(ineqA))] = -1
    ineqB = matrix(1, nrow = 2, ncol = 1)
    constraints[["ineqA"]] = rbind(constraints[["ineqA"]], ineqA)
    constraints[["ineqB"]] = rbind(constraints[["ineqB"]], ineqB)
    
    test = (ineqA[, grepl("theta_c\\.\\d+", colnames(ineqA))] %*% matrix(par[grepl("theta_c\\.\\d+", names(par))], ncol = 1) <= -ineqB)
    if(any(test == TRUE)){
      vals = par[grepl("theta_c\\.\\d+", names(par))]
      adj = (sum(vals) - 0.99)/length(vals)
      par[grepl("theta_c\\.\\d+", names(par))] = sign(vals)*(abs(vals) - adj)
    }
  }
  if(("sig_d" %in% names(par) | "sig_t" %in% names(par)) & unconstrained == FALSE & saturating_growth == FALSE){
    #The variance of the trend must be the smallest variance component
    ineqA = matrix(0, nrow = ("sig_e" %in% names(par) & det_obs == FALSE) +
                     ("sig_c" %in% names(par) & det_cycle == FALSE) +
                     (any(grepl("sig_s", names(par))) & det_seas == FALSE),
                   ncol = length(par), dimnames = list(NULL, names(par)))
    nr = 1
    if("sig_e" %in% names(par) & det_obs == FALSE){
      ineqA[nr, c(colnames(ineqA)[colnames(ineqA) %in% c("sig_d", "sig_t")], "sig_e")] = c(rep(-1, sum(colnames(ineqA) %in% c("sig_d", "sig_t"))), 1)
      nr = nr + 1
    }
    if("sig_c" %in% names(par) & det_cycle == FALSE){
      ineqA[nr, c(colnames(ineqA)[colnames(ineqA) %in% c("sig_d", "sig_t")], "sig_c")] = c(rep(-1, sum(colnames(ineqA) %in% c("sig_d", "sig_t"))), 1)
      nr = nr + 1
    }
    if(any(grepl("sig_s", names(par))) & det_seas == FALSE){
      ineqA[nr, c(colnames(ineqA)[colnames(ineqA) %in% c("sig_d", "sig_t")], colnames(ineqA)[grepl("sig_s", colnames(ineqA))])] = c(rep(-1, sum(colnames(ineqA) %in% c("sig_d", "sig_t"))), rep(1, sum(grepl("sig_s", colnames(ineqA)))))
    }
    ineqB = matrix(0, nrow = nrow(ineqA), ncol = 1)
    constraints[["ineqA"]] = rbind(constraints[["ineqA"]], ineqA)
    constraints[["ineqB"]] = rbind(constraints[["ineqB"]], ineqB)
    test_val = suppressWarnings(min(c(ifelse("sig_e" %in% names(par) & det_obs == FALSE, par["sig_e"], NA),
                     ifelse("sig_c" %in% names(par) & det_cycle == FALSE, par["sig_c"], NA),
                     ifelse(any(grepl("sig_s", names(par))) & det_seas == FALSE, par[grepl("sig_s", names(par))], NA)), na.rm = TRUE))
    if(is.finite(test_val)){
      if(sum(par[names(par) %in% c("sig_d", "sig_t")]) >= test_val){
        if("sig_t" %in% names(par) & det_trend == FALSE){
          par["sig_t"] = test_val/2.1
        }
        if("sig_d" %in% names(par) & det_drift == FALSE){
          par["sig_d"] = test_val/2.1
        }
      }
    }
  }
  
  if(nrow(ineqA) == 0){
    constraints = NULL
  }
  
  return(list(par = par, constraints = constraints))
}
