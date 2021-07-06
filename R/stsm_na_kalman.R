#' Missing Value Imputation by Kalman Smoothing and State Space Models
#' 
#' Simplified version taken from the `imputeTS` package.
#' Uses Kalman Smoothing on structural time series models for imputation. It uses "StructTS" to build a
#' "basic structural model" if the frequency of y is greater than 1. Otherwise, it uses a local trend model.
#' @param y Univariate time series
stsm_na_kalman = function(y){
  missindx = is.na(y)
  if(!any(is.na(y))){
    return(y)
  }
  if(sum(!missindx) < 3){
    y[missindx] = y[missindx - 1]
  }
  mod = stats::StructTS(y)$model0
  kal = stats::KalmanSmooth(y, mod, nit = -1)
  erg = kal$smooth
  karima = erg[missindx, , drop = FALSE] %*% as.matrix(mod$Z)
  y[missindx] = karima
  return(y)
}