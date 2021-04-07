#' Detect factors
#'
#' @param sig_level Significance level to determine statistically significant seasonal frequencies
#' @param models list of univariate model structures
#' @import data.table
#' @return 
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
#' seasonality = stsm_detect_seasonality(y = NA000334Q$y, freq = 4)
#' }
#' #NEED TO ADD at_export
# stsm_detect_factors = function(models, sig_level = 0.01){
# 
#   #For variables that have a unit root, check for cointegration among the trends
#   ur.vars = names(which(sapply(models, function(x){x[["ndiffs"]] > 0})))
#   if(length(ur.vars) > 0){
#     dt = as.data.table(lapply(ur.vars, function(x){
#       return(rowSums(models[[x]][["prior"]][, c("trend", "cycle", "remainder")]))
#     }))
#     #Select the appropriate lag length
#     K = table(vars::VARselect(dt, type = "const")$selection)
#     K = as.numeric(names(K))[which.max(K)][1]
#     
#     #Perform the Johansen cointegration test
#     vecm = urca::ca.jo(dt, type = "eigen", ecdet = "const", K = K, spec = "transitory")
#     coint.test = data.table(vecm@cval, keep.rownames = TRUE)[, "stat" := vecm@teststat]
#     coint.test[, "rn" := as.numeric(gsub("[[:alpha:]]|[[:punct:]]", "", rn))]
#     col = paste0(ifelse(sig_level <= 0.01, 1, ifelse(sig_level <= 0.05, 5, 10)), "pct")
#     coint.test = coint.test[order(rn), ]
#     coint.test[, "sig" := stat >= eval(parse(text = paste0("`", col, "`")))]
#     for(coint in coint.test$rn){
#       if(coint.test[rn == coint, ]$sig == FALSE){
#         break
#       }else if(coint == max(coint.test$rn) & coint.test[rn == coint, ]$sig == TRUE){
#         coint = 0
#       }
#     }
#     if(coint > 0){
#       vecm2 = urca::cajorls(vecm, r = coint)
#       rlm = coef(vecm2$rlm)
#       #Check that the ect coefficients are negative
#       #If any are positive, correct the sign #remove those as cointegrating relationships
#       for(i in 1:coint){
#         if(rlm[paste0("ect", i), i] > 0){
#           rlm[paste0("ect", i), i] = -1*rlm[paste0("ect", i), i]
#         }
#       }
#       # coint = coint - sum(sign(sapply(1:coint, function(x){rlm[paste0("ect", x), x]})) == 1)
#       # if(coint > 0){
#         # vecm2 = urca::cajorls(vecm, r = coint)
#         # alpha = coef(vecm2$rlm)
#         alpha = rlm
#         alpha = t(matrix(alpha[grepl("ect", rownames(alpha)), ], nrow = coint, 
#                        dimnames = list(paste0("ect", 1:coint), colnames(alpha))))
#         alpha0 = qr.Q(qr(alpha), complete = TRUE)[, (nrow(alpha) - (ncol(dt) - coint - 1)):nrow(alpha)]
#         
#         #Get the common trends and smooth them
#         common_trends = data.table(as.matrix(dt) %*% alpha0)
#         colnames(common_trends) = paste0("trend_", 1:ncol(common_trends))
#         common_trends = common_trends[, lapply(.SD, function(x){
#           stats::predict(stats::loess(x ~ t, data.table(x = x)[, "t" := 1:length(x)]))
#         })]
#         
#         #Get the common cycles
#         common_cycles = data.table(cbind(as.matrix(dt), 1) %*% vecm2$beta)
#         colnames(common_cycles) = paste0("cycle_", 1:ncol(common_cycles))
#         n_cycles = ncol(common_cycles)
#         n_trends = ncol(common_trends)
#       # }
#     }else{
#       coint = n_trends = 0
#     }
#   }else{
#     coint = n_trends = 0
#   }
#   
#   #Check for factors in the seasonally adjusted and detrended data
#   if(coint == 0){
#     dt = as.data.table(lapply(names(models), function(x){
#       return(rowSums(models[[x]][["prior"]][, c("cycle", "remainder")]))
#     }))
# 
#     #Calculate the relative eigenvalues and the p-value
#     cov = cov(scale(dt), use = "pairwise.complete.obs")
#     ev = abs(sort(eigen(cov)$values, decreasing = TRUE))
#     ev[ev <= .Machine$double.eps * max(ncol(dt)) * max(ev)] = 0
#     ev.test = data.table(p_value = rep(0, length(ev) - 1))
#     for (i in 1:nrow(ev.test)) {
#       val = 2 * nrow(dt) * log((ev[i] + ev[i + 1])/(2 * (ev[i] * ev[i + 1])^0.5))
#       ev.test$p_value[i] = stats::pchisq(q = val, df = 2, lower.tail = FALSE)
#     }
#     ev.test = ev.test[!is.na(p_value), ]
#     
#     #Perform a Johansen-style test for significant eigenvalues to determine the number of factors
#     ev.test[, "sig" := p_value <= 0.0001][, "rn" := 1:.N]
#     for(n_cycles in ev.test$rn){
#       if(ev.test[rn == n_cycles, ]$sig == FALSE){
#         n_cycles = n_cycles - 1
#         break
#       }
#     }
#     
#     #Retrieve the factors
#     if(n_cycles > 0 & n_cycles < ncol(dt)){
#       common_cycles = suppressMessages(suppressWarnings(data.table(psych::fa(dt, nfactors = n_cycles, rotate = "none")$scores)))
#       colnames(common_cycles) = paste0("cycle_", 1:ncol(common_cycles))
#     }else{
#       n_cycles = 0
#     }
#   }
#   
#   #Redefine the model priors
#   is_drift = any(sapply(models, function(x){grepl("drift", x$trend)}))
#   for(i in 1:length(models)){
#     if(sum(grepl("trend", colnames(models[[i]][["prior"]]))) > 0 & n_trends > 0){
#       models[[i]][["prior"]] = cbind(models[[i]][["prior"]][, "trend" := NULL], common_trends)
#       if(is_drift == TRUE){
#         models[[i]][["prior"]] = models[[i]][["prior"]][, "drift" := NULL]
#         models[[i]][["prior"]][, paste0("drift_", 1:n_trends) := lapply(.SD, function(x){
#           x - shift(x, type = "lag", n = 1)
#         }), .SDcols = paste0("trend_", 1:ncol(common_trends))]
#         models[[i]][["prior"]][1, paste0("drift_", 1:n_trends) := 0]
#       }
#     }
#     if(sum(grepl("cycle", colnames(models[[i]][["prior"]]))) > 0 & n_cycles > 0){
#       models[[i]][["prior"]] = cbind(models[[i]][["prior"]][, "cycle" := NULL], common_cycles)
#     }
#     
#     lm = lm(z ~ . - 1, data = cbind(z = dt[, c(i), with = FALSE][[1]], models[[i]][["prior"]]))
#     models[[i]][["prior"]][, "remainder" := lm$residuals]
#   }
#     
#   return(list(models = models, n_trends = n_trends, n_cycles = n_cycles))
# }
#  
