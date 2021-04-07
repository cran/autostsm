#' Create a multivariate state space model from a list of univariate state space models
#' 
#' @param ssms list of univariate state space models created from ssm
#' @param uc either "common" or "separate"
#' "common" will create a common unobserved component for each variable
#' "separate" will create a separate unobserved component for each variable
#' @param common_ucs named character vector corresponding the common unobserved components across variables.
#' The default is NULL which will set all unobserved components to be common.
#' @return list of matrices
#' #NEED TO ADD at_export
# stsm_multivar_ssm = function(ssms, uc = "separate", common_ucs = NULL){
#   
#   #uc = "separate"
#   #common_ucs = NULL
#   #ssms = lapply(1:length(models), function(x){
#   # stsm_ssm(models[[x]]$par, matrix(yt[x, ], nrow = 1), models[[x]]$decomp, models[[x]]$trend, models[[x]]$init)
#   #})
#   #names(ssms) = names(models)
#   if(is.character(uc)){
#     uc = tolower(uc)
#     if(!tolower(uc) %in% c("common", "separate")){
#       stop("uc must be 'common' or 'separate'")
#     }
#   }else{
#     stop("uc must be a string, either 'common' or 'separate'")
#   }
#   if(!is.null(common_ucs) & uc == "common"){
#     if(!all(unlist(lapply(ssms, function(x){all(common_ucs %in% unique(gsub("[[:digit:]]|\\.", "", dimnames(x[["Dm"]])[[1]])))})))){
#       stop("common_ucs must be present in the rownames of Dm, Fm, Qm, and B0 across all variables.")
#     }
#   }else if(is.null(common_ucs) & uc == "common"){
#     if(length(unique(lapply(ssms, function(x){dimnames(x[["Dm"]])[[1]]}))) != 1){
#       stop("All rownames of Dm, Fm, Qm, and B0 must be the same across all variables for a full common uc model.")
#     }
#   }
#   
#   #Stack Am
#   Am = lapply(names(ssms), function(x){
#     matrix(ssms[[x]][["Am"]], nrow = nrow(ssms[[x]][["Am"]]), ncol = ncol(ssms[[x]][["Am"]]), 
#            dimnames = list(x, colnames(ssms[[x]][["Am"]])))
#   })
#   Am = do.call("rbind", Am)
#   names(Am) = names(ssms)
#   
#   #Stack Hm
#   Hm = lapply(names(ssms), function(x){
#     matrix(ssms[[x]][["Hm"]], nrow = nrow(ssms[[x]][["Hm"]]), ncol = ncol(ssms[[x]][["Hm"]]), 
#            dimnames = list(x, colnames(ssms[[x]][["Hm"]])))
#   })
#   names(Hm) = names(ssms)
#   if(uc == "common"){
#     if(is.null(common_ucs)){
#       Hm = as.matrix(do.call(function(x, y){rbind(x, y, use.names = TRUE, fill = TRUE)}, lapply(names(Hm), function(x){
#         data.table(Hm[[x]])
#       })), rownames = names(Hm))
#     }else{
#       Hm = as.matrix(do.call(function(x, y){rbind(x, y, use.names = TRUE, fill = TRUE)}, lapply(names(Hm), function(x){
#         temp = data.table(Hm[[x]])
#         eqn = ssms[[x]]$Fm[paste0(common_ucs, ".0"),]
#         eqn = unique(gsub("[[:digit:]]|\\.", "", names(eqn)[eqn != 0]))
#         colnames(temp)[colnames(temp) %in% eqn] = paste0(x, ":", colnames(temp)[colnames(temp) %in% eqn])
#         return(temp)
#       })), rownames = names(Hm))
#     }
#     Hm[is.na(Hm)] = 0
#     Hvars = unique(sapply(colnames(Hm), function(x){strsplit(x, "\\.")[[1]][1]}))
#     col_order = unname(unlist(sapply(Hvars, function(x){sort(colnames(Hm)[grepl(x, colnames(Hm))])})))
#     Hm = matrix(Hm[, col_order], nrow = nrow(Hm), ncol = ncol(Hm), dimnames = list(rownames(Hm), col_order))
#   }else{
#     dimnames = lapply(1:length(Hm), function(x){list(rownames(Hm[[x]]), 
#                                                      paste0(names(Hm)[x], ":", colnames(Hm[[x]])))})
#     dimnames = list(unlist(lapply(dimnames, function(x){x[[1]]})), 
#                     unlist(lapply(dimnames, function(x){x[[2]]})))
#     Hm = as.matrix(bdiag(Hm))
#     dimnames(Hm) = dimnames
#   }
#   
#   #Stack Rm
#   Rm = lapply(names(ssms), function(x){
#     matrix(ssms[[x]][["Rm"]], nrow = nrow(ssms[[x]][["Rm"]]), ncol = ncol(ssms[[x]][["Rm"]]), 
#            dimnames = list(x, x))
#   })
#   names(Rm) = names(ssms)
#   names = names(Rm)
#   Rm = as.matrix(Matrix::bdiag(Rm))
#   dimnames(Rm) = list(names, names)
#   
#   #Stack Dm
#   Dm = lapply(names(ssms), function(x){
#     matrix(ssms[[x]][["Dm"]], nrow = nrow(ssms[[x]][["Dm"]]), ncol = ncol(ssms[[x]][["Dm"]]), 
#            dimnames = list(rownames(ssms[[x]][["Dm"]]), x))
#   })
#   names(Dm) = names(ssms)
#   if(uc == "common"){
#     if(is.null(common_ucs)){
#       Dm = as.matrix(do.call(function(x, y){merge(x, y, by = colnames(x)[1], all = TRUE)}, lapply(names(Dm), function(x){
#         data.table(Dm[[x]], keep.rownames = TRUE)
#       })), rownames = TRUE)
#     }else{
#       Dm = as.matrix(do.call(function(x, y){merge(x, y, by = colnames(x)[1], all = TRUE)}, lapply(names(Dm), function(x){
#         temp = data.table(Dm[[x]], keep.rownames = TRUE)
#         eqn = ssms[[x]]$Fm[paste0(common_ucs, ".0"),]
#         eqn = unique(gsub("[[:digit:]]|\\.", "", names(eqn)[eqn != 0]))
#         rownames(temp)[rownames(temp) %in% eqn] = paste0(x, ":", rownames(temp)[rownames(temp) %in% eqn])
#         return(temp)
#       })), rownames = TRUE)
#     }
#     Dm = matrix(rowSums(Dm, na.rm = TRUE), ncol = 1, nrow = nrow(Dm), dimnames = list(rownames(Dm), NULL))
#     Dm[is.na(Dm)] = 0
#     Dm = matrix(Dm[colnames(Hm), ], ncol = 1, nrow = nrow(Dm), dimnames = list(colnames(Hm), NULL))
#   }else if(uc == "separate"){
#     rnames = unlist(lapply(1:length(Dm), function(x){list(paste0(names(Dm)[x], ":", rownames(Dm[[x]])))}))
#     Dm = do.call("rbind", Dm)
#     rownames(Dm) = rnames
#   }
#   
#   #Stack Fm
#   Fm = lapply(names(ssms), function(x){
#     matrix(ssms[[x]][["Fm"]], nrow = nrow(ssms[[x]][["Fm"]]), ncol = ncol(ssms[[x]][["Fm"]]), 
#            dimnames = dimnames(ssms[[x]][["Fm"]]))
#   })
#   names(Fm) = names(ssms)
#   if(uc == "common"){
#     col_order = sapply(strsplit(rownames(Dm), "\\."), function(x){
#       paste0(x[[1]], ".", as.numeric(x[2]) + 1)
#     })
#     if(is.null(common_ucs)){
#       Fm = as.matrix(do.call(function(x, y){merge(x, y[, c(colnames(x)[1], colnames(y)[!colnames(y) %in% colnames(x)]), with = FALSE], by = colnames(x)[1], all = TRUE)}, lapply(names(Fm), function(x){
#         data.table(Fm[[x]], keep.rownames = TRUE)
#       })), rownames = TRUE)
#     }else{
#       Fm = as.matrix(do.call(function(x, y){merge(x, y[, c(colnames(x)[1], colnames(y)[!colnames(y) %in% colnames(x)]), with = FALSE], by = colnames(x)[1], all = TRUE)}, lapply(names(Fm), function(x){
#         temp = data.table(Fm[[x]], keep.rownames = TRUE)
#         eqn = ssms[[x]]$Fm[paste0(common_ucs, ".0"),]
#         eqn = unique(gsub("[[:digit:]]|\\.", "", names(eqn)[eqn != 0]))
#         rownames(temp)[rownames(temp) %in% eqn] = paste0(x, ":", rownames(temp)[rownames(temp) %in% eqn])
#         return(temp)
#       })), rownames = TRUE)
#     }
#     Fm[is.na(Fm)] = 0
#     Fm = Fm[rownames(Dm), col_order]
#     Fm[rownames(Fm) %in% colnames(Fm), colnames(Fm) %in% rownames(Fm)] = diag(sum(rownames(Fm) %in% colnames(Fm)))
#   }else if(uc == "separate"){
#     dimnames = lapply(1:length(Fm), function(x){list(paste0(names(Fm)[x], ":", rownames(Fm[[x]])), 
#                                                      paste0(names(Fm)[x], ":", colnames(Fm[[x]])))})
#     dimnames = list(unlist(lapply(dimnames, function(x){x[[1]]})), 
#                     unlist(lapply(dimnames, function(x){x[[2]]})))
#     Fm = as.matrix(bdiag(Fm))
#     dimnames(Fm) = dimnames
#   }
#   
#   #Stack Qm
#   Qm = lapply(names(ssms), function(x){
#     matrix(ssms[[x]][["Qm"]], nrow = nrow(ssms[[x]][["Qm"]]), ncol = ncol(ssms[[x]][["Qm"]]), 
#            dimnames = dimnames(ssms[[x]][["Qm"]]))
#   })
#   names(Qm) = names(ssms)
#   if(uc == "common"){
#     if(is.null(common_ucs)){
#       Qm = as.matrix(do.call(function(x, y){merge(x, y[, c(colnames(x)[1], colnames(y)[!colnames(y) %in% colnames(x)]), with = FALSE], by = colnames(x)[1], all = TRUE)}, lapply(names(Qm), function(x){
#         data.table(Qm[[x]], keep.rownames = TRUE)
#       })), rownames = TRUE)
#     }else{
#       Qm = as.matrix(do.call(function(x, y){merge(x, y[, c(colnames(x)[1], colnames(y)[!colnames(y) %in% colnames(x)]), with = FALSE], by = colnames(x)[1], all = TRUE)}, lapply(names(Qm), function(x){
#         temp = data.table(Qm[[x]], keep.rownames = TRUE)
#         eqn = ssms[[x]]$Fm[paste0(common_ucs, ".0"),]
#         eqn = unique(gsub("[[:digit:]]|\\.", "", names(eqn)[eqn != 0]))
#         rownames(temp)[rownames(temp) %in% eqn] = paste0(x, ":", rownames(temp)[rownames(temp) %in% eqn])
#         return(temp)
#       })), rownames = TRUE)
#     }
#     Qm[is.na(Qm)] = 0
#     Qm = Qm[rownames(Fm), rownames(Fm)]
#   }else if(uc == "separate"){
#     dimnames = lapply(1:length(Qm), function(x){list(paste0(names(Qm)[x], ":", rownames(Qm[[x]])), 
#                                                      paste0(names(Qm)[x], ":", colnames(Qm[[x]])))})
#     dimnames = list(unlist(lapply(dimnames, function(x){x[[1]]})), 
#                     unlist(lapply(dimnames, function(x){x[[2]]})))
#     Qm = as.matrix(bdiag(Qm))
#     dimnames(Qm) = dimnames
#   }
#   
#   #Stack B0
#   B0 = lapply(names(ssms), function(x){
#     matrix(ssms[[x]][["B0"]], nrow = nrow(ssms[[x]][["B0"]]), ncol = ncol(ssms[[x]][["B0"]]), 
#            dimnames = dimnames(ssms[[x]][["B0"]]))
#   })
#   names(B0) = names(ssms)
#   if(uc == "common"){
#     if(is.null(common_ucs)){
#       B0 = as.matrix(do.call(function(x, y){merge(x, y, by = colnames(x)[1], all = TRUE)}, lapply(names(B0), function(x){
#         data.table(B0[[x]], keep.rownames = TRUE)
#       })), rownames = TRUE)
#     }else{
#       B0 = as.matrix(do.call(function(x, y){merge(x, y, by = colnames(x)[1], all = TRUE)}, lapply(names(B0), function(x){
#         temp = data.table(B0[[x]], keep.rownames = TRUE)
#         eqn = ssms[[x]]$Fm[paste0(common_ucs, ".0"),]
#         eqn = unique(gsub("[[:digit:]]|\\.", "", names(eqn)[eqn != 0]))
#         rownames(temp)[rownames(temp) %in% eqn] = paste0(x, ":", rownames(temp)[rownames(temp) %in% eqn])
#         return(temp)
#       })), rownames = TRUE)
#     }
#     B0 = matrix(rowSums(B0, na.rm = TRUE), ncol = 1, nrow = nrow(B0), dimnames = list(rownames(B0), NULL))
#     B0[is.na(B0)] = 0
#     B0 = matrix(B0[colnames(Hm), ], ncol = 1, nrow = nrow(B0), dimnames = list(colnames(Hm), NULL))
#   }else if(uc == "separate"){
#     rnames = unlist(lapply(1:length(B0), function(x){list(paste0(names(B0)[x], ":", rownames(B0[[x]])))}))
#     B0 = do.call("rbind", B0)
#     rownames(B0) = rnames  
#   }
#   
#   #Stack P0
#   P0 = lapply(names(ssms), function(x){
#     matrix(ssms[[x]][["P0"]], nrow = nrow(ssms[[x]][["P0"]]), ncol = ncol(ssms[[x]][["P0"]]), 
#            dimnames = dimnames(ssms[[x]][["P0"]]))
#   })
#   names(P0) = names(ssms)
#   if(uc == "common"){
#     if(is.null(common_ucs)){
#       P0 = as.matrix(do.call(function(x, y){merge(x, y[, c(colnames(x)[1], colnames(y)[!colnames(y) %in% colnames(x)]), with = FALSE], by = colnames(x)[1], all = TRUE)}, lapply(names(P0), function(x){
#         data.table(P0[[x]], keep.rownames = TRUE)
#       })), rownames = TRUE)
#     }else{
#       P0 = as.matrix(do.call(function(x, y){merge(x, y[, c(colnames(x)[1], colnames(y)[!colnames(y) %in% colnames(x)]), with = FALSE], by = colnames(x)[1], all = TRUE)}, lapply(names(P0), function(x){
#         temp = data.table(P0[[x]], keep.rownames = TRUE)
#         eqn = ssms[[x]]$Fm[paste0(common_ucs, ".0"),]
#         eqn = unique(gsub("[[:digit:]]|\\.", "", names(eqn)[eqn != 0]))
#         rownames(temp)[rownames(temp) %in% eqn] = paste0(x, ":", rownames(temp)[rownames(temp) %in% eqn])
#         return(temp)
#       })), rownames = TRUE)
#     }
#     P0[is.na(P0)] = 0
#     P0 = P0[rownames(Qm), rownames(Qm)]
#   }else if(uc == "separate"){
#     dimnames = lapply(1:length(P0), function(x){list(paste0(names(P0)[x], ":", rownames(P0[[x]])), 
#                                                      paste0(names(P0)[x], ":", colnames(P0[[x]])))})
#     dimnames = list(unlist(lapply(dimnames, function(x){x[[1]]})), 
#                     unlist(lapply(dimnames, function(x){x[[2]]})))
#     P0 = as.matrix(bdiag(P0))
#     dimnames(P0) = dimnames
#   }
#   
#   #Stack beta
#   beta = lapply(names(ssms), function(x){
#     matrix(ssms[[x]][["beta"]], nrow = nrow(ssms[[x]][["beta"]]), ncol = ncol(ssms[[x]][["beta"]]), 
#            dimnames = list(x, colnames(ssms[[x]][["beta"]])))
#   })
#   names(beta) = names(ssms)
#   beta = do.call("rbind", beta)
#   
#   return(list(B0 = B0, P0 = P0, Am = Am, Dm = Dm, Hm = Hm, Fm = Fm, Rm = Rm, Qm = Qm, beta = beta))
# }
