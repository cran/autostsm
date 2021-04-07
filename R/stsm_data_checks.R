#' Data check for input y
#'
#' Checks for proper input of the table y
#' @param y input data y
#' @import data.table
#' @return none
stsm_check_y = function(y){
  if(!(is.data.frame(y) | is.data.table(y) | stats::is.ts(y))){
    stop("y must be a data frame or data table")
  }else if(is.data.frame(y) | is.data.table(y)){
    y = as.data.table(y)
    if(ncol(y) != 2){
      stop("y must have two columns")
    }else if(!all(unlist(lapply(colnames(y), function(x){class(y[, c(x), with = FALSE][[1]])})) %in% 
                  c("Date", "yearmon", "POSIXct", "POSIXt", "POSIXlt", "numeric"))){
      stop("y must have a date and numeric column")
    }
  }
}

#' Data check for input exo
#' 
#' Checks for proper input of the table exo
#' @param exo exo datagenous
#' @param y input data y
#' @import data.table
#' @return none
stsm_check_exo = function(exo, y){
  if(!is.null(exo)){
    if(!(is.data.frame(exo) | is.data.table(exo))){
      stop("exo must be a data frame or data table")
    }else{
      exo = as.data.table(exo)
      if(nrow(exo) != nrow(y)){
        stop("exo and y must have the same number of rows")
      }else if(!all(unlist(lapply(colnames(exo), function(x){class(exo[, c(x), with = FALSE][[1]])})) %in% 
                    c("Date", "yearmon", "POSIXct", "POSIXt", "POSIXlt", "numeric"))){
        stop("exo must have a date and numeric column")
      }
    }
  }
}

#' Data check for input exo.fc
#' 
#' Checks for proper input of the table exo.fc
#' @param exo.fc exogenous forecast data
#' @param n.ahead forecast periods
#' @import data.table
#' @return none
stsm_check_exo_fc = function(exo.fc, n.ahead){
  if(!is.null(exo.fc)){
    if(!(is.data.frame(exo.fc) | is.data.table(exo.fc))){
      stop("exo must be a data frame or data table")
    }else{
      exo.fc = as.data.table(exo.fc)
      if(nrow(exo.fc) != n.ahead){
        stop("exo.fc must have the same number of rows as n.ahead")  
      }else if(!all(unlist(lapply(colnames(exo.fc), function(x){class(exo.fc[, c(x), with = FALSE][[1]])})) %in% 
                   c("Date", "yearmon", "POSIXct", "POSIXt", "POSIXlt", "numeric", "integer"))){
        stop("exo.fc must have a date and numeric column")
      }
    }
  }
}

#' Format exo
#' 
#' Format the exo table
#' @param exo exogenous data
#' @param dates dates vector
#' @param range range of data to include
#' @import data.table
#' @return a data table
stsm_format_exo = function(exo, dates, range){
  if(!is.null(exo)){
    exo = as.data.table(exo)
    exo = exo[range[1]:range[length(range)], ]
    exo = merge.data.table(exo, data.table(date = dates), all = TRUE)
    exo = exo[, names(which(unlist(exo[, lapply(.SD, is.numeric)]))), with = FALSE]
  }
  return(exo)
}

