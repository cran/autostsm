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
    col_classes = unique(unlist(lapply(colnames(y), function(x){class(y[, c(x), with = FALSE][[1]])})))
    if(ncol(y) != 2){
      stop("y must have two columns")
    }else if(length(col_classes) < 2){
      stop("y must habe a date and numeric column")
    }else if(!all(col_classes %in% c("Date", "yearmon", "POSIXct", "POSIXt", "POSIXlt", "numeric"))){
      stop("y must have a date and numeric column")
    }
  }
}

#' Data check for input exo
#' 
#' Checks for proper input of the table exo
#' @param exo matrix of exogenous data
#' @param y input data y
#' @import data.table
#' @return none
stsm_check_exo = function(exo, y){
  if(!is.null(exo)){
    if(!(is.data.frame(exo) | is.data.table(exo))){
      stop("exo must be a data frame or data table")
    }else{
      exo = as.data.table(exo)
      col_classes = unique(unlist(lapply(colnames(exo), function(x){class(exo[, c(x), with = FALSE][[1]])})))
      if(nrow(exo) != nrow(y)){
        stop("exo and y must have the same number of rows")
      }else if(length(col_classes) < 2){
        stop("exo must habe a date and numeric column")
      }else if(!all(col_classes %in% c("Date", "yearmon", "POSIXct", "POSIXt", "POSIXlt", "numeric"))){
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
#' @param exo_obs exogenous observation data
#' @param exo_state exogenous state data
#' @param dates dates vector
#' @param range range of data to include
#' @import data.table
#' @return a data table
stsm_format_exo = function(exo_obs, exo_state, dates, range){
  if(!is.null(exo_obs)){
    col_classes = unlist(lapply(colnames(exo_obs), function(x){class(exo_obs[, c(x), with = FALSE][[1]])}))
    date_obs = colnames(exo_obs)[col_classes %in% c("Date", "yearmon", "POSIXct", "POSIXt", "POSIXlt")]
    colnames(exo_obs)[colnames(exo_obs) != date_obs] = paste0("obs.", colnames(exo_obs)[colnames(exo_obs) != date_obs])
  }
  if(!is.null(exo_state)){
    col_classes = unlist(lapply(colnames(exo_state), function(x){class(exo_state[, c(x), with = FALSE][[1]])}))
    date_state = colnames(exo_state)[col_classes %in% c("Date", "yearmon", "POSIXct", "POSIXt", "POSIXlt")]
    colnames(exo_state)[colnames(exo_state) != date_state] = paste0("state.", colnames(exo_state)[colnames(exo_state) != date_state])
  }
  if(!is.null(exo_obs) & !is.null(exo_state)){
    col_classes = unlist(lapply(colnames(exo_obs), function(x){class(exo_obs[, c(x), with = FALSE][[1]])}))
    date_obs = colnames(exo_obs)[col_classes %in% c("Date", "yearmon", "POSIXct", "POSIXt", "POSIXlt")]
    col_classes = unlist(lapply(colnames(exo_state), function(x){class(exo_state[, c(x), with = FALSE][[1]])}))
    date_state = colnames(exo_state)[col_classes %in% c("Date", "yearmon", "POSIXct", "POSIXt", "POSIXlt")]
    exo = merge(exo_obs, exo_state, by.x = date_obs, by.y = date_state, all = TRUE)
  }else if(!is.null(exo_obs)){
    exo = exo_obs 
  }else if(!is.null(exo_state)){
    exo = exo_state
  }else{
    exo = NULL
  }
  
  if(!is.null(exo)){
    exo = as.data.table(exo)
    exo = merge.data.table(exo, data.table(date = dates), all = TRUE)
    exo = exo[, names(which(unlist(exo[, lapply(.SD, is.numeric)]))), with = FALSE]
  }
  return(exo)
}

