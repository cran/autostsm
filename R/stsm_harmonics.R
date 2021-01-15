#' Get harmonics for the seasons
#'
#' @param freq Frequency of the data
#' @param decomp Decomposition model ("tend-cycle-seasonal", "trend-seasonal", "trend-cycle", "trend-noise")
#' @param seasons The seasonal periods: i.e. c(365.25, 7 if yearly and weekly seasonality). Default is NULL and will be estimated via wavelet analysis.
#' @param wavelet Whether seasons were detected using wavelet analysis
#' @import data.table
#' @return list containing the initial values for the Kalman filter
stsm_harmonics = function(freq, decomp = "", seasons = NULL, wavelet = FALSE){
  #Define the default harmonics allowed
  #yearly, semi-yearly, quarterly, every other month, monthly, twice per month, weekly, the work week, every other day
  harmonic_defaults = c(1, 2, 4, 6, 12, 24, 365.25/7, 365.25/5, 365.25/2)
  
  #Define a function to build the harmonics
  get_harmonics = function(seasons, freq, harmonic_defaults){
    return(sort(unique(unlist(lapply(seasons, function(s){
      h = 1:floor(s/2)
      if(h[length(h)] == floor(s/2)){
        h[length(h)] = s/2
      }
      for(i in harmonic_defaults){
        if(any(h == floor(i))){
          h[h == floor(i)] = i
        }
      }
      return(h[h %in% harmonic_defaults & h <= freq/2])
    })))))
  }
  
  #Get seasonal harmonics if not already set
  if(grepl("seasonal", decomp)){
    #Define the harmonics
    if(wavelet == FALSE){
      if(is.null(seasons)){
        harmonics = harmonic_defaults[harmonic_defaults <= freq/2]
      }else{
        harmonics = get_harmonics(seasons, freq, harmonic_defaults)
      }
    }else{
      harmonics = freq/seasons
    }
    return(sort(unique(harmonics)))
  }else{
    return(NULL)
  }
}
  