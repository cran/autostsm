#' Cox-Stuart Test
#' 
#' Taken from the `tsutils` package. Performs the Cox-Stuart test for trend, deviation, or dispersion
#' @param y input data
#' @param type Type of test: "trend", "deviation", or "dispersion"
#' If type = "trend", test for changes in trend
#' If type = "deviation", test for changes in deviation
#' If type = "dispersion", test for changes in dispersion (range)
#' @param sig_level Significance level to determine statistically significant seasonal frequencies
#' @return list describing the results
stsm_coxstuart = function (y, type = c("trend", "deviation", "dispersion"), sig_level = 0.01){
  split = function(y){
    # Helper function
    n = length(y)
    
    # Find K according to Cox Stuart guidelines
    if(n < 48){
      k = 2
    }else if(n < 64){
      k = 3
    }else if(n < 90){
      k = 4
    }else{
      k = 5
    }
    
    # Split time series to subsamples
    rem = n %% k
    kk = floor(n/k)
    idx = array(1:(n-rem),c(k,kk))
    idx[,(round(kk/2)+1):kk] = idx[,(round(kk/2)+1):kk] + rem
    z = array(y[idx], c(k,kk))
    return(z)
  }
  
  #Replicated from tsutils::coxstuart
  type = match.arg(type, c("trend", "deviation", "dispersion"))
  switch(type, trend = {
    z = y
  }, deviation = {
    z = split(y)
    z = apply(z, 2, stats::sd)
  }, dispersion = {
    z = split(y)
    z = apply(z, 2, function(x) {
      diff(range(x))
    })
  })
  n = length(z)
  k = ceiling(n/2)
  idx1 = 1:(n - k)
  idx2 = (1 + k):n
  pair = z[idx1] - z[idx2]
  stat = min(c(sum(pair > 0), sum(pair < 0)))
  if(sum(pair) != 0) {
    p = stats::pbinom(stat, k, 0.5)
  }else{
    p = 1
  }
  if(p <= sig_level/2){
    H = 1
    txt = "H1: There is change (upwards or downwards)."
  }else {
    H = 0
    txt = "H0: There is no change present."
  }
  return(list(H = H, p.value = p, Htxt = txt))
}


