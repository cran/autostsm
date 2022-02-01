
# autostsm

<!-- badges: start -->
<!-- badges: end -->

`autostsm` (Automatic Structural Time Series Model) is designed to automatically detect the appropriate decomposition for a univariate time series into trend, seasonal, and cycle components using a state space approach. The package also has the functionality perform structural interpolation: i.e. interpolated a quarterly series into a monthly series using either end of period, period average, or period sum methods. The unobserved components are estimated with the Kalman filter and all model parameters are estimated using maximum likelihood. The Kalman filter and smoother are implemented in `Rcpp` so it is reasonably fast. This package is designed to handle time series of any frequency (standard or not), series with multiple seasonalities, seasonalities with fractional periodicities, missing data imputation, and irregularly spaced dates (i.e. daily data with missing data due to weekends and holidays, etc.).

## Installation

You can install the development version of `autostsm` like so:

```{r, eval = FALSE}
devtools::install_git("https://bitbucket.org/ajhubb/autostsm.git")
```

## Example

This is a basic example which shows you how to use the package on some example data. See the vignette for more examples:

```{r, eval = FALSE}
library(autostsm)

##### Unemployment rate #####
#Not seasonally adjusted
data("UNRATENSA", package = "autostsm") #From FRED
UNRATENSA = data.table(UNRATENSA, keep.rownames = TRUE)
colnames(UNRATENSA) = c("date", "y")
UNRATENSA[, "date" := as.Date(date)]
UNRATENSA[, "y" := as.numeric(y)]
stsm = stsm_estimate(UNRATENSA, verbose = TRUE)
stsm_fc = stsm_forecast(stsm, y = UNRATENSA, n.ahead = floor(stsm$freq)*3, plot = TRUE)
stsm_fc = merge(stsm_fc, 
                stsm_detect_anomalies(stsm, y = UNRATENSA, plot = TRUE), 
                by = "date", all = TRUE)
stsm_fc = merge(stsm_fc, 
                stsm_detect_breaks(stsm, y = UNRATENSA, plot = TRUE, show_progress = TRUE), 
                by = "date", all = TRUE)
```

