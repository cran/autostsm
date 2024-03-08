#autostsm 3.1.4

## Minor changes

* bug fix to structural break detection with exogenous data

#autostsm 3.1.3

## Minor changes

* bug fix to weekly seasonal detection

#autostsm 3.1.2

## Minor changes

* using sentinel "_PACKAGE"

#autostsm 3.1.1

## Major changes

updated to fix issue calculating unobserved series at time 0 when there are lags

## Minor changes

* bug fixes involving exogenous inputs

#autostsm 3.1

## Minor changes

* bug fixes involving exogenous inputs

#autostsm 3.0.4

## Minor changes

* updated required version of kalmanfilter package

#autostsm 3.0.3

## Minor changes

* updated to use the kalmanfilter package rather than have its own internal Kalman filter routine, 
the kalmanfilter package is actually the Kalman filter routine that was built in this package

#autostsm 3.0.2

## Minor Changes

* updated Rcpp code to handle deprecation of << assignment operator

# autostsm 3.0.1

## Bug fixes

* setting deterministic models wrongly set other parameters to fixed
* build issues forced package to no longer depend on data.table but imports it

# autostsm 3.0

## Major changes

* Added ability to use exogenous data in the observation and/or state equations with options to specify what variables enter the equation for each unobserved components using R's formula syntax

## Bug fixes

* check for stationary cycle didn't work as expected
* an error occurred when no constraints were being used
* an error occurred during seasonal detection for a specific case
* setting cycle to trig didn't work as expected
* fixed some typos in the vignette
* updated vignette for exogenous observation and state data

# autostsm 2.1

## Major changes

* added structural break detection for the drift and trend selection
* added more user specified settings including setting cycle with NULL, FALSE, trig, arma with user specified arma order

## Bug fixes

* arma and cycle selection didn't work as expected

# autostsm 2.0

## Major changes

* Added structural interpolation models to interpolate lower frequency time series to a higher frequency (i.e quarterly to monthly)
