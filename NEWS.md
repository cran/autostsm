# autostsm 3.0

## Major changes

* Added ability to use exogenous date in the observation and/or state equations with options to specify what variables enter the equation for each unobserved components using R's formula syntax

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