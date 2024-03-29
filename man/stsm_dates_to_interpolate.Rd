% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stsm_dates_to_interpolate.R
\name{stsm_dates_to_interpolate}
\alias{stsm_dates_to_interpolate}
\title{Create dates to interpolate}
\usage{
stsm_dates_to_interpolate(y, dates, exo = NULL, interpolate)
}
\arguments{
\item{y}{Univariate time series of data values.}

\item{dates}{Vector of date values for y}

\item{exo}{Matrix of exogenous variables. Can be used to specify regression effects or other seasonal effects like holidays, etc.}

\item{interpolate}{Character string of how to interpolate}
}
\value{
List of the data, dates, and exo
}
\description{
Create dates to interpolate
}
\examples{
\dontrun{
#GDP Not seasonally adjusted
library(autostsm)
data("NA000334Q", package = "autostsm") #From FRED
NA000334Q = data.table(NA000334Q, keep.rownames = TRUE)
colnames(NA000334Q) = c("date", "y")
NA000334Q[, "date" := as.Date(date)]
NA000334Q[, "y" := as.numeric(y)]
NA000334Q = NA000334Q[date >= "1990-01-01", ]
dates_interp = stsm_dates_to_interpolate(y = NA000334Q$y, dates = NA000334Q$date, 
interpolate = "monthly")
}
}
