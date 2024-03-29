% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stsm_detect_breaks.R
\name{stsm_detect_breaks}
\alias{stsm_detect_breaks}
\title{Detect Structural Breaks}
\usage{
stsm_detect_breaks(
  model,
  y,
  components = c("trend", "cycle", "seasonal"),
  freq = NULL,
  exo_obs = NULL,
  exo_state = NULL,
  sig_level = 0.01,
  ci = 0.8,
  smooth = TRUE,
  plot = FALSE,
  cores = NULL,
  show_progress = FALSE
)
}
\arguments{
\item{model}{Structural time series model estimated using stsm_estimate.}

\item{y}{Univariate time series of data values. May also be a 2 column data frame containing a date column.}

\item{components}{Vector of components to test for structural breaks}

\item{freq}{Frequency of the data (1 (yearly), 4 (quarterly), 12 (monthly), 365.25/7 (weekly), 365.25 (daily)), default is NULL and will be automatically detected}

\item{exo_obs}{Matrix of exogenous variables to be used in the observation equation.}

\item{exo_state}{Matrix of exogenous variables to be used in the state matrix.}

\item{sig_level}{Significance level to determine statistically significant anomalies}

\item{ci}{Confidence interval, value between 0 and 1 exclusive.}

\item{smooth}{Whether or not to use the Kalman smoother}

\item{plot}{Whether to plot everything}

\item{cores}{Number of cores to use for break detection}

\item{show_progress}{Whether to show progress bar}
}
\value{
data table (or list of data tables) containing the dates of detected anomalies from the filtered and/or smoothed series
}
\description{
Detect structural breaks using the estimated structural time series model
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
stsm = stsm_estimate(NA000334Q)
breaks = stsm_detect_breaks(model = stsm, y = NA000334Q, plot = TRUE, cores = 2)
}
}
