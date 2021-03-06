% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stsm_detect_multiplicative.R
\name{stsm_detect_multiplicative}
\alias{stsm_detect_multiplicative}
\title{Detect if log transformation is best}
\usage{
stsm_detect_multiplicative(y, freq, sig_level = 0.01, prior = NULL)
}
\arguments{
\item{y}{an object created from stsm_detect_frequency}

\item{freq}{Frequency of the data}

\item{sig_level}{Significance level to determine statistically significant seasonal frequencies}

\item{prior}{A data table created by stsm_prior}
}
\value{
a logical indicating if the model should be multiplicative or not
}
\description{
Detect if log transformation is best
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
multiplicative = stsm_detect_multiplicative(y = NA000334Q$y, freq = 4)
}
}
