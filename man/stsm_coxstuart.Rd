% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stsm_stat_tests.R
\name{stsm_coxstuart}
\alias{stsm_coxstuart}
\title{Cox-Stuart Test}
\usage{
stsm_coxstuart(
  y,
  type = c("trend", "deviation", "dispersion"),
  sig_level = 0.01
)
}
\arguments{
\item{y}{input data}

\item{type}{Type of test: "trend", "deviation", or "dispersion"
If type = "trend", test for changes in trend
If type = "deviation", test for changes in deviation
If type = "dispersion", test for changes in dispersion (range)}

\item{sig_level}{Significance level to determine statistically significant seasonal frequencies}
}
\value{
list describing the results
}
\description{
Taken from the `tsutils` package. Performs the Cox-Stuart test for trend, deviation, or dispersion
}
