% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ae_summary_table_helpers.R
\name{data_available}
\alias{data_available}
\title{Data availability}
\usage{
data_available(hist_trt, hist_ctr, trt_current, ctr_current)
}
\arguments{
\item{hist_trt}{Historical treatment data}

\item{hist_ctr}{Historical control data}

\item{trt_current}{current treatment data}

\item{ctr_current}{current control data}
}
\value{
A vector to inform where information is available
}
\description{
Gives as output 4 T/F values whether input data is available, ix stands for index
}
