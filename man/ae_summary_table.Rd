% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tables.R
\name{ae_summary_table}
\alias{ae_summary_table}
\title{Posterior Summary Statistics for Safety Topics of Interest}
\usage{
ae_summary_table(
  input_data,
  cb_list_ctrl,
  cb_list_trt,
  saf_topic,
  explore = FALSE,
  seed = NULL
)
}
\arguments{
\item{input_data}{dataframe from data_table_prep() including the new trial}

\item{cb_list_ctrl}{expects List with control Arm indicators}

\item{cb_list_trt}{expects List with treatment Arm indicators}

\item{saf_topic}{the safety topic}

\item{explore}{for exploring prior data conflict, expected by posterior_dist()}

\item{seed}{a seed}
}
\value{
a table of posterior summary statistics for each event type for rates and proportions
}
\description{
Creates a table of posterior summary statistics for each
safety topic of interest. Creates the arrays fot thr pdf.
}
