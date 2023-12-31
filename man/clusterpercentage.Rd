% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MicrogliaMorphologyR_functions.R
\name{clusterpercentage}
\alias{clusterpercentage}
\title{What are your morphology cluster percentages across variables of interest?}
\usage{
clusterpercentage(data, clustercol, ...)
}
\arguments{
\item{data}{is your input data frame}

\item{clustercol}{is the name of your column which contains cluster IDs. Make sure to put this in quotes.}

\item{...}{list out your variables of interest without combining (e.g., Cohort, Sex, Treatment)}
}
\description{
'clusterpercentage' groups your data by variables of interest then calculates percentage of each morphology cluster within those groups.
}
