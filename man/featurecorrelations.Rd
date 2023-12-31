% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MicrogliaMorphologyR_functions.R
\name{featurecorrelations}
\alias{featurecorrelations}
\title{Correlation heatmap across morphology features}
\usage{
featurecorrelations(data, featurestart, featureend, rthresh, pthresh, title)
}
\arguments{
\item{data}{is your input data frame}

\item{featurestart}{is first column number of morphology measures}

\item{featureend}{is last column number of morphology measures}

\item{rthresh}{is cutoff threshold for significant correlation values}

\item{pthresh}{is cutoff threshold for significant p-values}

\item{title}{is what you want to name your heatmap}
}
\description{
'featurecorrelations' allows you to generate a heatmap depicting significant correlations across features
}
