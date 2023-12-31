% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MicrogliaMorphologyR_functions.R
\name{metadata_columns}
\alias{metadata_columns}
\title{Tidy up metadata}
\usage{
metadata_columns(data, metadata, sep)
}
\arguments{
\item{data}{is your tidied dataframe containing all of your numerical morphology data}

\item{metadata}{is your metadata listed out in order that they appear in 'Name' column: e.g., c("Cohort","Sex","MouseID")}

\item{sep}{is the character, "_" or "-", that separates your metadata}
}
\description{
'metadata_columns' allows you to extract individual metadata columns from the UniqueID column.
Each piece of metadata should be separated by a common deliminator, 'sep' parameter
}
