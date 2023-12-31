% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MicrogliaMorphologyR_functions.R
\name{stats_morphologymeasures.animal}
\alias{stats_morphologymeasures.animal}
\title{Stats analysis: individual morphology measures}
\usage{
stats_morphologymeasures.animal(data, model, posthoc1, posthoc2, adjust)
}
\arguments{
\item{data}{is your input data frame}

\item{model}{is your linear mixed model (e.g., Value ~ Treatment*Sex + (1|MouseID))}

\item{posthoc1}{is your posthoc comparisons (e.g., when considering sex: ~Treatment|Sex)}

\item{posthoc2}{is your posthoc comparisons (e.g., when not considering sex: ~Treatment)}

\item{adjust}{is your method of multiple test correction (from \code{emmeans} package: "tukey","scheffe","sidak","dunnettx","mvt","holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr","none"). See "P-value adjustments" section under ?emmeans::summary.emmGrid for more information.}
}
\description{
Linear mixed model to statistically assess how your experimental variables of interest
influence each morphology measure, at the animal level
}
\details{
The stats_morphologymeasures.animal function fits a linear model using the \code{lm} function
for each morphology measure individually within your dataset. Posthocs are run for each morphology measure
individually and bound together into the final dataframe that is output by this function.
}
