% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/s3_generics.R
\name{clusterProbability}
\alias{clusterProbability}
\title{Get the probability of the cluster assignments generated during a celda run.}
\usage{
clusterProbability(celda.mod, counts, log = FALSE)
}
\arguments{
\item{celda.mod}{A celda model object}

\item{counts}{The count matrix used to generate the model}

\item{log}{If FALSE, then the normalized conditional probabilities will be returned. If TRUE, then the unnormalized log probabilities will be returned.}
}
\description{
Get the probability of the cluster assignments generated during a celda run.
}
