% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/final_loadings_families.R
\name{find_neurons_family}
\alias{find_neurons_family}
\title{Predict the families for a list of neurons}
\usage{
find_neurons_family(listneurons, computedens = FALSE)
}
\arguments{
\item{listneurons}{Either names of flycircuit neurons or a
\code{\link[nat]{neuronlist}} object containing actual neurons (which will
be passed to \code{\link{voxel_dens}}).}

\item{computedens}{Whether or not to compute densities (which would be
necessary for non-flycircuit neurons)}
}
\description{
Predict the families for a list of neurons
}
\details{
This will return
}
\seealso{
Other find-family: \code{\link{find_scores_family}},
  \code{\link{neurons_against_fam}}
}

