% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/final_loadings_families.R
\name{list_scores_neurons_cv_fun}
\alias{list_scores_neurons_cv_fun}
\title{Find the list of scores of all the neurons for the cross validation approach}
\usage{
list_scores_neurons_cv_fun(listneurons = familynblast::fc_neuron_typec,
  computedens = FALSE, zeroscore = -100)
}
\arguments{
\item{listneurons}{Either names of flycircuit neurons or a
\code{\link[nat]{neuronlist}} object containing actual neurons (which will
be passed to \code{\link{voxel_dens}}).}

\item{computedens}{Whether or not to compute densities (which would be
necessary for non-flycircuit neurons)}

\item{zeroscore}{Log score to use when a neuron occupies a supervoxel that is
not overlapped by any family member.}
}
\description{
Find the list of scores of all the neurons for the cross validation approach
}

