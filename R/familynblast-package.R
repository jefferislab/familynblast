#' NBLAST to identify neuron families and subfamilies
#'
#' This package implements two strategies for predicting the class of neurons
#' based on a training set of neurons. The first strategy using a representation
#' based on dividing the brain into a set of 7065 supervoxels. This strategy is
#' targetd at broad \bold{families} of neurons (think Kenyon cells or the three
#' main classes of Kenyon cell).
#'
#' The second strategy combines the first strategy with the additional use of
#' NBLAST scores for the similarity between test neurons and a training set that
#' has been divided into subfamilies (think a/B core vs surface Kenyon cells).
#'
#' A collection of pre-computed families is distributed based on annotation work
#' carried out on the FlyCircuit single neuron collection.
#'
#' @section Family Models: You can evaluate neurons against an existing family
#'   model using \code{\link{neurons_against_fam}}. You can compute scores
#'   against all families using \code{\link{find_scores_family}} or predict the
#'   family (including the prediction that there is no matching family) for a
#'   list of neurons using \code{\link{find_neurons_family}}.
#'
#'   You can generate new family models using the
#'   \code{\link{create_probab_sv_knowing_fam}} function.
#'
#'   There is pre-computed data distributed with the package for 137 existing
#'   models based on FlyCircuit neurons registered against the FCWB template
#'   brain. See \code{\link{probability_sv_knowing_family}}.
#'
#' @section Subfamily Models: FIXME add details
#'
#' @name familynblast-package
#' @aliases familynblast
#' @docType package
#' @keywords package
#' @encoding UTF-8
#' @references Developed as part of a research internship by Mélina Durande
#'   (ENS Lyon) in the lab of Gregory Jefferis (Neurobiology Division, MRC LMB,
#'   Cambridge) in May-July 2016. (FIXME Link to report?).
#'
#'   The supervoxel parcellation is unpublished work due to Kristin Branson
#'   (HHMI Janelia research Campus).
#'
#'   The registered data and annotations for the FlyCircuit neurons were
#'   described in:
#'
#'   Costa, M., Ostrovsky, A.D., Manton, J.D., Prohaska, S., and Jefferis,
#'   G.S.X.E. (2014). NBLAST: Rapid, sensitive comparison of neuronal structure
#'   and construction of neuron family databases. Biorxiv preprint.
#'   \href{http://dx.doi.org/10.1101/006346}{doi: 10.1101/006346}. The original
#'   FlyCircuit data were obtained from:
#'
#'   Chiang A.S., Lin C.Y., Chuang C.C., Chang H.M., Hsieh C.H., Yeh C.W., Shih
#'   C.T., Wu J.J., Wang G.T., Chen Y.C., Wu C.C., Chen G.Y., Ching Y.T., Lee
#'   P.C., Lin C.Y., Lin H.H., Wu C.C., Hsu H.W., Huang Y.A., Chen J.Y., et al.
#'   (2011). Three-dimensional reconstruction of brain-wide wiring networks in
#'   Drosophila at single-cell resolution. Curr Biol 21 (1), 1--11.
NULL
