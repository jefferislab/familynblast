#' All flycircuit neurons with the correct labels we are working on
#' @name fc_neuron_typec
#' @docType data
NULL

#' All 137 families of flycircuit neurons
#' @details We still have to get rid of the ones that contains one neuron
#' @name correct_families
#' @docType data
NULL

#' The prior probabilities for all the families in correct_families
#' @name probability_correct_families
#' @seealso \code{\link{correct_families}}
#' @docType data
#' @examples
#' hist(probability_correct_families)
NULL

#' The numbers of neurons crossing each supervoxel, matrix of 7065
#' @name number_neurons_fromfam_insv
#' @docType data
NULL

#' The numbes of neurons crossing each supervoxel, matrix of 7065
#' @name number_neurons_fromfam_insv
#' @docType data
NULL

#' The probabilities of all the supervoxels in each family, matrix of nrow = 7065
#' @name probability_sv_knowing_family
#' @docType data
NULL

#' A table that contains the number of voxels in the supervoxels
#' @name svoxels.fcwb.table
#' @docType data
NULL

#' Supervoxel label field in FCWB space
#'
#' This was originally supplied by Kristin Branson in the JFRC2 template space
#' and then bridged to the FCWB space containing the flycircuit neurons. It
#' contains 7065 non-zero levels defining supervoxels covering nc82 positive
#' regions in the brain (i.e. neuropil).
#' @name svoxels.fcwb
#' @docType data
NULL

#' Precalculated super voxel density representation of flycircuit neurons
#'
#' This object is a sparse \code{\link[Matrix]{Matrix}} since it has many zeros.
#' It has 7066 columns (all supervoxels including catch-all 0 level) and 16129
#' rows for each flycircuit neuron. The rows are named with the flycircuit
#' \code{gene_name} - see \code{\link[flycircuit]{fc_gene_name}}. CHECK
#' @name voxel_dens_allneurons
#' @docType data
NULL

#' Selection of flycircuit Kenyon cells used for subfamily model
#'
#' @description \code{kcs.subfam.training} consists of 180 neurons in 6
#'   sub-families, while \code{kcs.subfam.test} are the remaining
#'   \code{flycircuit} KCs.
#'
#'
#' @name kcs.subfam.training
#' @aliases kcs.subfam.test
#' @docType data
NULL

#' Probability matrices for flycircuit Kenyon cell subfamily model
#'
#' @description \code{probability_sv_knowing_subfamily} is a 7065 x 6 matrix
#'   describing the probability of each supervoxel being crossed by one of the
#'   neurons in the \code{\link{kcs.subfam.training}} training set. Computed by
#'   \code{\link{prior_prob_svscores}}.
#'
#'   \code{probability_subfamily_kcs} contains the prior probability for each of
#'   the 6 subfamilies in the complete flycircuit dataset. Computed by
#'   \code{\link{prior_prob_subfam}}.
#'
#' @description \code{probabilities_nblastscores_kcs} contains the probabilities
#'   for NBLAST scores to fall into a certain bin (11 bins in range -1 to +1 )
#'   for neurons in the 6 different KC families. It is a 3D array with
#'   dimensions. Computed by \code{\link{prior_prob_nblastscores}}.
#'
#' @seealso \code{\link{prior_prob_svscores}}, \code{\link{prior_prob_subfam}},
#'   \code{\link{prior_prob_nblastscores}}
#' @name probability_sv_knowing_subfamily
#' @aliases probability_subfamily_kcs probabilities_nblastscores_kcs
#' @docType data
#' @examples
#' # dimensions of NBLAST probabilities by subfamily
#' dim(probabilities_nblastscores_kcs)
#' dimnames(probabilities_nblastscores_kcs)
NULL

