#' Computes the number of points from a setofneurons in each supervoxel
#'
#' @param setofneurons \code{\link[nat]{neuronlist}} object containing set
#' @param svoxels A \code{\link[nat]{im3d}} object containing defining a set of
#'   supervoxels. Defaults to \code{\link{svoxels.fcwb}}.
#'
#' @param svoxels.table A \code{\link{table}} indicating how many pixels are
#'   present in each supervoxel. Defaults to \code{\link{svoxels.fcwb.table}}.
#'
#' @return A \code{matrix} with length(setofneurons) rows and
#'   \code{length(svoxels.table)} columns.
#' @export
#' @importFrom nat xyzmatrix coord2ind as.im3d
#' @examples
#' # convert a sample set of neurons to supervoxel representation
#' kcdens=voxel_dens(nat::kcs20)
voxel_dens = function(setofneurons,
                      svoxels = familynblast::svoxels.fcwb,
                      svoxels.table = familynblast::svoxels.fcwb.table) {
  svoxel.output.neu = matrix(0,
                             nrow = length(setofneurons),
                             ncol = length(svoxels.table))
  colnames(svoxel.output.neu) = names(svoxels.table)
  rownames(svoxel.output.neu) = names(setofneurons)
  for (n in names(setofneurons)) {
    message("working on neuron:", n)
    # find 1d indices of 3D coords into FCWB space
    idxs = coord2ind(xyzmatrix(setofneurons[[n]]), imdims = svoxels)
    # then get label values for those coords
    svoxel.ids = svoxels[idxs]
    # compute density i.e. number of points in each domain
    tsvoxel.ids = table(svoxel.ids)
    svoxel.output.neu[n, names(tsvoxel.ids)] = tsvoxel.ids
  }
  return(svoxel.output.neu)
}


#' Computing the score of one neuron against a particular family, given the
#' index of the name of the family in correct_families
#'
#' @param listneurons Either names of flycircuit neurons or a
#'   \code{\link[nat]{neuronlist}} object containing actual neurons (which will
#'   be passed to \code{\link{voxel_dens}}).
#' @param familyind Index of the family we want to test neurons against (either
#'   numeric or character vector, see examples)
#' @param computedens Whether or not to compute densities (which would be
#'   necessary for non-flycircuit neurons)
#' @param zeroscore Log score to use when a neuron occupies a supervoxel that is
#'   not overlapped by any family member.
#'
#' @return
#' @export
#' @importFrom nat is.neuronlist
#' @family find-family
#'
#' @examples
#' # One neuron against a specific family
#' neurons_against_fam(kcs20[1], "gamma Kenyon cell")
#' # multiple neurons against a specific family
#' neurons_against_fam(kcs20, "gamma Kenyon cell")
#'
#' # find the KC families that we know about
#' kcfams=grep("Kenyon", names(correct_families), value = T)
#' # compare all neurons against all families
#' res=sapply(kcfams, function(fam) neurons_against_fam(kcs20, fam))
#' # make a prediction using the class of max score
#' pred=colnames(res)[apply(res,1, which.max)]
#' maxscore=apply(res, 1,max)
#' # compare with manually defined type, NB there is no family definition
#' # for alpha'/beta' KCs at the moment so these are predicted as a/B
#' # but have slightly lower scores than the real a/B neurons
#' cbind(kcs20[,c("Name","type")], pred, max=maxscore)
neurons_against_fam  = function(listneurons,familyind, computedens = FALSE, zeroscore = -100){
  if (is.neuronlist(listneurons)){
    voxel_dens_allneurons = voxel_dens(listneurons)
    # drop the 0 level (i.e. not in a well-defined supervoxel)
    if("0" %in% colnames(voxel_dens_allneurons))
      voxel_dens_allneurons=voxel_dens_allneurons[,-1, drop=FALSE]
    listneurons=names(listneurons)
    rownames(voxel_dens_allneurons) = listneurons
  } else {
    voxel_dens_allneurons=familynblast::voxel_dens_allneurons
    # if the character vector of input neurons is named then we assume that we
    # want those names not the input elements
    if(!is.null(names(listneurons))) listneurons=names(listneurons)
  }
  Score =c()
  for (k in seq_along(listneurons)){                                        ## for each neuron
    print(k)
    neuron = listneurons[k]
    sv_neurons = names(voxel_dens_allneurons[neuron,voxel_dens_allneurons[neuron,]>0])
    prob = familynblast::probability_sv_knowing_family[sv_neurons,familyind]
    zeros = sum(prob==0)
    prob2 = prob[which(prob!=0)]
    score_f = sum(log(prob2))
    score_f = score_f+zeroscore*zeros
    Score = c(Score,score_f)
    }
    Score = Score/length(sv_neurons)+log(familynblast::probability_correct_families[familyind])
  return(Score)
}


#' Find scores of list of neurons against all families
#'
#' @inheritParams neurons_against_fam
#' @return
#' @export
#' @family find-family
#'
#' @examples
find_scores_family = function(listneurons, computedens = FALSE, zeroscore = -100){
  matrix_scores_neurons = matrix(0,nrow =listneurons, ncol = length(familynblast::correct_families))
  for(k in 1:length(familynblast::correct_families)){
    matrix_scores_neurons[,k] = neurons_against_fam(listneurons,familyind = k,computedens=FALSE, zeroscore=-100)
  }

  list_scores_neurons = list()
  for(l in 1:length(listneurons)){
    list_scores_neurons = append(list_scores_neurons,matrix_scores_neurons[l])
  }
  names(list_scores_neurons) = names(listneurons)
  for(i in seq_along(list_scores_neurons)){
    names(list_scores_neurons[[i]])=names(familynblast::correct_families)
  }
  return(list_scores_neurons)
}


#' Predict the families for a list of neurons
#'
#' @details This will return
#' @inheritParams neurons_against_fam
#'
#' @return
#' @export
#'
#' @family find-family
find_neurons_family = function(listneurons,computedens = FALSE){
  familiesof_listofneurons = c()
  list_scores_neurons = find_scores_family(listneurons, computedens=FALSE)
  for(i in seq_along(listneurons)){
    if(sum(is.na(list_scores_neurons[[i]])) != 137 & rev(sort(list_scores_neurons[[i]]))[1]>-12){
        familiesof_listofneurons = c(familiesof_listofneurons,names(which.max(list_scores_neurons[[i]])))
      }else{
        familiesof_listofneurons = c(familiesof_listofneurons,"no family for this neuron")
      }

  }
  names(familiesof_listofneurons) = names(listneurons)
  return(familiesof_listofneurons)
}


#' Find the percentage of correct hits within the nth highest scores
#'
#' @param listneurons listneurons is the list of neurons names we want to find the family
#' @param nb The number of high scores to consider
#'
#' @return The percentage of correct labelling when we keep the first top score,
#' the top two scores, top nb scores.
#' @export
find_percentage_correct_hits = function(listneurons, nb = 3){
  Percents = c()
  scorestouse = list_scores_neurons_cv_fun(listneurons)
  for(j in 1:nb){
    correct_hits=0
    for(l in seq_along(listneurons)){
      name_hit=names(scorestouse[[l]])[rev(order(scorestouse[[l]]))[1:j]]
      if(listneurons[l] %in% name_hit){
        print("well done !!")
        correct_hits = correct_hits+1
      }else{
        print("Oups, you're wrong")
      }
      percentage_correct = correct_hits/length(scorestouse)
    }
    Percents = c(Percents, percentage_correct)
  }
  return(Percents)
}


#' Find the list of scores of all the neurons for the cross validation approach
#'
#' @inheritParams neurons_against_fam
#'
#' @return
#' @export
list_scores_neurons_cv_fun = function(listneurons = familynblast::fc_neuron_typec, computedens = FALSE,zeroscore = -100 ){
  if (computedens == TRUE){
    voxel_dens_allneurons = voxel_dens(listneurons)
    rownames(voxel_dens_allneurons) = listneurons
  }
  list_scores_neurons_cv = list()
  for (k in seq_along(listneurons)){                                        ## for each neuron
    print(k)
    neuron = listneurons[k]
    ### Find the neuron and take it out
    fc_neuron_typecb = listneurons[which(names(listneurons)!=names(neuron))]
    correct_familiesb = familynblast::correct_families
    correct_familiesb[[neuron[[1]]]] = familynblast::correct_families[[neuron[[1]]]][names(familynblast::correct_families[[neuron[[1]]]])!=names(neuron)]

    ### Process with the families
    probability_correct_families_b = familynblast::probability_correct_families
    indx = which(names(correct_familiesb)==neuron[[1]])
    probability_correct_families_b[,indx]=length(correct_familiesb[[indx]])/length(fc_neuron_typecb)

    ### Process to compute the number of neuron in the family
    number_neurons_fromfam_insv_b = familynblast::number_neurons_fromfam_insv

    if(length(ncol(voxel_dens_allneurons[names(correct_familiesb[[indx]]),]!=0))==0){
      number_neurons_fromfam_insv_b[,indx] = sum(voxel_dens_allneurons[names(correct_familiesb[[indx]]),]!=0)
    }else{
      number_neurons_fromfam_insv_b[,indx] = colSums(voxel_dens_allneurons[names(correct_familiesb[[indx]]),]!=0)
    }

    probability_sv_knowing_familyb = familynblast::probability_sv_knowing_family
    probability_sv_knowing_familyb[,indx] = number_neurons_fromfam_insv_b[,indx]/length(familynblast::correct_families[[indx]])


    sv_neurons_b = names(voxel_dens_allneurons[names(neuron),voxel_dens_allneurons[names(neuron),]>0])

    Score =c()
    for(l in 1:length(correct_familiesb)){
      prob = probability_sv_knowing_familyb[sv_neurons_b,l]
      zeros = sum(prob==0)
      prob2 = prob[which(prob!=0)]
      score_f = sum(log(prob2))
      score_f = score_f -100*zeros
      Score = c(Score,score_f)
    }
    Score = Score/length(sv_neurons_b)
    for(l in 1:length(correct_familiesb)){
      Score[l] = Score[l]+log(probability_correct_families_b[l])
    }
    Score= list(Score)
    list_scores_neurons_cv = append(list_scores_neurons_cv,Score)
  }
  names(list_scores_neurons_cv) = names(listneurons)

  for(i in seq_along(list_scores_neurons_cv)){
    names(list_scores_neurons_cv[[i]])=names(familynblast::correct_families)
  }
  return(list_scores_neurons_cv)
}

#' Calculate prior probability for a set of families
#'
#' The set of families should be of the form: list of families and in each list
#' the list of neurons in this family. An example is
#' \code{\link{correct_families}}.
#'
#' @param setoffamilies is a list of neurons from the different families.
#'
#' @return
#' @export
#' @seealso \code{\link{correct_families}}
create_probab_families = function(setoffamilies){
  probability_correct_families = matrix(0,nrow=1,ncol=length(setoffamilies))
  for (i in seq_along(setoffamilies)){
    allneu = sum(sapply(1:length(setoffamilies), function(n) length(setoffamilies[[n]])))
    probability_correct_families[1,i]= length(setoffamilies[[i]])/allneu
  }
return(probability_correct_families)
}


#' Create probability matrix with supervoxels as rows and families as columns
#' for a setoffamilies
#'
#' @param setoffamilies A named list of families, with each list item containing
#'   a character vector of ids of neurons within that family. The ids are
#'   assumed to be FlyCircuit ids if the \code{neurons} argument is missing
#' @param db A \code{\link[nat]{neuronlist}} containing neuron structures. If
#'   present it is assumed that you want \code{\link{voxel_dens}} to compute the
#'   supervoxel densities of the neurons named in \code{setoffamilies}.
#' @param ... Additional arguments passed to \code{\link{voxel_dens}}
#' @return A matrix with supervoxels as rows and families as columns (and
#'   appropriate row and column names). The 0 level will always be dropped.
#' @details if \code{neurons} argument is missing then it is assumed that
#'   neurons are FlyCircuit neurons present in the pre-computed
#'   \code{\link{voxel_dens_allneurons}}; there will be an error if this is not
#'   the case.
#' @export
#' @importFrom flycircuit fc_gene_name
#' @seealso \code{\link{voxel_dens}},
#' @examples
#' kcs20.setoffamilies=by(names(kcs20), kcs20[,'type'], as.character)
#' kcs20.prob=create_probab_sv_knowing_fam(kcs20.setoffamilies, db=kcs20)
create_probab_sv_knowing_fam = function(setoffamilies, db=NULL, ...){
  neuron_ids=unlist(setoffamilies, use.names=F)
  if (is.null(db)){
    # convert any kind of FlyCircuit id to the preferred gene_name form
    neuron_ids=fc_gene_name(neuron_ids)
    # check ids are present in pre-computed data
    voxel_dens_allneurons=familynblast::voxel_dens_allneurons
    if(!all(neuron_ids%in%rownames(voxel_dens_allneurons)))
      stop("Some neuron ids in setoffamilies are not present in precomputed voxel_dens_allneurons.\n",
           "Check ids or supply neurons argument to allow computation from scratch!")
  } else {
    message("Computing voxel densities for ", length(neuron_ids), " neurons")
    voxel_dens_allneurons = voxel_dens(db[neuron_ids], ...)
  }
  # number of super voxels (excluding 0)
  nsvoxels=ncol(voxel_dens_allneurons)-1
  probability_sv_knowing_family = matrix(0,nrow=nsvoxels,ncol=length(setoffamilies))
  for (i in seq_along(setoffamilies)){
    print(paste("moving to family",i))
    setofneurons = setoffamilies[[i]]   ### ??
    neurons_per_sv=colSums(voxel_dens_allneurons[setofneurons,-1, drop=FALSE]!=0)
    probability_sv_knowing_family[,i] = neurons_per_sv/length(setofneurons)
  }
  colnames(probability_sv_knowing_family)=names(setoffamilies)
  rownames(probability_sv_knowing_family)=colnames(voxel_dens_allneurons)[-1]
  return(probability_sv_knowing_family)
}
