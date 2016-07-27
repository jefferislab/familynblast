# Function find_family_and_subfamily --------------------------------------
#####  find_family_and_subfamily  is a function that takes for arguments a list of neurons listneurons, a subfamily
##### computens is set to FALSE by default but is required if the neurons are not in the dps neuronlist
##### This functionsfinds the family of a neuron and if the family contains subfamilies, it finds the subfamily of the neuron.
##### To find the family, the approach is the probabilistic approach, what is the probability to have a certain series of supervoxel if we are from the same family, and for the subfamily,
##### probability to have a series of supervoxels and series of nblast scores.
##### returns a matrix wit the family in first row, subfamily in second row and name the name of the neuron. index_subf should give the index of the family, and the index of
##### the subfamily within the family

## Just to test i create the list of the subfamilies :
#kcs.subfam.training # is a list of the subfamilies and in it, the neurons from the data set of the subfamily.
# subfamilies in our examples will just be names_kcsm
# The pre-computed data that could be add as arguments/ re-computed if necessary are probabilities_nblastscores_kcs and probabily_subfamily_kcs






#' compute_score_subfamily, computes all the scores of the list of neurons
#' against theneurons of the subfamilie
#'
#' @param listofneurons A list of neurons to test, list of neurons types named
#'   by flycircuit identifiers
#' @param zeronbl The value to penalize for zero in nblast
#' @param zerosv The value to penalize for zero in supervoxels
#' @param subfamilies A list of neurons from the training set, neurons type
#'   named by flycircuit identifiers
#'
#' @return returns a matrix with in rows the neurons and in columns the
#'   subfamilies. We have the scores of the neurons against the subfamilies
#' @export
#' @importFrom nat.nblast nblast
#' @examples
#' # example using actual neuron objects in a neuronlist distributed with nat
#' # package
#' compute_score_subfamily(kcs20)
#'
#' # example using named flycircuit neurons - must have dps object of flycircuit
#' # neurons loaded
#' compute_score_subfamily(familynblast::kcs.subfam.test[25:35])
compute_score_subfamily = function(listofneurons,zeronbl = -0.9, zerosv = -100,subfamilies = familynblast::kcs.subfam.training){
  scores_neurons_families_nblast = matrix(0,nrow = length(listofneurons), ncol = length(unique(familynblast::kcs.subfam.training)))
  rownames(scores_neurons_families_nblast) = names(listofneurons)
  colnames(scores_neurons_families_nblast) = unique(subfamilies)

  scorecut=seq(-1,1,0.2)

  for(i in seq_along(listofneurons)){
    print(paste("working on neuron",i))
    # FIXME Deal with situation where we have neurons that are not part of dps!
    neuron = dps[names(listofneurons)[i]]                                             ## Take the neuron
    #-------------------------------------------------------------------  Compute the nblast scores against the set of 180 neurons

    maxscore=max(nat.nblast::smat.fcwb)*nrow(neuron[[1]]$points)
    scores_neuron = c()   ## List of the nblast scores of the neuron against the 180 neurons of reference
    for(j in 1:length(subfamilies)){

      scores_neuron = c(scores_neuron,nblast(neuron,dps[names(subfamilies)[j]])/maxscore)
    }
    names(scores_neuron) = names(subfamilies)
    #-------------------------------------------------------------------  Find the list of supervoxels crossed

    sv_neuron = names(voxel_dens_allneurons[names(listofneurons)[i],voxel_dens_allneurons[names(listofneurons)[i],]>0])
    # drop catch-all 0 voxel
    sv_neuron=setdiff(sv_neuron, "0")

    #-------------------------------------------------------------------  Compute the score
    for(l in seq_along(unique(subfamilies))){              ## Now compute the score on each family
      names3 = subfamilies[which(subfamilies ==  unique(subfamilies)[l])]
      for(k in names(names3)){          ## Part on nblast
        # FIXME this needs to be made more generic so it can work with neurons
        # other than KCs
        if (familynblast::probabilities_nblastscores_kcs[k,l,findInterval(scores_neuron[k],scorecut)]==0){
          scores_neurons_families_nblast[names(listofneurons)[i],l] =  scores_neurons_families_nblast[names(listofneurons)[i],l] +zeronbl
        }else{
          # message("k=",k," l=",l," fi=", findInterval(scores_neuron[k],scorecut))
          logp=log(familynblast::probabilities_nblastscores_kcs[k,l,findInterval(scores_neuron[k],scorecut)])
          scores_neurons_families_nblast[names(listofneurons)[i],l] =  scores_neurons_families_nblast[names(listofneurons)[i],l] +
            +logp
        }
      }
      scores_neurons_families_nblast[names(listofneurons)[i],l] = scores_neurons_families_nblast[names(listofneurons)[i],l]/length(subfamilies)

      ## Part on supervoxels
      prob = familynblast::probability_sv_knowing_subfamily[sv_neuron,l]
      zeros = sum(prob==0)
      prob2 = prob[which(prob!=0)]
      score_f = sum(log(prob2))
      score_f = score_f+zerosv*zeros
      scores_neurons_families_nblast[names(listofneurons)[i],l] =  scores_neurons_families_nblast[names(listofneurons)[i],l] + score_f/length(sv_neuron) +
        log(familynblast::probability_subfamily_kcs[l])
    }
  }
  return(scores_neurons_families_nblast)
}

#' Compute the prior probabilities of the nblast scores
#'
#' @param subfamilies list of neurons from the training set, neurons type named by flycircuit identifiers
#' @param bins bins to cut the nblast scores
#' @return
#' @export
#'
#' @examples
#' @importFrom nat.nblast nblast
prior_prob_nblastscores = function(subfamilies = familynblast::kcs.subfam.training, bins=seq(-1,1,0.2)){
  probabilities_nblastscores = array(0,dim=c(length(subfamilies),length(unique(subfamilies)),length(bins)))
  dimnames(probabilities_nblastscores) = list(names(subfamilies),unique(subfamilies),bins)
  for(i in seq_along(subfamilies)){              ## For every neuron of the 1080 that will make the set
    print(i)
    neuron = dps[names(subfamilies)[i]]
    names2 = setdiff(names(subfamilies),names(neuron))
    ## Get rid of the query neuron from the rest of the set of 180
    for (j in seq_along(unique(subfamilies))){     ## For every family
      names3 = subfamilies[names2]
      names3 = names3[which(names3 ==  unique(subfamilies)[j])]
      for(l in seq_along(names3)){
        scoren = nblast(neuron,dps[names(names3[l])])/(max(nat.nblast::smat.fcwb)*nrow(neuron[[1]]$points))   ## Score between the neuron and the neuron l of the family j, normalized
        probabilities_nblastscores[i,j,findInterval(scoren,bins)] = probabilities_nblastscores[i,j,findInterval(scoren,bins)]+1
      }
    }
  }
  probabilities_nblastscores = probabilities_nblastscores/length(subfamilies)
  return(probabilities_nblastscores)
}



#' Compute the prior probabilities of the subfamilies
#'
#' @param subfamiliesall list of neurons neurons type named by flycircuit identifiers, all neurons from the subfamilies
#'
#' @return
#' @export
#'
#' @examples
## subfamiliesall should be the list of all the neurons in
prior_prob_subfam = function(subfamiliesall = familynblast::kcs.subfam.test){
  probabily_subfamily = c()
  for(i in unique(subfamiliesall)){
    print(i)
    probabily_subfamily = c(probabily_subfamily, table(subfamiliesall)[i]/sum(table(subfamiliesall)))
  }
  names(probabily_subfamily) = unique(subfamiliesall)
  return(probabily_subfamily)
}

#' Compute the prior probabilities of the supervoxels
#'
#' @param subfamilies list of neurons from the training set, neurons type named by flycircuit identifiers
#' @param computedens if neurons are not from flycircuit dataset
#'
#' @return
#' @export
#'
#' @examples
prior_prob_svscores = function(subfamilies = familynblast::kcs.subfam.training,computedens = FALSE){
  if (computedens == TRUE){
    voxel_dens_allneurons = voxel_dens(names(subfamilies))
    rownames(voxel_dens_allneurons) = names(subfamilies)
  } else {
    voxel_dens_allneurons = familynblast::voxel_dens_allneurons
  }
  number_neurons_fromsubfam_insv = matrix(0,nrow=7065,ncol=length(unique(subfamilies)))
  for (i in seq_along(unique(subfamilies))){
    print(paste("moving to subfamily",i))
    setofneurons = subfamilies[which(subfamilies ==unique(subfamilies)[i])]
    for(j in 1:7065){
      number_neurons_fromsubfam_insv[j,i] = sum(voxel_dens_allneurons[names(setofneurons),j]!=0)
    }
  }
  # Filling up the matrix of probabilities for supervoxels knowing the family --------
  ### We have to determine the probability of having sj knowing we are in the family i
  probability_sv_knowing_subfamily = matrix(0,nrow=7065,ncol=length(unique(subfamilies)))
  for(l in seq_along(unique(subfamilies))){
    for (k in 1:7065){
      probability_sv_knowing_subfamily[k,l] = number_neurons_fromsubfam_insv[k,l]/length(subfamilies)   ##### BE CAREFUL, WE DIVIDE HERE BY length(subfamilies)  BUT IT
      #### WILL HAVE TO BE THE NUMBER OF NEURONS IN THE FAMILY LINKED TO THE SUBFAMILY
    }
  }
  return(probability_sv_knowing_subfamily)

}
