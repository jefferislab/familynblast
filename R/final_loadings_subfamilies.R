load(file="/Volumes/JData5/JPeople/Melina/Branson/data/probabilities_nblastscores_kcs") ### Contains the array that has, the families in columns,
#the neurons in lines and in depth the value of the score
load(file="/Volumes/JData5/JPeople/Melina/Branson/data/probability_sv_knowing_subfamily")
rownames(probability_sv_knowing_subfamily) = colnames(voxel_dens_allneurons)
load("/Volumes/JData5/JPeople/Melina/Branson/data/probabily_subfamily_kcs")
# test_set_kenyoncells ----------------------------------------------------
kcs0m=fc_neuron_type(regex="Kenyon")
kc1m = kcs0m[which(kcs0m=="gamma Kenyon cell")]
set.seed(32)
kc1m=sample(kc1m,30)
kc2m = kcs0m[which(kcs0m=="alpha'/beta' Kenyon cell")]
set.seed(32)
kc2m=sample(kc2m,30)
kc3m = kcs0m[which(kcs0m=="alpha/beta posterior Kenyon cell")]
set.seed(32)
kc3m=sample(kc3m,30)
kc4m = kcs0m[which(kcs0m=="alpha/beta core Kenyon cell")]
set.seed(32)
kc4m=sample(kc4m,30)
kc5m = kcs0m[which(kcs0m=="alpha/beta surface Kenyon cell")]
set.seed(32)
kc5m=sample(kc5m,30)
kc6m = kcs0m[which(kcs0m=="gamma dorsal Kenyon cell")]
set.seed(32)
kc6m=sample(kc6m,30)
kcslistm = list(kc1m,kc2m,kc3m,kc4m,kc5m,kc6m)
kcsm = unlist(list(kc1m,kc2m,kc3m,kc4m,kc5m,kc6m))
names_kcsm = list(names(kc1m),names(kc2m),names(kc3m),names(kc4m),names(kc5m),names(kc6m))
names(names_kcsm) = c("gamma Kenyon cell","alpha'/beta' Kenyon cell","alpha/beta posterior Kenyon cell","alpha/beta core Kenyon cell","alpha/beta surface Kenyon cell",
                      "gamma dorsal Kenyon cell")

# Rest of the data --------------------------------------------------------

kcs0=fc_neuron_type(regex="Kenyon")
kc1 = kcs0[which(kcs0=="gamma Kenyon cell")]
#set.seed(32)
#kc1=sample(kc1,30)
kc2 = kcs0[which(kcs0=="alpha'/beta' Kenyon cell")]
#set.seed(32)
#kc2=sample(kc2,30)

kc3 = kcs0[which(kcs0=="alpha/beta posterior Kenyon cell")]
#set.seed(32)
#kc3=sample(kc3,30)
kc4 = kcs0[which(kcs0=="alpha/beta core Kenyon cell")]
#set.seed(32)
#kc4=sample(kc4,30)
kc5 = kcs0[which(kcs0=="alpha/beta surface Kenyon cell")]
#set.seed(32)
#kc5=sample(kc5,30)
kc6 = kcs0[which(kcs0=="gamma dorsal Kenyon cell")]
#set.seed(32)
#kc6=sample(kc6,30)
kcslist = list(kc1,kc2,kc3,kc4,kc5,kc6)
kcs = unlist(list(kc1,kc2,kc3,kc4,kc5,kc6))
names_kcs = list(names(kc1),names(kc2),names(kc3),names(kc4),names(kc5),names(kc6))
names(names_kcs) = c("gamma Kenyon cell","alpha'/beta' Kenyon cell","alpha/beta posterior Kenyon cell","alpha/beta core Kenyon cell","alpha/beta surface Kenyon cell",
                     "gamma dorsal Kenyon cell")

### dotprops list
kcs.dps = dps[names(kcs)]

Lab = c("G","ApBp","ABp","ABc","ABs","Gd")
names(Lab) = names(names_kcsm)

# Function find_family_and_subfamily --------------------------------------
#####  find_family_and_subfamily  is a function that takes for arguments a list of neurons listneurons, a subfamily
##### computens is set to FALSE by default but is required if the neurons are not in the dps neuronlist
##### This functionsfinds the family of a neuron and if the family contains subfamilies, it finds the subfamily of the neuron.
##### To find the family, the approach is the probabilistic approach, what is the probability to have a certain series of supervoxel if we are from the same family, and for the subfamily,
##### probability to have a series of supervoxels and series of nblast scores.
##### returns a matrix wit the family in first row, subfamily in second row and name the name of the neuron. index_subf should give the index of the family, and the index of
##### the subfamily within the family

## Just to test i create the list of the subfamilies :
kcsm # is a list of the subfamilies and in it, the neurons from the data set of the subfamily.
# subfamilies in our examples will just be names_kcsm
# The pre-computed data that could be add as arguments/ re-computed if necessary are probabilities_nblastscores_kcs and probabily_subfamily_kcs
compute_score_subfamily = function(listofneurons,zeronbl = -9, zerosv = -100,subfamilies = kcsm){
  scores_neurons_families_nblast = matrix(0,nrow = length(listofneurons), ncol = length(unique(kcsm)))
  rownames(scores_neurons_families_nblast) = names(listofneurons)
  colnames(scores_neurons_families_nblast) = unique(subfamilies)


  # FIXME it should have been this
  scorecut=seq(-1,1,0.2)
  # but probabilities_nblastscores_kcs was computed with this in error
  scorecut=seq(-9,1,1)
  probabilities_nblastscores_kcs = probabilities_nblastscores_kcsbad
  #
  for(i in seq_along(listofneurons)){
    print(paste("working on neuron",i))
    neuron = dps[names(listofneurons)[i]]                                             ## Take the neuron
    #-------------------------------------------------------------------  Compute the nblast scores against the set of 180 neurons

    maxscore=max(smat.fcwb)*nrow(neuron[[1]]$points)
    scores_neuron = c()   ## List of the nblast scores of the neuron against the 180 neurons of reference
    for(j in 1:length(subfamilies)){

      scores_neuron = c(scores_neuron,nblast(neuron,dps[names(subfamilies)[j]])/maxscore)
    }
    names(scores_neuron) = names(subfamilies)
    #-------------------------------------------------------------------  Find the list of supervoxels crossed

    sv_neuron = names(voxel_dens_allneurons[names(listofneurons)[i],voxel_dens_allneurons[names(listofneurons)[i],]>0])

    #-------------------------------------------------------------------  Compute the score
    for(l in seq_along(unique(subfamilies))){              ## Now compute the score on each family
      names3 = subfamilies[which(subfamilies ==  unique(subfamilies)[l])]
      for(k in names(names3)){          ## Part on nblast
        #browser()
        if (probabilities_nblastscores_kcs[k,l,findInterval(scores_neuron[k],scorecut)]==0){
          scores_neurons_families_nblast[names(listofneurons)[i],l] =  scores_neurons_families_nblast[names(listofneurons)[i],l] -0.9
        }else{
          scores_neurons_families_nblast[names(listofneurons)[i],l] =  scores_neurons_families_nblast[names(listofneurons)[i],l] +
            log(probabilities_nblastscores_kcs[k,l,
                                               findInterval(scores_neuron[k],scorecut)])
        }
      }
      scores_neurons_families_nblast[names(listofneurons)[i],l] = scores_neurons_families_nblast[names(listofneurons)[i],l]/length(subfamilies)

      ## Part on supervoxels
      prob = probability_sv_knowing_subfamily[sv_neuron,l]
      zeros = sum(prob==0)
      prob2 = prob[which(prob!=0)]
      score_f = sum(log(prob2))
      score_f = score_f -100*zeros
      scores_neurons_families_nblast[names(listofneurons)[i],l] =  scores_neurons_families_nblast[names(listofneurons)[i],l] + score_f/length(sv_neuron)+ log(probabily_subfamily_kcs[l])
    }
  }
  return(scores_neurons_families_nblast)
}

prior_prob_nblastscores = function(subfamilies = kcsm, bins=seq(-1,1,0.2)){
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
        scoren = nblast(neuron,dps[names(names3[l])])/(max(smat.fcwb)*nrow(neuron[[1]]$points))   ## Score between the neuron and the neuron l of the family j, normalized
        probabilities_nblastscores[i,j,findInterval(scoren,bins)] = probabilities_nblastscores[i,j,findInterval(scoren,bins)]+1
      }
    }
  }
  probabilities_nblastscores = probabilities_nblastscores/length(subfamilies)
  return(probabilities_nblastscores)
}

## subfamiliesall should be the list of all the neurons in
prior_prob_subfam = function(subfamiliesall = kcs){
  probabily_subfamily = c()
  for(i in unique(subfamiliesall)){
    print(i)
    probabily_subfamily = c(probabily_subfamily, table(subfamiliesall)[i]/sum(table(subfamiliesall)))
  }
  names(probabily_subfamily) =  unique(subfamiliesall)
  return(probabily_subfamily)
}

prior_prob_svscores = function(subfamilies = kcsm,computedens = FALSE){
  if (computedens == TRUE){
    voxel_dens_allneurons = voxel_dens(names(subfamilies))
    rownames(voxel_dens_allneurons) = names(subfamilies)
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
