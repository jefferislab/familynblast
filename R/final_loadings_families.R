library("nat.flybrains")
library(nat)
library(flycircuit)
library(nat.nblast)

### Pour trouver jet.colors il faut regarder la défintion dans colorampalette 
#----------------------------------
# fc_download_data('http://flybrain.mrc-lmb.cam.ac.uk/si/nblast/flycircuit/allbyallblastcv4.5.ff',
#                  type='ff')
# # set that as default all by all score matrix
# options('flycircuit.scoremat'="allbyallblastcv4.5.ff")
# # load neuron list
# # the actual neuron data will be downloaded and cached to your machine on demand
# dps<-read.neuronlistfh("http://flybrain.mrc-lmb.cam.ac.uk/si/nblast/flycircuit/dpscanon.rds",
#                        localdir=getOption('flycircuit.datadir'))
# remotesync(dps,download.missing=T)

# set default neuronlist
options('nat.default.neuronlist'='dps')
#load("/Volumes/JData5/JPeople/Melina/Branson/data/all_indices_voxels")                  ### All the indices of the supervoxels
load("/Volumes/JData5/JPeople/Melina/Branson/data/voxel_dens_allneurons")               ### All the number of points from the neurons crossing the supervoxels
voxel_dens_allneurons = voxel_dens_allneurons[,2:7066]
load("/Volumes/JData5/JPeople/Melina/Branson/data/fc_neuron_typec")                     ### All the neurons with the correct labels we are working on
load("/Volumes/JData5/JPeople/Melina/Branson/data/correct_families")                    ### All 137 families, we still have to get rid of the ones that contains one neuron
load("/Volumes/JData5/JPeople/Melina/Branson/data/probability_correct_families")        ### Probability to have one particular family 
load("/Volumes/JData5/JPeople/Melina/Branson/data/number_neurons_fromfam_insv")         ### The numbes of neurons crossing each supervoxel, matrix of 706
load("/Volumes/JData5/JPeople/Melina/Branson/data/probability_sv_knowing_family")       ### The probabilities of all the supervoxels in each family, matrix of nrow = 7065
load(file="/Volumes/JData5/JPeople/Melina/Branson/data/names_svoxels") # Names of the supervoxels
probability_sv_knowing_family=round(probability_sv_knowing_family,3)
rownames(probability_sv_knowing_family) = names_svoxels[-1]
svoxels.fcwb=read.im3d("/Users/Melina/Documents/Stage/Scripts/Funs/FCWB_AnatomySubCompartments20150108_ms7065centers.nrrd")
quicktable <- function(x) {
  xname=deparse(substitute(x))
  tt=tabulate(x+1)
  levels=seq.int(from=0, length.out = length(tt))
  nz=tt!=0L
  
  structure(tt[nz], .Dim = sum(nz), 
            .Dimnames = structure(list(as.character(levels[nz])), .Names = xname),
            class = "table")
}
svoxels.fcwb.table=quicktable(svoxels.fcwb)
nsvoxels=length(svoxels.fcwb.table)


### voxel dens computes the number of points from a setofneurons in each supervoxel
voxel_dens= function(setofneurons){
  svoxel.output.neu = matrix(0, nrow = length(setofneurons), ncol = nsvoxels)
  colnames(svoxel.output.neu)=names(svoxels.fcwb.table)
  rownames(svoxel.output.neu)=names(setofneurons)
  for(n in names(setofneurons)) {
    message("working on neuron:", n)
    # find 1d indices of 3D coords into FCWB space
    idxs=coord2ind(xyzmatrix(setofneurons[[n]]), imdims = as.im3d(FCWB))
    # then get label values for those coords
    svoxel.ids=svoxels.fcwb[idxs]
    # compute density i.e. number of points in each domain
    tsvoxel.ids=table(svoxel.ids)
    svoxel.output.neu[n,names(tsvoxel.ids)]=tsvoxel.ids
  }
  return(svoxel.output.neu)  
}

### computing the score of one neuron against a particular family, given the index of the name of the family in correct_families
neurons_against_fam  = function(listneurons,familyind, computedens = FALSE,zeroscore = -100){
  if (computedens == TRUE){
    voxel_dens_allneurons = voxel_dens(listneurons)
    rownames(voxel_dens_allneurons) = listneurons
  }
  Score =c()
  for (k in seq_along(listneurons)){                                        ## for each neuron
    print(k)
    neuron = listneurons[k]                  
    sv_neurons = names(voxel_dens_allneurons[names(neuron),voxel_dens_allneurons[names(neuron),]>0])    
    prob = probability_sv_knowing_family[sv_neurons,familyind]
    zeros = sum(prob==0)
    prob2 = prob[which(prob!=0)]
    score_f = sum(log(prob2))
    score_f = score_f+zeroscore*zeros
    Score = c(Score,score_f)
    }
    Score = Score/length(sv_neurons)+log(probability_correct_families[familyind])
  return(Score)
}

### find_scores_family intermediate function that contains all the scores of the list ofneurons against the families
find_scores_family = function(listneurons,computedens = FALSE,zeroscore = -100){
  matrix_scores_neurons = matrix(0,nrow =listneurons, ncol = length(correct_families))
  for(k in 1:length(correct_families)){
    matrix_scores_neurons[,k] = neurons_against_fam(listneurons,familyind = k,computedens=FALSE, zeroscore=-100)
  }
  
  list_scores_neurons = list()
  for(l in 1:length(listneurons)){
    list_scores_neurons = append(list_scores_neurons,matrix_scores_neurons[l])
  }
  names(list_scores_neurons) = names(listneurons)
  for(i in seq_along(list_scores_neurons)){
    names(list_scores_neurons[[i]])=names(correct_families)
  }
  return(list_scores_neurons)
}

#### find_neurons_family contains the list of families associated to the list of neurons given
find_neurons_family = function(listneurons,save = FALSE, path="",computedens = FALSE){
  familiesof_listofneurons = c()
  list_scores_neurons = find_scores_family(listneurons,save = FALSE, path="",computedens=FALSE)
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

### find the percentage of correct hits within the first nb highest scores 
find_percentage_correct_hits = function(listneurons,nb = 3){
  Percents = c()
  scorestouse = list_scores_neurons_cv_fun 
  for(j in 1:nb){
  for(l in seq_along(scorestouse)){
      name_hit=names(scorestouse[[l]])[rev(order(scorestouse[[l]]))[1:j]]
      if(length(intersect(neuron[[1]],name_hit))>0){
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

### find the list of scores of all the neurons for the cross validation approach !! 
list_scores_neurons_cv_fun = function(listneurons = fc_neuron_typec, computedens = FALSE,zeroscore = -100 ){
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
    correct_familiesb = correct_families
    correct_familiesb[[neuron[[1]]]] = correct_families[[neuron[[1]]]][names(correct_families[[neuron[[1]]]])!=names(neuron)]
    
    ### Process with the families
    probability_correct_families_b = probability_correct_families
    indx = which(names(correct_familiesb)==neuron[[1]]) 
    probability_correct_families_b[,indx]=length(correct_familiesb[[indx]])/length(fc_neuron_typecb)
    
    ### Process to compute the number of neuron in the family 
    number_neurons_fromfam_insv_b = number_neurons_fromfam_insv
    
    if(length(ncol(voxel_dens_allneurons[names(correct_familiesb[[indx]]),]!=0))==0){
      number_neurons_fromfam_insv_b[,indx] = sum(voxel_dens_allneurons[names(correct_familiesb[[indx]]),]!=0)
    }else{
      number_neurons_fromfam_insv_b[,indx] = colSums(voxel_dens_allneurons[names(correct_familiesb[[indx]]),]!=0)
    }
    
    probability_sv_knowing_familyb = probability_sv_knowing_family
    probability_sv_knowing_familyb[,indx] = number_neurons_fromfam_insv_b[,indx]/length(correct_families[[indx]])
    
    
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
  names(list_scores_neurons_cv) = names(fc_neuron_typec)
  
  for(i in seq_along(list_scores_neurons_cv)){
    names(list_scores_neurons_cv[[i]])=names(correct_families)
  }
  return(list_scores_neurons_cv)
}

### Create the probability matrices 
### The set of families should be of the form: list of families and in each list the list of neurons in this family. An example is the "correct_families"
create_probab_families = function(setoffamilies){
  probability_correct_families = matrix(0,nrow=1,ncol=length(setoffamilies)) 
  for (i in seq_along(setoffamilies)){
    allneu = sum(sapply(1:length(setoffamilies), function(n) length(setoffamilies[[n]])))
    probability_correct_families[1,i]= length(setoffamilies[[i]])/allneu
  } 
return(probability_correct_families)
}

### create the probability matrix that contains the supervoxels in rows and families in columns for a setoffamilies
create_probab_sv_knowing_fam = function(setoffamilies,computedens=FALSE){
  if (computedens == TRUE){
    voxel_dens_allneurons = voxel_dens(listneurons)
    rownames(voxel_dens_allneurons) = listneurons
  }
  number_neurons_fromfam_insv = matrix(0,nrow=7065,ncol=length(setoffamilies))
  for (i in seq_along(setoffamilies)){
    print(paste("moving to family",i))
    setofneurons = setoffamilies[[i]]   ### ??
    for(j in 1:7065){
      number_neurons_fromfam_insv[j,i] = sum(voxel_dens_allneurons[names(setofneurons),j]!=0)
    }
  }
  # Filling up the matrix of probabilities for supervoxels knowing the family --------
  ### We have to determine the probability of having sj knowing we are in the family i 
  probability_sv_knowing_family = matrix(0,nrow=7065,ncol=length(setoffamilies)) 
  for(l in seq_along(setoffamilies)){
    for (k in 1:7065){
      probability_sv_knowing_family[k,l] = number_neurons_fromfam_insv[k,l]/length(correct_families[[l]])
    }
  }
  return(probability_sv_knowing_family)
}

