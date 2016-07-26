### Contains the array that has, the families in columns,
load(file="/Volumes/JData5/JPeople/Melina/Branson/data/probabilities_nblastscores_kcs")
#the neurons in lines and in depth the value of the score
load(file="/Volumes/JData5/JPeople/Melina/Branson/data/probability_sv_knowing_subfamily")
rownames(probability_sv_knowing_subfamily) = colnames(voxel_dens_allneurons)[-1]
load("/Volumes/JData5/JPeople/Melina/Branson/data/probabily_subfamily_kcs")
# test_set_kenyoncells ----------------------------------------------------
library(flycircuit)
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
kcsm = unlist(kcslistm)
names_kcsm=lapply(kcslistm, names)
names(names_kcsm) = sapply(kcslistm, unique)

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
kcs = unlist(kcslist)
names_kcs = lapply(kcslist, names)
names(names_kcs) = sapply(kcslist, unique)
Lab = c("G","ApBp","ABp","ABc","ABs","Gd")
names(Lab) = names(names_kcsm)


kcs.subfam.training=kcsm
kcs.subfam.test=kcs
use_data(kcs.subfam.training)
use_data(kcs.subfam.test)

# Load dps object with all flycircuit neurons
# see https://gist.github.com/jefferis/bbaf5d53353b3944c090

library(devtools)
devtools::source_gist("bbaf5d53353b3944c090")

# compute priors
probabilities_nblastscores_kcs = prior_prob_nblastscores()
use_data(probabilities_nblastscores_kcs)

probability_sv_knowing_subfamily =prior_prob_svscores()
use_data(probability_sv_knowing_subfamily)

probability_subfamily_kcs = prior_prob_subfam()
use_data(probability_subfamily_kcs)
