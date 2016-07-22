
load("/Volumes/JData5/JPeople/Melina/Branson/data/fc_neuron_typec")                     ### All the neurons with the correct labels we are working on
load("/Volumes/JData5/JPeople/Melina/Branson/data/correct_families")                    ### All 137 families, we still have to get rid of the ones that contains one neuron
load("/Volumes/JData5/JPeople/Melina/Branson/data/probability_correct_families")        ### Probability to have one particular family
load("/Volumes/JData5/JPeople/Melina/Branson/data/number_neurons_fromfam_insv")         ### The numbes of neurons crossing each supervoxel, matrix of 706
load("/Volumes/JData5/JPeople/Melina/Branson/data/probability_sv_knowing_family")       ### The probabilities of all the supervoxels in each family, matrix of nrow = 7065

# we don't really need this since it will be evident from other objects
load(file="/Volumes/JData5/JPeople/Melina/Branson/data/names_svoxels") # Names of the supervoxels
probability_sv_knowing_family=round(probability_sv_knowing_family,3)
rownames(probability_sv_knowing_family) = names(svoxels.fcwb.table)[-1]
colnames(probability_sv_knowing_family)=names(correct_families)

probability_correct_families=c(probability_correct_families)
names(probability_correct_families)=names(correct_families)

library(devtools)
use_data(probability_sv_knowing_family, overwrite = T)
use_data(number_neurons_fromfam_insv)
use_data(probability_correct_families, overwrite = T)
use_data(correct_families)
use_data(fc_neuron_typec)
