library(Matrix)
load("/Volumes/JData5/JPeople/Melina/Branson/data/voxel_dens_allneurons")
voxel_dens_allneurons=Matrix(voxel_dens_allneurons, sparse=T)

library(devtools)
use_data(voxel_dens_allneurons, overwrite = T)
