library(Matrix)
voxel_dens_allneurons=Matrix(voxel_dens_allneurons, sparse=T)

library(devtools)
use_data(voxel_dens_allneurons)
