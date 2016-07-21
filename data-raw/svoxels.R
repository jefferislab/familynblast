library(nat)
svoxels.fcwb=read.im3d("/Volumes/JData5/JPeople/Melina/Branson/JFR/FCWB_AnatomySubCompartments20150108_ms7065centers.nrrd")
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
library(devtools)
use_data(svoxels.fcwb.table)
