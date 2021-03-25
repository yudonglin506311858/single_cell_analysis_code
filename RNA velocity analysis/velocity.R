conda create -n velocyto
conda activate velocyto
conda install R
conda install -c conda-forge openmp
conda install -c anaconda boost

#R environment
install.packages("devtools")

library(devtools)
install_github("velocyto-team/velocyto.R")
library(velocyto.R)
