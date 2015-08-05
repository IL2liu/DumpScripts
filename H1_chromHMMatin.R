##############################
# run HMM to color chromatin #
##############################

# library used
library(chromHMMatin)

# set wd
setwd("3000-mappedProfiles_wNA/")

# load the dataset
d_ <- read.delim("AllDiscretized_bigTable.txt")
