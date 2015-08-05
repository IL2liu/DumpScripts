#!/usr/bin/env Rscript

#########################################################
# input:                                                #
# - dataframe: in a bed format                          #
# - e.g chr start end   profile1 profile2 profile..     #
#        1  1     3000   20       20       20           #
# output:                                               #
# - Robjects consists of                                #
# - 1. f : input dataframe                              #
# - 2. pca_f : pca object                               #
#########################################################

# commandline parser
suppressPackageStartupMessages(library("optparse"))

option_list <- list(
      make_option(c("-r", "--remove"), action="store_true", default=FALSE,
                  help="Remove rows that contain only zeros [default]"),
      make_option(c("-o", "--fname"), type="character",
                  help="add the file name for your PCA object")

      )

parser <- OptionParser(usage = "%prog [options] matrix", option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)

opt <- arguments$options

if(length(arguments$args) != 1) {

cat("Incorrect number of required positional arguments\n\n")

print_help(parser)

stop()

} else {

file <- arguments$args

}

if( file.access(file) == -1) {

stop(sprintf("Specified file ( %s ) does not exist", file))

} else {

f <- read.delim(file)
# extract the matrix without coordinates
# chr start end prof1....
mat_f <- f[,4:ncol(f)]

}

if(opt$remove) {
   mat_f[which(rowSums(mat_f)==0),] <- NA

   #remove the genomic regions that contain NA                   
   mat_f <- mat_f[complete.cases(mat_f),]

} 

# logify the matrix
log_mat <- f <- log(mat_f + 1)

# pca analysis
pca_f <- prcomp(log_mat, scale.=T)

# save the R objects
save(pca_f,f, file=opt$fname)
