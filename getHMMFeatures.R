#!/usr/bin/env Rscript

###################################################
# get the HMM features from the HMM in .rda files #
###################################################

# commandline parser
suppressPackageStartupMessages(library("optparse"))                            

option_list <- list(
     make_option(c("-f", "--feature"), type="character",
                 help="add the HMM feature you'd like to extract"),
     make_option(c("-i", "--input"), type="character",
                 help="input directory"), 
     make_option(c("-p", "--pattern"), type="character",
                 help="pattern for files"), 
     make_option(c("-o", "--output"), type="character",
                 help="output directory")
     )

parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options

getFeature <- function(rda, feature){
   env <- new.env()
   nm <- load(rda, env)
   HMM_feature <- env[[nm]][[feature]]
   return (HMM_feature)
}

for (fname in list.files(path = opt$input, glob2rx(pattern= opt$pattern))){
   feature <- getFeature(rda = paste(opt$input, fname, sep='/'), feature = opt$feature)
   # get the basename without extension
   f_name <- sub("^([^.]*).*", "\\1", fname)
   out_f <- paste(opt$output, f_name, sep='/')
   write.table(feature, file=paste(out_f, '.ls', sep=''), 
               col.names=f_name, row.names=FALSE, 
               quote=FALSE )
}
