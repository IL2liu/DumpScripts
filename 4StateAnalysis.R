

# set global chunk options
# for figures
opts_chunk$set(fig.path='figs/', fig.align='center', fig.show='hold',dev='CairoPDF', out.width='.4\\linewidth')
# replacing "=" into "->" to make it R thing
options(replace.assign=TRUE,width=90)
# caching chunks
opts_chunk$set(cache.extra = R.version,cache.path='cache/')
opts_chunk$set(cache.extra = rand_seed)



load('input/fit41366702728.rda')
ls()
print (names(fit4))



## em <- fit4$B
## head(em)
## colnames(em) <- c('Red', 'Pink', 'Yellow', 'Black')
## 
## library(gplots)
## get_heatMat <- function(mat_, margins, Colv, Rowv, cutZero){
##   # construct heatmap, given a matrix
##   heatmap.2(mat_, col=cm.colors(256),
##             Colv=Colv, Rowv=Rowv,
##             dendrogram="none", trace="none",
##             density.info = 'none',scale="none",
##              keysize=0.8, margins = margins,
##             symbreaks=cutZero)
## }
## 
## pdf('figs/em_matrix.pdf', useDingbats = FALSE)
## HeatMapEm <- get_heatMat(em, margins=c(5,10), Colv=NA, Rowv=TRUE, cutZero=FALSE)
## dev.off()
## 



## trans <- fit4$Q
## colnames(trans) <- c('Red', 'Pink', 'Yellow', 'Black')
## rownames(trans) <- c('Red', 'Pink', 'Yellow', 'Black')
## 
## pdf('figs/trans_matrix.pdf', useDingbats=FALSE)
## HeatMapTrans <- get_heatMat(trans, margins=c(5,5), Colv=NA, Rowv=NA, cutZero=FALSE)
## dev.off()
## 



## 
## get_enrichmentPerGenomeAverage <- function(em_mat, fit_object){
## 
##   # input: emmision matrix
##   # fit_object: fit object run from BaumWelch() function in chromHMMatin
##   averagePerState <- table(fit_object$vPath)/sum(table(fit_object$vPath))
##   marksPerGenomeAverage <- em_mat %*% averagePerState
##   mat_perGenomeAverage <- sweep(em_mat, 1,  marksPerGenomeAverage, "/" )
##   return(mat_perGenomeAverage)
## 
## }
## 
## em_GenomeAverage <- get_enrichmentPerGenomeAverage(em, fit4)
## get_heatMat(log(em_GenomeAverage), Colv=NA, Rowv=TRUE,margins=c(5,10), cutZero=TRUE)
## 
## #thresholding the enrichment for 10 profiles
## em_GenomeAverage_thresholded <- em_GenomeAverage
## em_GenomeAverage_thresholded[em_GenomeAverage_thresholded < 1e-6] <- 1e-6
## 
## pdf('figs/enrich_matrix.pdf', useDingbats=FALSE)
## get_heatMat(log(em_GenomeAverage_thresholded), Colv=NA, Rowv=TRUE,margins=c(5,10), cutZero=TRUE)
## dev.off()
## 



## coverage_data <- read.delim('input/discretized_with_NAs.txt.gz')
## head(coverage_data)
## coverage_data[1:5,1:5]
## dim(coverage_data)
## # get the color state
## state <- fit4$vPath
## # exclude the chr column
## coverage_data <- coverage_data[, c(2:69)]
## # combine it with the color state
## coverage_data <- cbind(coverage_data, state)
## # remove the NAs
## coverage_data <- coverage_data[complete.cases(coverage_data), ]
## # get only the profiles without the color state
## coverage_mat <- coverage_data[, c(1:68)]
## # compute the genome-wide average of each profile
## GW <- apply(coverage_mat, 2, mean)
## 
## library(plyr)
## 
## # compute the average of each state for each profile
## mean_state <- aggregate(.~state, data=coverage_data, FUN=mean)
## head(mean_state)
## mean_state[1:4,1:4]
## # exclude the first column -> the state's info
## coverage_state <- mean_state[,2:ncol(mean_state)]
## # combine it with the genome-wide average
## coverage_profiles <- rbind(coverage_state, GW)
## head(coverage_profiles)
## rownames(coverage_profiles) <- c('Red','Pink', 'Yellow','Black',"GW")
## coverage_profiles <- t(coverage_profiles)
## 
## pdf('figs/coverage_matrix.pdf', useDingbats=FALSE)
## get_heatMat(coverage_profiles, Colv=NA, Rowv=TRUE, margins=c(4,10), cutZero=FALSE)
## dev.off()



sessionInfo()



library(knitr)
knit("4StateAnalysis.Rnw" ) # compile to tex
purl("4StateAnalysis.Rnw", documentation = 0) # extract R code only
knit2pdf("4StateAnalysis.Rnw")


