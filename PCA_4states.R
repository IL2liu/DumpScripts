

# set global chunk options
# for figures
opts_chunk$set(fig.path='figs/', fig.align='center', fig.show='hold',dev='CairoPDF', out.width='.4\\linewidth')
# replacing "=" into "->" to make it R thing
options(replace.assign=TRUE,width=90)
# caching chunks
opts_chunk$set(cache.extra = R.version,cache.path='cache/')
opts_chunk$set(cache.extra = rand_seed)



## H1_colorProfiles <- read.delim('input/4states_H1Profiles.bed')
## 
## # subset 200000 bins randomly
## Rand_RowIdx <- sample(nrow(H1_colorProfiles), size=200000)
## head(Rand_RowIdx)
## subsetH1_colorProfiles <- H1_colorProfiles[Rand_RowIdx, ]
## subsetH1_colorProfiles[1:3, 1:5]
## dim(subsetH1_colorProfiles)
## 



## maps <- subsetH1_colorProfiles[, c(1,2,3)]
## colorProfiles <- subsetH1_colorProfiles[,4:162]
## colorStates <- subsetH1_colorProfiles[,163]
## 
## # log the reads
## Log_colorProfiles <- log(colorProfiles + 1)



## pca_colorProfiles <- prcomp(Log_colorProfiles, scale.=T)
## library(rgl)
## colorStates <- c('red', 'yellow', 'pink', 'black')[as.factor(colorStates)]
## plot3d(pca_colorProfiles$x[, 1:3], col=colorStates)
## 



animate_pca <- function(dir){
  
  dir.create(dir, showWarnings = FALSE)
  f=spin3d(c(1,1,1), rpm=.5)
  for (i in 1:120) {
    view3d(userMatrix=f(i)$userMatrix)
    filename <- paste("4states",formatC(i,digits=3,flag="0"),".png",sep="")
    rgl.snapshot(file.path(dir,filename))
  }
  clear3d()
}

get_3dPlot <- function(pca_scores, col){
  rgl.open()
  rgl.points(scale(pca_scores), col=col, bg='white')
  box3d(col='black')
  axes3d(c())
}

#animate PCA from 91 profiles
get_3dPlot(pca_colorProfiles$x[, 1:3], col=colorStates)
animate_pca("output_animation")





sessionInfo()



library(knitr)
knit("PCA_4states.Rnw" ) # compile to tex
purl("PCA_4states.Rnw", documentation = 0) # extract R code only
knit2pdf("PCA_4states.Rnw")


