

# set global chunk options
# for figures
opts_chunk$set(fig.path='figs/', fig.align='center', fig.show='hold',dev='CairoPDF', out.width='.4\\linewidth')
# replacing "=" into "->" to make it R thing
options(replace.assign=TRUE,width=90)
# caching chunks
opts_chunk$set(cache.extra = R.version,cache.path='cache/')
opts_chunk$set(cache.extra = rand_seed)



## color4Dataset <- read.delim('input/H1_hESC_4ColorStates.txt')
## dim(color4Dataset)
## head(color4Dataset)[,1:5]
## tail(color4Dataset)[,1:5]
## mat4Dataset <- color4Dataset[,2:70]
## head(mat4Dataset)



## mat4DatasetNoNa <- mat4Dataset[complete.cases(mat4Dataset),]
## head(mat4DatasetNoNa)
## dim(mat4DatasetNoNa)
## nrow(mat4DatasetNoNa)/nrow(mat4Dataset) #85% of the whole genomic loci could be used



## mat4DatasetNoNa <- as.data.frame(mat4DatasetNoNa)
## # select for color 0
## 
## selectSet <- function(input, state, sampling=FALSE, sample.no){
##   State <- input[input$Color==state, ]
##   NonState <- input[input$Color!=state, ]
##   NonState$Color <- state+1
##   if (sampling){
##     State <- State[sample(nrow(State), size = sample.no, replace=FALSE),]
##     NonState <- NonState[sample(nrow(NonState), size = sample.no, replace=FALSE),]
##     Dataset <- rbind(State, NonState)
##     Dataset$Color <- as.factor(Dataset$Color)
##     return (Dataset)
##   }
##   Dataset <- rbind(State, NonState)
##   Dataset$Color <- as.factor(Dataset$Color)
##   return (Dataset)
## }
## 
## set0 <- selectSet(mat4DatasetNoNa, state = 0, sampling=TRUE, sample.no = 1000)
## set1 <- selectSet(mat4DatasetNoNa, state = 1, sampling=TRUE, sample.no = 1000)
## set2 <- selectSet(mat4DatasetNoNa, state = 2, sampling=TRUE, sample.no = 1000)
## set3 <- selectSet(mat4DatasetNoNa, state = 3, sampling=TRUE, sample.no = 1000)
## 
## dim(set3)
## head(set3)
## 
## 



library(randomForest)
set.seed(19)


RF_Profiles <- function(dataset){

  stateClassifier <- randomForest(x = dataset[, 1:68], y = dataset[,69], importance=TRUE, proximity=TRUE)
  profile_importance <- (stateClassifier$importance)
  profile_importance <- profile_importance[, "MeanDecreaseAccuracy"]
  profiles <- names(sort(profile_importance, decreasing =T)[1:10]) 
  return (profiles)

}

Profiles_set0 <- RF_Profiles(set0)
Profiles_set1 <- RF_Profiles(set1)
Profiles_set2 <- RF_Profiles(set2)
Profiles_set3 <- RF_Profiles(set3)

table(mat4DatasetNoNa$Color)/sum(table(mat4DatasetNoNa$Color))
Profiles_set0
Profiles_set1
Profiles_set2
Profiles_set3

# state 0 versus state2
State1_vs_State2Set <- function(input, state1, state2, sampling=FALSE, sample.no){
  State <- input[input$Color==state1, ]
  NonState <- input[input$Color==state2, ]
  if (sampling){
    State <- State[sample(nrow(State), size = sample.no, replace=FALSE),]
    NonState <- NonState[sample(nrow(NonState), size = sample.no, replace=FALSE),]
    Dataset <- rbind(State, NonState)
    Dataset$Color <- as.factor(Dataset$Color)
    return (Dataset)
  }
  Dataset <- rbind(State, NonState)
  Dataset$Color <- as.factor(Dataset$Color)
  return (Dataset)
}

PromEnhan <- State1_vs_State2Set(mat4DatasetNoNa, state1=0, state2=2, sampling=TRUE, sample.no=1000)
Profiles_PromEnhan <- RF_Profiles(PromEnhan)
Profiles_PromEnhan

PromEnhanClassifier <- randomForest(x = PromEnhan[, 1:68], y = PromEnhan[,69], importance=TRUE, proximity=TRUE)
varImpPlot(PromEnhanClassifier, main="Variable Importance:RF")




sessionInfo()



library(knitr)
knit("RandomForestClassifier.Rnw" ) # compile to tex
purl("RandomForestClassifier.Rnw", documentation = 0) # extract R code only
knit2pdf("RandomForestClassifier.Rnw")


