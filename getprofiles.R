### R code from vignette source 'getprofiles.Rnw'

###################################################
### code chunk number 1: getprofiles.Rnw:21-27
###################################################
dir("data")
load("data/158_profiles.rda")
ls()
head(profiles[,1:6])
key <- read.delim("data/key.txt")
head(key)


###################################################
### code chunk number 2: getprofiles.Rnw:43-49
###################################################
reads <- as.matrix(profiles[,4:ncol(profiles)])
maps <- profiles[,1:3]
rm(profiles)
sum(reads)
blacklist <- rowSums(reads) == 0
reads[blacklist,] <- NA


###################################################
### code chunk number 3: getprofiles.Rnw:58-59
###################################################
superinput <- rowSums(reads[,c("X008", "X144", "X195")])


###################################################
### code chunk number 4: getprofiles.Rnw:64-86
###################################################
load('output/fitlist_mixture.rda')
colrs = c('#1A1A1A7F', '#ED12887F')
for (feature in c('CHD1', 'H3K9me3', 'CTCF')) {
   fit = fitlist[[feature]]
   pdf(paste(feature, ".pdf", sep=""), width=10, height=5)
   HMMdata <- NULL
   HMMdata <- data.frame(seqname=maps[,1], input=superinput,
      reads[,colnames(reads) %in%
      key$ID[which(key$replicate == feature)]])
   for (i in 3:ncol(HMMdata)) {
      plot(.003*(1:5000), HMMdata[1:5000,i], type="h",
         col=colrs[fit$vPath+1], main=feature,
         ylab="Read count", xlab="Position on chr1 (Mb)")
   }
   dev.off()
}
