load("data/158_profiles.rda")
key <- read.delim("data/key.txt")
reads <- as.matrix(profiles[,4:ncol(profiles)])
maps <- profiles[,1:3]
rm(profiles)
blacklist <- rowSums(reads) == 0
reads[blacklist,] <- NA
superinput <- rowSums(reads[,c("X008", "X144", "X195")])
library(HummingBee)
fitlist <- list()
#for (feature in unique(key$replicate)) {
for (feature in 'SUZ12') {
   HMMdata <- data.frame(seqname=maps[,1], input=superinput,
      reads[,colnames(reads) %in%
      key$ID[which(key$replicate == feature)]])
   print(mlNB(HMMdata[,2]))
}
