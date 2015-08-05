### R code from vignette source 'discretize_H1.Rnw'

###################################################
### code chunk number 1: discretize_H1.Rnw:21-27
###################################################
dir("data")
load("data/158_profiles.rda")
ls()
head(profiles[,1:6])
key <- read.delim("data/key.txt")
head(key)


###################################################
### code chunk number 2: discretize_H1.Rnw:43-49
###################################################
reads <- as.matrix(profiles[,4:ncol(profiles)])
maps <- profiles[,1:3]
rm(profiles)
sum(reads)
blacklist <- rowSums(reads) == 0
reads[blacklist,] <- NA


###################################################
### code chunk number 3: discretize_H1.Rnw:58-59
###################################################
superinput <- rowSums(reads[,c("X008", "X144", "X195")])
save(superinput, file='superinput.rda')
