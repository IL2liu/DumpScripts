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
   pdf(paste("~/Dropbox/FastBW/", feature, ".pdf", sep=""), width=12)
   HMMdata <- NULL
   HMMdata <- data.frame(seqname=maps[,1], input=superinput,
      reads[,colnames(reads) %in%
      key$ID[which(key$replicate == feature)]])
   if (ncol(HMMdata) < 3) {
      dev.off()
      next
   }
#   Q <- matrix(c(
#      0.892385,0.023884,0.060823,
#      0.070055,0.915551,0.073261,
#      0.037560,0.060565,0.865915), nrow=3)
#   alpha <- 0.995256
#   beta <- 44.064042
#   gammas <- matrix(c(
#      0.177337,0.056239,0.247651,
#      0.557039,0.205209,0.168496,
#      4.232230,1.366844,0.753653), nrow=3)
   Q <- matrix(c(
      0.946357,0.031118,0.022849,
      0.028076,0.940256,0.033188,
      0.025566,0.028626,0.943963), nrow=3)
   alpha = 0.911266
   beta = 60.524137
   gammas = matrix(c(
      0.264054,0.086798,0.059893,
      0.552719,0.187703,0.156656,
      4.187121,1.406957,0.719832), nrow=3)
   fit <- BaumWelch.NB(HMMdata, Q=Q, alpha=alpha, beta=beta, gammas=gammas)
   if (is.null(fit)) {
      dev.off()
      next
   }
   for (i in 3:ncol(HMMdata)) {
      plot(HMMdata[1:5000,i], type="h", col=fit$vPath+1,
         main=colnames(HMMdata)[i])
   }
   dev.off()
   fitlist[[feature]] <- fit
}
save(fitlist, file='output/fitlist.rda')
