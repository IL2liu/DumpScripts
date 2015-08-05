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
done = scan('data/done_already.txt', what="")
for (feature in unique(key$replicate)) {
   if (feature %in% done) next
   cat(file=stderr(), paste(feature, "\n"))
   pdf(paste("output/", feature, ".pdf", sep=""), width=12)
   HMMdata <- NULL
   HMMdata <- data.frame(seqname=maps[,1], input=superinput,
      reads[,colnames(reads) %in%
      key$ID[which(key$replicate == feature)]])
   if (ncol(HMMdata) < 3) {
      dev.off()
      next
   }
   fit <- HummingBee(HMMdata)
   for (i in 3:ncol(HMMdata)) {
      plot(HMMdata[1:5000,i], type="h", col=fit$vPath+1,
         main=colnames(HMMdata)[i])
   }
   dev.off()
   fitlist[[feature]] <- fit
}
save(fitlist, file='output/fitlist.rda')

#for (feature in unique(key$replicate)) {
#for (feature in 'RXRA') {
   #pdf(paste("~/Dropbox/FastBW/", feature, ".pdf", sep=""), width=12)
#   HMMdata <- NULL
#   HMMdata <- data.frame(seqname=maps[,1], input=superinput,
#      reads[,colnames(reads) %in%
#      key$ID[which(key$replicate == feature)]])
#   write.table(HMMdata, file='RXRA.txt', sep='\t', quote=F,
#   row.names=F)
   #if (ncol(HMMdata) < 3) {
   #   dev.off()
   #   next
   #}
   #Q <- matrix(.025, nrow=3, ncol=3)
   #diag(Q) <- .95
   #fit <- HummingBee(HMMdata)
   #fit <- BaumWelch.NB(HMMdata, Q=fit$Q, alpha=fit$alpha,
   #   beta=fit$beta, gammas=fit$gamma)
   #if (is.null(fit)) {
   #   dev.off()
   #   next
   #}
   #for (i in 3:ncol(HMMdata)) {
   #   plot(HMMdata[1:5000,i], type="h", col=fit$vPath+1,
   #      main=colnames(HMMdata)[i])
   #}
   #dev.off()
   #fitlist[[feature]] <- fit
#}
