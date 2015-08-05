HMMdata = read.delim('discretized_with_NAs.txt.gz')
library(chromHMMatin)
fit5 = BaumWelch(HMMdata, m=5)
save(fit5, file=paste('fit5', format(Sys.time(), "%s"), ".rda", sep=''))
