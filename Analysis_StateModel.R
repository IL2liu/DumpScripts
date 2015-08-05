

# set global chunk options
# for figures
opts_chunk$set(fig.path='figs/', fig.align='center', fig.show='hold',dev='CairoPDF', out.width='.4\\linewidth')
# replacing "=" into "->" to make it R thing
options(replace.assign=TRUE,width=90)
# caching chunks
opts_chunk$set(cache.extra = R.version,cache.path='cache/')
opts_chunk$set(cache.extra = rand_seed)



LogLik_4state <- read.delim('outputLogLik_4states/loglik_4states.ls')
LogLik_5state <- read.delim('outputLogLik_5states/loglik_5states.ls')

# get the first 5 runs
state4_top5 <- as.numeric(LogLik_4state)[1:5]
state5_top5 <- as.numeric(LogLik_5state)[1:5]

box_state <- rbind(state4_top5, state5_top5)
box_state <- t(box_state)

# draw the histograms
h <- hist(box_state, plot = FALSE)

# color each bin according to the box_state(4-state or 5-state)
latest.ob_4state <- box_state[, 1]
latest.ob_5state <- box_state[,2]
bin_4state <- cut(latest.ob_4state, h$breaks)
bin_5state <- cut(latest.ob_5state, h$breaks)
clr <- rep("white", length(h$counts))
clr[bin_4state] <- "green"
clr[bin_5state] <- "red"

# producing the figure
pdf("figs/log_likelihood.pdf", useDingbats=FALSE)
plot(h, col=clr, xlab="", main="Profile Log-Likelihood")
legend('topright', fill=c('green', 'red'), c('Four-state HMM', 'Five-state HMM'),
       bg='transparent', box.lwd=0, box.col="transparent")
dev.off()




vPath_4state <- read.delim('output_4states/vPath_4states.ls')
freq_table4 <- apply(vPath_4state, 2, function(x)table(x))
library(xtable)
freq_table4 <- as.data.frame(freq_table4)
freq_table4_5 <- freq_table4[, c(1:5)]
freq_table4_3 <- freq_table4[, c(6:8)]
print(xtable(freq_table4_5 , caption='Four-state HMM'))
print(xtable(freq_table4_3 , caption='Four-state HMM'))




vPath_5state <- read.delim('output_5states/vPath_5states.ls')
freq_table5 <- apply(vPath_5state, 2, function(x)table(x))
freq_table5 <- as.data.frame(freq_table5)
freq_table5_5 <- freq_table5[, c(1:5)]
freq_table5_10 <- freq_table5[, c(6:10)]
freq_table5_12 <- freq_table5[, c(11:12)]
print(xtable(freq_table5_5 , caption='Five-state HMM'))
print(xtable(freq_table5_10 , caption='Five-state HMM'))
print(xtable(freq_table5_12 , caption='Five-state HMM'))




print(xtable(freq_table4_5))
print(xtable(freq_table5_5))



sessionInfo()



library(knitr)
knit("Analysis_StateModel.Rnw" ) # compile to tex
purl("Analysis_StateModel.Rnw", documentation = 0) # extract R code only
knit2pdf("Analysis_StateModel.Rnw")


