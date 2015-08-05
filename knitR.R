

# set global chunk options
opts_chunk$set(fig.path='figs/', fig.align='center', fig.show='hold',dev='CairoPDF', out.width='.4\\linewidth')
options(replace.assign=TRUE,width=90)
opts_chunk$set(cache.extra = R.version,cache.path='cache/')
opts_chunk$set(cache.extra = rand_seed)



set.seed(1121)
(x=rnorm(20))
mean(x);var(x)
source("foo.R")



## two plots side by side (option fig.show='hold')
par(mar=c(4,4,.1,.1),cex.lab=.95,cex.axis=.9,mgp=c(2,.7,0),tcl=-.3,las=1)
boxplot(x)
hist(x,main='')



library(ggbio)



plot(1)         # high-level plot
abline(0, 1)    # low-level change
## many low-level changes in a loop (a single R expression)
for(i in 1:10) {
    abline(v = i, lty = 2)
}



library(knitr)
knit("knitR.Rnw" ) # compile to tex
purl("knitR.Rnw", documentation = 0) # extract R code only
knit2pdf("knitR.Rnw")


