

# set global chunk options
# for figures
opts_chunk$set(fig.path='figs/', fig.align='center', fig.show='hold',dev='CairoPDF', out.width='.4\\linewidth')
# replacing "=" into "->" to make it R thing
options(replace.assign=TRUE,width=90)
# caching chunks
opts_chunk$set(cache.extra = R.version,cache.path='cache/')
opts_chunk$set(cache.extra = rand_seed)



# Create 3D pie chart for color chromatin in human
H1color_states <- read.table("input_color/domain_states_hESC.txt", header=T)
head(H1color_states)
colnames(H1color_states) <- c("chr", "start", "end", "tag")
H1color_states <- within(H1color_states, tag <- factor(tag, 
                                     labels=c("Red", "Pink", "Yellow", "Black")))

coverage_color <- tapply(INDEX=H1color_states$tag, 
                         X=as.numeric(H1color_states$end)-as.numeric(H1color_states$start),
                         sum)
percent_coverage_color <- coverage_color/sum(as.numeric(coverage_color))*100
percent_coverage_color

# reorder the color states
percent_coverage_color
library(plotrix)
color_states =c("red", "deeppink2", "gold2", "black")
pdf("figs/3dPizza.pdf")
pie3D(percent_coverage_color,col = color_states) 
dev.off()



sessionInfo()



library(knitr)
knit("get3Dpizza.Rnw" ) # compile to tex
purl("get3Dpizza.Rnw", documentation = 0) # extract R code only
knit2pdf("get3Dpizza.Rnw")


