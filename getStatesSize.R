

# set global chunk options
# for figures
opts_chunk$set(fig.path='figs/', fig.align='center', fig.show='hold',dev='CairoPDF', out.width='.4\\linewidth')
# replacing "=" into "->" to make it R thing
options(replace.assign=TRUE,width=90)
# caching chunks
opts_chunk$set(cache.extra = R.version,cache.path='cache/')
opts_chunk$set(cache.extra = rand_seed)



color_domains <- read.table("input_color/domain_states_hESC.txt",
                            header=TRUE)
head(color_domains)
colnames(color_domains) <- c("chr", "start", "end", "tag")
color_domains$size <- color_domains$end - color_domains$start
color_size <- color_domains[,c(4,5)]
head(color_size)
colnames(color_size) <- c("state", "size")
color_size <- within(color_size, state <- factor(state,
                                                 labels=c("Red", "Pink",
                                                          "Yellow", "Black")))
color_size$state <- factor(color_size$state, 
                           levels=c("Red", "Pink", 
                                    "Yellow", "Black"))

library(ggplot2)

pdf("figs/domain_size.pdf")
DomainsizePlot <- ggplot(color_size, aes(state, size, fill=state)) + 
                  geom_boxplot(outlier.shape = NA, fatten=5) + 
                  coord_cartesian(ylim = c(0,50000))+ 
                  scale_fill_manual(name = "Color States", 
                  values=c("red", "deeppink2", "gold2","black" )) 
DomainsizePlot
dev.off()




sessionInfo()



library(knitr)
knit("getStatesSize.Rnw" ) # compile to tex
purl("getStatesSize.Rnw", documentation = 0) # extract R code only
knit2pdf("getStatesSize.Rnw")


