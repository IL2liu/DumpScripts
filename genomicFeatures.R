

# set global chunk options
# for figures
opts_chunk$set(fig.path='figs/', fig.align='center', fig.show='hold',dev='CairoPDF', out.width='.4\\linewidth')
# replacing "=" into "->" to make it R thing
options(replace.assign=TRUE,width=90)
# caching chunks
opts_chunk$set(cache.extra = R.version,cache.path='cache/')
opts_chunk$set(cache.extra = rand_seed)



## CpG <- read.delim("output/CpGIsland.h1Color.bed", header=FALSE)
## Exon <- read.delim("output/RefSeqExon.h1Color.bed", header=FALSE)
## Gene <- read.delim("output/RefSeqGene.h1Color.bed", header=FALSE)
## TES <- read.delim("output/RefSeqTES.h1Color.bed", header=FALSE)
## TSS <- read.delim("output/RefSeqTSS.h1Color.bed", header=FALSE)
## TSS2kb <- read.delim("output/RefSeqTSS2kb.h1Color.bed", header=FALSE)
## OpenChromatin <- read.delim("output/ValidateOpenChromatin.h1Color.bed", header=FALSE)
## laminB1 <- read.delim("output/laminB1lads.h1Color.bed", header=FALSE)




getCoverageGenomicFeature <- function(x){
  coverage <- tapply(INDEX=x$V4, X=x$V3-x$V2, sum)
  percentage <- coverage/sum(coverage)
  return(percentage)
}

EnrichCpG <- getCoverageGenomicFeature(CpG)
EnrichExon <- getCoverageGenomicFeature(Exon)
EnrichGene <- getCoverageGenomicFeature(Gene)
EnrichTES <- getCoverageGenomicFeature(TES)
EnrichTSS <- getCoverageGenomicFeature(TSS)
EnrichTSS2kb <- getCoverageGenomicFeature(TSS2kb)
EnrichLamin <- getCoverageGenomicFeature(laminB1)
EnrichOpenChromatin <- getCoverageGenomicFeature(OpenChromatin)



pdf("figs/FeatureCoverageH1.pdf", height=6, width=6)
par(mfrow=c(2,4))
color_states =c("red", "deeppink2", "gold2", "black")
pie(EnrichCpG, col=color_states, main="CpG", labels="")
pie(EnrichExon, col=color_states, main="Exon", labels="" )
pie(EnrichGene, col=color_states, main="Gene", labels="")
pie(EnrichTES, col=color_states, main="TES", labels="")
pie(EnrichTSS, col=color_states, main="TSS", labels="")
pie(EnrichTSS2kb, col=color_states, main="TSS2kb", labels="")
pie(EnrichLamin, col=color_states, main="LaminB1", labels="")
pie(EnrichOpenChromatin, col=color_states, main="Open Chromatin", labels="")
dev.off()



library(plotrix)
pdf("figs/3DFeatureCoverageH1.pdf", height=6, width=6)
par(mfrow=c(2,4))
color_states =c("red", "deeppink2", "gold2", "black")
pie3D(EnrichCpG, col=color_states, main="CpG", labels="",explode=0.1)
pie3D(EnrichExon, col=color_states, main="Exon", labels="",explode=0.1)
pie3D(EnrichGene, col=color_states, main="Gene", labels="",explode=0.1)
pie3D(EnrichTES, col=color_states, main="TES", labels="",explode=0.1)
pie3D(EnrichTSS, col=color_states, main="TSS", labels="",explode=0.1)
pie3D(EnrichTSS2kb, col=color_states, main="TSS2kb", labels="",explode=0.1)
pie3D(EnrichLamin, col=color_states, main="Lamin B1", labels="",explode=0.1)
pie3D(EnrichOpenChromatin, col=color_states, main="Open Chromatin", labels="",explode=0.1)
dev.off()




library(ggplot2)
library(reshape)
H1_allGenomicFeatures <- data.frame(
  
      OpenChromatin = EnrichOpenChromatin,
      #TSS = EnrichTSS,
      #TSS2kb = EnrichTSS2kb,
      #TES = EnrichTES,
      Exon = EnrichExon,
      Gene = EnrichGene,
      CpG = EnrichCpG,
      Lamin = EnrichLamin
)

head(H1_allGenomicFeatures)
H1_allGenomicFeatures$category <- row.names(H1_allGenomicFeatures)
mH1_allGenomicFeatures <- melt(H1_allGenomicFeatures, id.vars="category")
head(mH1_allGenomicFeatures)
dim(mH1_allGenomicFeatures)
factor()

# attempt to reorganize the order of color states

#head(mH1_allGenomicFeatures)
#myfunkyvector = c(1,3,5,2,4) + 5*rep(0:7, each=5)
#mH1_allGenomicFeatures[myfunkyvector,]
#mH1_allGenomicFeatures <- mH1_allGenomicFeatures[myfunkyvector,]
#rownames(mH1_allGenomicFeatures) <- NULL
#head(mH1_allGenomicFeatures)


# reorganize the order of color states
mH1_allGenomicFeatures$category <- factor(mH1_allGenomicFeatures$category, 
                                         levels=c(0,1,2,3))

p <- ggplot(data = mH1_allGenomicFeatures, 
            aes(x = variable, y = value, fill = factor(category)))+
     geom_bar(colour="black", stat='identity') +
     coord_flip() +
     scale_fill_manual(name = "Color States", values=c("red", "deeppink2", "gold2", "black")) +
     ylab('Coverage')+
     xlab('') + 
     theme(text = element_text(size=20),
           axis.text.y = element_text(angle=0, vjust=1, colour="black")) +
     theme(legend.position = "none") 

ggsave(filename='figs/AllFeatureCoverage.pdf', plot=p, height=6)




sessionInfo()



library(knitr)
knit("genomicFeatures.Rnw" ) # compile to tex
purl("genomicFeatures.Rnw", documentation = 0) # extract R code only
knit2pdf("genomicFeatures.Rnw")


