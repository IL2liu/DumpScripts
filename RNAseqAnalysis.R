

# set global chunk options
# for figures
opts_chunk$set(fig.path='figs/', fig.align='center', fig.show='hold',dev='CairoPDF', out.width='.4\\linewidth')
# replacing "=" into "->" to make it R thing
options(replace.assign=TRUE,width=90)
# caching chunks
opts_chunk$set(cache.extra = R.version,cache.path='cache/')
opts_chunk$set(cache.extra = rand_seed)



## 
## tss_rep1 <- read.table("input_RNAseq/wgEncodeCaltechRnaSeqH1hescR2x75Il200TSSRep1V3.gtf.gz")
## tss_rep2 <- read.table("input_RNAseq/wgEncodeCaltechRnaSeqH1hescR2x75Il200TSSRep2V3.gtf.gz")
## tss_rep3 <- read.table("input_RNAseq/wgEncodeCaltechRnaSeqH1hescR2x75Il200TSSRep3V3.gtf.gz")
## tss_rep4 <- read.table("input_RNAseq/wgEncodeCaltechRnaSeqH1hescR2x75Il200TSSRep4V3.gtf.gz")
## 
## head(tss_rep1)
## head(tss_rep2)
## head(tss_rep3)
## head(tss_rep4)
## 
## dim(tss_rep1)
## dim(tss_rep2)
## dim(tss_rep3)
## dim(tss_rep4)



## library(plyr)
## # subset each replicate for [chr, start, gene_id, transcript_id, FPKM]
## tss_rep1 <- tss_rep1[, c(1, 4, 10, 13, 16, 25)]
## tss_rep2 <- tss_rep2[, c(1, 4, 10, 13, 16, 25)]
## tss_rep3 <- tss_rep3[, c(1, 4, 10, 13, 16, 25)]
## tss_rep4 <- tss_rep4[, c(1, 4, 10, 13, 16, 25)]
## 
## colnames(tss_rep1) <- c('chr', 'start', 'geneID', 'transcriptID', 'FPKM1', 'npIDR1')
## colnames(tss_rep2) <- c('chr', 'start', 'geneID', 'transcriptID', 'FPKM2', 'npIDR2')
## colnames(tss_rep3) <- c('chr', 'start', 'geneID', 'transcriptID', 'FPKM3', 'npIDR3')
## colnames(tss_rep4) <- c('chr', 'start', 'geneID', 'transcriptID', 'FPKM4', 'npIDR4')
## 
## head(tss_rep1)
## head(tss_rep2)
## head(tss_rep3)
## head(tss_rep4)
## 
## tss_reps <- Reduce(function(...) merge(...,by=c('chr', 'start', 'geneID', 'transcriptID'), all=TRUE, suffixes=c(1,2,3,4)), list(tss_rep1, tss_rep2, tss_rep3, tss_rep4))
## 
## # obtained FPKM from 4 replicates and npIDR from rep1 and rep3
## # note npIDR1 and npIDR3 are computed as IDR of rep1-rep2 and IDR or rep3-rep4 respectively.
## tss_reps <- tss_reps[, c(1,2,3,4,5,7,9,11,6,10)]
## 
## dim(tss_reps)
## head(tss_reps)



## 
## tss_bed <- tss_reps[, c(1,2)]
## tss_bed$end <- tss_bed$start + 1
## head(tss_bed)
## dim(tss_bed)
## 
## # all the common TSS in these four replicates
## write.table(tss_bed, "output/RNAseq_AllTSS.hg19.bed",
##             quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
## 



## expressed_tss <- tss_reps[tss_reps$npIDR1 & tss_reps$npIDR3 < 0.1,]
## NONexpressed_tss <- tss_reps[tss_reps$npIDR1 & tss_reps$npIDR3 > 0.1,]
## expressed_tss <- expressed_tss[complete.cases(expressed_tss), ]
## NONexpressed_tss <- NONexpressed_tss[complete.cases(NONexpressed_tss), ]
## 
## expressedTSS_bed <- expressed_tss[,c(1,2)]
## NONexpressedTSS_bed <- NONexpressed_tss[,c(1,2)]
## expressedTSS_bed$end <- expressedTSS_bed$start + 1
## NONexpressedTSS_bed$end <- NONexpressedTSS_bed$start + 1
## 
## # expressed and non-expressed common TSS in these four replicates
## write.table(expressedTSS_bed, "output/expressedRNAseq.hg19.bed",
##             quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE)
## write.table(NONexpressedTSS_bed, "output/NONexpressedRNAseq.hg19.bed",
##             quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE)
## 



## allTSS <- read.delim('output_TSScolor/RNAseq_AllTSS.h1.color.bed', header=FALSE)
## expressedTSS <- read.delim('output_TSScolor/expressedRNAseq.h1.color.bed', header=FALSE)
## NONexpressedTSS <- read.delim('output_TSScolor/NONexpressedRNAseq.h1.color.bed', header=FALSE)
## head(allTSS)
## 
## getCoverageGenomicFeature <- function(x) {
##       coverage <- tapply(INDEX = x$V4, X = x$V3 - x$V2, sum)
##       percentage <- coverage/sum(coverage)
##       return(percentage)
## }
## 
## coverageALLTSS <- getCoverageGenomicFeature(allTSS)
## coverageExpressedTSS <- getCoverageGenomicFeature(expressedTSS)
## coverageNONexpressedTSS <- getCoverageGenomicFeature(NONexpressedTSS)
## 
## # combine into one dataframe
## library(ggplot2)
## library(reshape)
## H1_tss <- data.frame(
##             All = coverageALLTSS,
##             Expressed = coverageExpressedTSS,
##             NonExpressed = coverageNONexpressedTSS
##             )
## 
## head(H1_tss)
## 
## H1_tss$category <- rownames(H1_tss)
## mH1_tss <- melt(H1_tss, id.vars="category")
## head(mH1_tss)
## 
## plot_npIDR <- ggplot(data=mH1_tss, aes(x=variable, y=value, fill=factor(category)))+
##       geom_bar(colour="black", stat='identity') +
##       scale_fill_manual(name="Color States",
##                         values=c("red","deeppink2", "gold2", "black")) +
##       xlab('') +
##       ylab('Coverage')+
##       theme(text = element_text(size=20),
##                   axis.text.y = element_text(angle=0, vjust=1, colour="black")) +
##       theme(legend.position = "none")
## 
## ggsave(filename='figs/(Non)ExpressedTSSnpIDR.pdf', plot=plot_npIDR)
## 



## tss_FPKM <- tss_reps[, c(1,2,5,6,7,8)]
## head(tss_FPKM)
## tss_FPKM$Log.mean <- log(rowMeans(tss_FPKM[,c(3,4, 5,6)])+1)
## tss_FPKM <- tss_FPKM[, c(1,2,7)]
## 
## colnames(allTSS) <- c("chr", "start", "end", "color")
## head(allTSS)
## 
## allTSS_FPKM <- merge(allTSS, tss_FPKM, by = c("chr", "start"))
## head(allTSS_FPKM)
## 
## library(ggplot2)
## plot_FPKM <- qplot(allTSS_FPKM$Log.mean, binwidth = 0.5,
##                   fill = factor(allTSS_FPKM$color))+ xlab('Log FPKM')+
##                   scale_fill_manual(name = "Color States", values = c("red", "deeppink2", "gold2", "black")
##                    , labels = c("Red", "Pink", "Yellow", "black"))
## ggsave(filename='figs/(Non)ExpressedTSSFPKM.pdf', plot=plot_FPKM)



## GO_IDs <- merge(allTSS, tss_reps, by=c("chr", "start"))
## 
## geneID <- GO_IDs[, c("geneID", "color")]
## transcriptID <- GO_IDs[, c("geneID", "color")]
## 
## geneID <- within(geneID, color <- factor(color, labels = c("Red","Pink", "Yellow", "Black")))
## transcriptID <- within(transcriptID, color <- factor(color, labels = c("Red","Pink", "Yellow", "Black")))
## 
## write.table(geneID, 'output_GOid/geneID', quote=FALSE, row.names=F, col.names=F)
## write.table(transcriptID, 'output_GOid/transcriptID', quote=FALSE, row.names=F, col.names=F)
## 



sessionInfo()



library(knitr)
knit("RNAseqAnalysis.Rnw" ) # compile to tex
purl("RNAseqAnalysis.Rnw", documentation = 0) # extract R code only
knit2pdf("RNAseqAnalysis.Rnw")


