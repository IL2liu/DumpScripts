library(Biobase)
library(GEOquery)
library(plyr)

setwd("/TEMP_DDN/users/gfilion/rlim/Desktop/Projects/CHIPdatasets/HumanAlignmentData/fasq-cells:H1-hESC/Ranalysis/EnrichmentFeaturesOnChromatinColor/H1_RNAseq/H1_RNAseq")

tss_rep1 <- read.table("input/wgEncodeCaltechRnaSeqH1hescR2x75Il200TSSRep1V3.gtf.gz")
tss_rep2 <- read.table("input/wgEncodeCaltechRnaSeqH1hescR2x75Il200TSSRep2V3.gtf.gz")
tss_rep3 <- read.table("input/wgEncodeCaltechRnaSeqH1hescR2x75Il200TSSRep3V3.gtf.gz")
tss_rep4 <- read.table("input/wgEncodeCaltechRnaSeqH1hescR2x75Il200TSSRep4V3.gtf.gz")

head(tss_rep1)
head(tss_rep2)
head(tss_rep3)
head(tss_rep4)
  
dim(tss_rep1)
dim(tss_rep2)
dim(tss_rep3)
dim(tss_rep4)

rep1_tss <- tss_rep1[,3:5]
rep2_tss <- tss_rep2[,3:5]
rep3_tss <- tss_rep3[,3:5]
rep4_tss <- tss_rep4[,3:5]
head(rep1_tss)
head(rep2_tss)
head(rep3_tss)
head(rep4_tss)

identical(rep1_tss,rep2_tss)
identical(rep4_tss, rep3_tss)
# replicates 1 and 2 are indetical & replicates 3 and 4 are identical in their Tx


# npIDR < 0.1 is expressed
# ref: http://authors.library.caltech.edu/35106/2/nature11233-s1.pdf
identical(tss_rep1$V25, tss_rep2$V25) 
identical(tss_rep3$V25, tss_rep4$V25)
#npIDRs rep1 and rep2 are identical, use either one, we used rep1 for subsetting
tss_expressed <- tss_rep1[tss_rep1$V25<0.1, c(1,4,5)]
tss_nonExpressed <- tss_rep1[tss_rep1$V25>0.1, c(1,4,5)]
tss_all <- tss_rep1[,c(1,4,5)]

head(tss_expressed)
head(tss_nonExpressed)
head(tss_all)

dim(tss_expressed)
dim(tss_nonExpressed)
dim(tss_all)

# to have the end (after 1 nucleotide Tx start site)
tss_expressed$V5 <- tss_expressed$V5 + 1
tss_nonExpressed$V5 <- tss_nonExpressed$V5+1

# output
write.table(tss_expressed, "output/RNAseqTSSH1ExpressedGene.hg19.bed", 
            quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

write.table(tss_nonExpressed, "output/RNAseqTSSH1NonExpressedGene.hg19.bed",
            quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

write.table(tss_all, "output/RNAseqTSSAllH1.hg19.bed",
            quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

