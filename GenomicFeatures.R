### R code from vignette source 'GenomicFeatures.Rnw'

###################################################
### code chunk number 1: loadGenomicFeatures
###################################################
library("GenomicFeatures")


###################################################
### code chunk number 2: supportedUCSCtables
###################################################
supportedUCSCtables()[1:4, ]


###################################################
### code chunk number 3: makeTranscriptDbFromUCSC (eval = FALSE)
###################################################
## mm9KG <- makeTranscriptDbFromUCSC(genome = "mm9", tablename = "knownGene")


###################################################
### code chunk number 4: discoverChromNames
###################################################
head(getChromInfoFromUCSC("hg19"))


###################################################
### code chunk number 5: makeTranscriptDbFromBiomart (eval = FALSE)
###################################################
## mmusculusEnsembl <-
##    makeTranscriptDbFromBiomart(biomart = "ensembl",
##                               dataset = "mmusculus_gene_ensembl")


###################################################
### code chunk number 6: saveFeatures 1 (eval = FALSE)
###################################################
## saveFeatures(mm9KG, file="fileName.sqlite")


###################################################
### code chunk number 7: loadFeatures-1 (eval = FALSE)
###################################################
## mm9KG <- loadFeatures("fileName.sqlite")


###################################################
### code chunk number 8: transcripts0
###################################################
samplefile <- system.file("extdata", "UCSC_knownGene_sample.sqlite",
                          package="GenomicFeatures")
txdb <- loadFeatures(samplefile)
txdb


###################################################
### code chunk number 9: transcripts1
###################################################
GR <- transcripts(txdb)
GR[1:3]


###################################################
### code chunk number 10: transcripts2
###################################################
GR <- transcripts(txdb, vals <- list(tx_chrom = "chr1", tx_strand = "+"))
length(GR)
unique(strand(GR))


###################################################
### code chunk number 11: transcriptsBy
###################################################
GRList <- transcriptsBy(txdb, by = "gene")
length(GRList)
names(GRList)[10:13]
GRList[11:12]


###################################################
### code chunk number 12: exonsBy
###################################################
GRList <- exonsBy(txdb, by = "tx")
length(GRList)
names(GRList)[10:13]
GRList[[12]]


###################################################
### code chunk number 13: internalID
###################################################
tx_ids <- names(GRList)
vals <- list(tx_id=tx_ids)
txs <- transcripts(txdb, vals, columns = c("tx_id", "tx_name"))
head(values(txs))


###################################################
### code chunk number 14: introns-UTRs
###################################################
length(intronsByTranscript(txdb))
length(fiveUTRsByTranscript(txdb))
length(threeUTRsByTranscript(txdb))


###################################################
### code chunk number 15: RNASEQData
###################################################
gr <- GRanges(
    seqnames = rep("chr5",4),
    ranges = IRanges(start = c(244620, 244670, 245804, 247502),
                     end = c(244652, 244702, 245836, 247534)),
    strand = rep("+", 4))


###################################################
### code chunk number 16: transcriptsByOverlaps
###################################################
transcriptsByOverlaps(txdb, gr)


###################################################
### code chunk number 17: exonsGroupedByTx
###################################################
annotGr <- exonsBy(txdb, "tx")


###################################################
### code chunk number 18: findOverlaps
###################################################
OL <- findOverlaps(query = annotGr, subject = gr)


###################################################
### code chunk number 19: getAnswers
###################################################
tdata <- annotGr[unique(queryHits(OL)),]
tdata
length(tdata)


###################################################
### code chunk number 20: SessionInfo
###################################################
sessionInfo()


