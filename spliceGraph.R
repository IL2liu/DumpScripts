### R code from vignette source 'spliceGraph.Rnw'

###################################################
### code chunk number 1: spliceGraph.Rnw:32-33
###################################################
reCnt <- FALSE


###################################################
### code chunk number 2: loadTxdb
###################################################
library("TxDb.Hsapiens.UCSC.hg18.knownGene")
txdb <- TxDb.Hsapiens.UCSC.hg18.knownGene


###################################################
### code chunk number 3: deactChr
###################################################
activeChr <- ! isActiveSeq(txdb)
activeChr[names(activeChr) == "chr19"] <- TRUE
isActiveSeq(txdb) <- activeChr


###################################################
### code chunk number 4: loadGenomicFeatures
###################################################
library("GenomicFeatures")
exsByEdges <- spliceGraph(txdb)


###################################################
### code chunk number 5: spliceGraph.Rnw:218-219
###################################################
exsByEdges[1]


###################################################
### code chunk number 6: origExToDJEx
###################################################
exsByEdges.flat <- unlist(exsByEdges, use.names=FALSE)
orig.Ex <- values(exsByEdges.flat)[["exon_ids"]]
disJ.Ex <- values(exsByEdges.flat)[["disJ_exon_id"]]

newExNames <- rep(disJ.Ex, elementLengths(orig.Ex))
origExToNewEx <- data.frame(orig.Ex=unlist(orig.Ex), 
                            disJ.Ex=newExNames, 
                            row.names=NULL)

head(origExToNewEx)


###################################################
### code chunk number 7: getGeneToEdge
###################################################
gnIDs <- values(exsByEdges)[["gene_id"]]
edgeIDs <- names(exsByEdges)
gnToEdge <- data.frame(gnIDs, edgeIDs, row.names=edgeIDs)

head(gnToEdge)


###################################################
### code chunk number 8: countingReads
###################################################
if(reCnt) {
  require("Rsamtools")
  bamPath <- "/shared/labs/EDI/users/mfitzgib/Solexa"
  fls <- sub(".bai$", "", 
             list.files(bamPath, recursive=TRUE,
                        pattern="accepted_hits.*bai$", 
                        full=TRUE))
  fls <- fls[c(1:3, 62:64)]
  bfs <- BamFileList(fls)
  names(bfs) <- 
    gsub("/shared/labs/EDI/users/mfitzgib/Solexa/tophat_",
         "",fls)
}


###################################################
### code chunk number 9: spliceGraph.Rnw:311-326
###################################################
if(reCnt) {
  library(parallel)
  resEx.byEdge <- 
    summarizeOverlaps(features = exsByEdges, reads = bfs,
                      mc.cores = getOption("mc.cores", 3L))
  
  cD.exByEdge <- assays(resEx.byEdge)$counts
  colnames(cD.exByEdge) <- 
    sub("/accepted_hits.bam", "", colnames(cD.exByEdge))
  save(cD.exByEdge, file="cD.exByEdge-SG-Vig.Rda")
} else {
  fn <- system.file("extdata", "cD.exByEdge-SG-Vig.Rda", 
                    package="GenomicFeatures")
  load(fn)
}


###################################################
### code chunk number 10: spliceGraph.Rnw:341-351
###################################################
exsByGenes <- exonsBy(txdb, "gene")
gnNames <- rep(names(exsByGenes), elementLengths(exsByGenes))
exsByGenes.flat <- unlist(exsByGenes, use.names=FALSE)

# remove duplicates
notDupl <- !duplicated(values(exsByGenes.flat)[["exon_id"]])
exsByGenes.flat <- exsByGenes.flat[notDupl]
names(exsByGenes.flat) <- values(exsByGenes.flat)[["exon_id"]]
gnNames <- gnNames[notDupl]
names(gnNames) <- names(exsByGenes.flat)


###################################################
### code chunk number 11: spliceGraph.Rnw:360-374
###################################################
if(reCnt) {
  res.exsByGenes <- 
    summarizeOverlaps(features = exsByGenes.flat, reads = bfs,
                      mc.cores = getOption("mc.cores", 3L))

  cD.exsByGenes <- assays(res.exsByGenes)$counts
  colnames(cD.exsByGenes) <- 
    sub("/accepted_hits.bam", "", colnames(cD.exsByGenes))
  save(cD.exsByGenes, file="cD.exsByGenes-SG-Vig.Rda")
} else {
  fn <- system.file("extdata", "cD.exsByGenes-SG-Vig.Rda", 
                    package="GenomicFeatures")
  load(fn)
}


###################################################
### code chunk number 12: spliceGraph.Rnw:399-410
###################################################
edIds <- names(exsByEdges)
gnIds.edges <- gnToEdge[rownames(cD.exByEdge),]$gnIDs
nrOfEdgesPerGns <- 
  elementLengths(split(edIds, factor(gnIds.edges)))
sglEdgeGns <- nrOfEdgesPerGns == 1

exIds <- names(exsByGenes.flat)
gnIds.exs <- gnNames[names(exsByGenes.flat)]
nrOfExsPerGns <- 
  elementLengths(split(exIds, factor(gnIds.exs)))
sglExGns <- nrOfExsPerGns == 1


###################################################
### code chunk number 13: spliceGraph.Rnw:420-421
###################################################
cD.exByEdge <- cD.exByEdge[, colnames(cD.exsByGenes)]


###################################################
### code chunk number 14: spliceGraph.Rnw:431-437
###################################################
benign <- grepl("Benign",colnames(cD.exByEdge ) )
design <- factor(ifelse(benign, "Benign", "SOC"), 
                 levels=c("Benign", "SOC"))
names(design) <- ifelse(benign, "Benign", "SOC")

design


###################################################
### code chunk number 15: spliceGraph.Rnw:450-472
###################################################
exsByEdges.flat <- unlist(exsByEdges, use.names=FALSE)
edgeIDs <- rep(names(exsByEdges), elementLengths(exsByEdges))
rle <- strand(exsByEdges.flat)
start <- start(exsByEdges.flat)
end <- end(exsByEdges.flat)

temp <- data.frame(start, end, edgeIDs, 
                   strand=rep(rle@values, rle@lengths),
                   chr=rep(seqnames(exsByEdges.flat)@values, 
                     seqnames(exsByEdges.flat)@lengths))


temp <- temp[order(temp$edgeIDs, -temp$end), ]
end <- temp$end[!duplicated(temp$edgeIDs)]
names(end) <- temp$edgeIDs[!duplicated(temp$edgeIDs)]

temp <- temp[order(temp$edgeIDs, temp$start), ]
T <- temp[!duplicated(temp$edgeIDs),]
rownames(T) <- T$edgeIDs
T$end <- end[rownames(T)]

T <- T[order(gnToEdge[T$edgeIDs,]$gnIDs, T$start), ]


###################################################
### code chunk number 16: spliceGraph.Rnw:481-492
###################################################
library(DEXSeq)

eData.Edge <- 
  newExonCountSet(countData = cD.exByEdge[T$edgeIDs,] ,
                  design = design,
                  geneIDs = gnToEdge[T$edgeIDs,]$gnIDs,
                  exonIDs = T$edgeIDs,
                  exonIntervals = T)
eData.Edge <- 
  eData.Edge[! gnToEdge[T$edgeIDs,]$gnIDs %in% 
             names(sglEdgeGns[sglEdgeGns]), ]


###################################################
### code chunk number 17: spliceGraph.Rnw:498-505
###################################################
T <- data.frame(chr=rep(seqnames(exsByGenes.flat)@values, 
                  seqnames(exsByGenes.flat)@lengths), 
                start=start(exsByGenes.flat), 
                end=end(exsByGenes.flat), 
                strand=rep(strand(exsByGenes.flat)@values, 
                  strand(exsByGenes.flat)@lengths),
                row.names=values(exsByGenes.flat)[["exon_id"]])


###################################################
### code chunk number 18: spliceGraph.Rnw:514-523
###################################################
eData.Ex <- 
  newExonCountSet(countData = cD.exsByGenes,
                  design = design,
                  geneIDs = gnNames[rownames(cD.exsByGenes)],
                  exonIDs = rownames(cD.exsByGenes),
                  exonIntervals=T[rownames(cD.exsByGenes),])
eData.Ex <- 
  eData.Ex[! gnNames[rownames(cD.exsByGenes)] %in% 
           names(sglExGns[sglExGns]), ]


###################################################
### code chunk number 19: spliceGraph.Rnw:536-542
###################################################
eData.Ex <- estimateSizeFactors(eData.Ex)
eData.Ex <- estimateDispersions(eData.Ex)
eData.Ex <- fitDispersionFunction(eData.Ex)
eData.Ex <- testForDEU(eData.Ex)
eData.Ex <- estimatelog2FoldChanges(eData.Ex)
tt.Ex <- DEUresultTable(eData.Ex)


###################################################
### code chunk number 20: spliceGraph.Rnw:547-553
###################################################
eData.Edge <- estimateSizeFactors(eData.Edge)
eData.Edge <- estimateDispersions(eData.Edge)
eData.Edge <- fitDispersionFunction(eData.Edge)
eData.Edge <- testForDEU(eData.Edge)
eData.Edge <- estimatelog2FoldChanges(eData.Edge)
tt.Edge <- DEUresultTable(eData.Edge)


###################################################
### code chunk number 21: spliceGraph.Rnw:560-577
###################################################
plotDisp <- function(eData, case="") {
  meanvalues <- rowMeans(counts(eData))
  plot(meanvalues, fData(eData)$dispBeforeSharing, 
       main = paste(case, " mean vs CR dispersion", sep=":"), 
       frame.plot=FALSE, pch=16, cex=0.8, log = "xy",
       ylab="Dispersion", xlab="Mean counts")
  x <- 0.01:max(meanvalues)
  y <- eData@dispFitCoefs[1] + eData@dispFitCoefs[2]/x
  lines(x, y, col = "purple", lwd=2)
}

fn <- "Disp-plots.png"
png(fn,  width=11, height=5.5, res=300, units="in")
par(mfrow=c(1,2))
plotDisp(eData.Ex, "Exons")
plotDisp(eData.Edge, "Edges")
dev.off()


###################################################
### code chunk number 22: spliceGraph.Rnw:598-606
###################################################
library(DESeq)
cds.Ex <- newCountDataSet(countData=cD.exsByGenes, 
                          conditions=design)
cds.Ex <- cds.Ex[sglExGns,]

cds.Edge <- newCountDataSet(countData=cD.exByEdge[sglEdgeGns,], 
                            conditions=design)
cds.Edge <- cds.Edge[sglExGns,]


###################################################
### code chunk number 23: spliceGraph.Rnw:615-618
###################################################
cds.Ex <- estimateSizeFactors( cds.Ex )
cds.Ex <- estimateDispersions( cds.Ex )
ttSingle.Ex <- nbinomTest( cds.Ex,"SOC", "Benign" )


###################################################
### code chunk number 24: spliceGraph.Rnw:623-626
###################################################
cds.Edge <- estimateSizeFactors( cds.Edge )
cds.Edge <- estimateDispersions( cds.Edge )
ttSingle.Edge <- nbinomTest( cds.Edge,"SOC", "Benign" )


###################################################
### code chunk number 25: spliceGraph.Rnw:634-648
###################################################
int.DEXSeq <- c("exonID", "pvalue", "padjust", 
                "meanBase", "log2fold(SOC/Benign)")
int.DESeq <- c("id", "pval", "padj", "baseMean", 
               "log2FoldChange")

temp <- tt.Edge[, int.DEXSeq]
colnames(temp) <- int.DESeq
finTT.Edge <- rbind(ttSingle.Edge[, int.DESeq], temp)
finTT.Edge$gnID <- gnToEdge[finTT.Edge$id,]$gnIDs

temp <- tt.Ex[, int.DEXSeq]
colnames(temp) <- int.DESeq
finTT.Ex <- rbind(ttSingle.Ex[, int.DESeq], temp)
finTT.Ex$gnID <- gnNames[finTT.Ex$id]


###################################################
### code chunk number 26: spliceGraph.Rnw:671-674
###################################################
sum(cD.exByEdge)
sum(cD.exsByGenes)
sum(cD.exsByGenes*100)/sum(cD.exByEdge)


###################################################
### code chunk number 27: spliceGraph.Rnw:694-706
###################################################
fn <- "MA-plots.png"
png(fn,  width=11, height=5.5, res=300, units="in")
par(mfrow=c(1,2))
plot(finTT.Ex$baseMean, finTT.Ex$log2FoldChange, log = "x",
     col = ifelse(finTT.Ex$padj < 0.1, "purple", "black"), 
     ylim = c(-5, 5), main = "Exons MvsA", pch=16, cex=0.7, 
     frame.plot=FALSE)
plot(finTT.Edge$baseMean, finTT.Edge$log2FoldChange,
     col = ifelse(finTT.Edge$padj < 0.1, "purple", "black"), 
     ylim = c(-5, 5), main = "Edges MvsA", pch=16, cex=0.7, 
     log = "x", frame.plot=FALSE)
dev.off()


###################################################
### code chunk number 28: spliceGraph.Rnw:723-731
###################################################
ex.test <- nrow(finTT.Ex)
edge.test <- nrow(finTT.Edge)

sig.diff.ex <- sum(finTT.Ex$padj < 0.05, na.rm=TRUE)
sig.diff.edge <- sum(finTT.Edge$padj < 0.05, na.rm=TRUE)

sig.diff.ex*100/ex.test
sig.diff.edge*100/edge.test


###################################################
### code chunk number 29: spliceGraph.Rnw:739-764
###################################################
p.cuts <- 10^seq( 0, -6, length.out=100 )

ps.Ex <- sapply(p.cuts, function(p.cut) {
  sum(finTT.Ex$pval < p.cut, na.rm=TRUE)*100/
    length(finTT.Ex$pval)
})

ps.Edge <- sapply(p.cuts, function(p.cut) {
  sum(finTT.Edge$pval < p.cut, na.rm=TRUE)*100/
    length(finTT.Edge$pval)
})


fn <- "Comparison-edges-exons.png"
png(fn,  width=7, height=7, res=300, units="in")
plot(ps.Edge, p.cuts, type="l", col="purple",
     xlab="% of differentially expressed elements",
     ylab="P-value cut off", lwd=2, frame.plot=FALSE, 
     log="xy", main="Exons versus edges")
lines(ps.Ex, p.cuts, lty=2, lwd=2)
grid(lwd=2)
legend("bottomright", legend=c("Edges", "Exons"), 
       col=c("purple", "black"), lty=c(1,2), lwd=2,
       border=NA, box.col=NA, bg=NA)
dev.off()


###################################################
### code chunk number 30: spliceGraph.Rnw:790-794
###################################################
library(xtable)
T <- xtable(head(finTT.Ex[order(finTT.Ex$padj),], n=10), 
            caption="Top table of the exon model", label="table:ex")
print(T, include.rownames=FALSE, table.placement="H")


###################################################
### code chunk number 31: spliceGraph.Rnw:797-801
###################################################
library(xtable)
T <- xtable(head(finTT.Edge[order(finTT.Edge$padj),], n=10), 
            caption="Top table of the edge model", label="table:edge")
print(T, include.rownames=FALSE, table.placement="H")


###################################################
### code chunk number 32: spliceGraph.Rnw:808-827
###################################################
plotExp <- function(fn, eDat) {
  png(fn,  width=11, height=5.5, res=300, units="in")
  par(mfrow=c(1,2))
  COL <- c("#3399FF", "#FF3333")
  plotDEXSeq(eDat, "90522", cex.axis = 1.2, cex = 1.3,
             lwd = 2, legend = TRUE, color=COL, 
             color.samples=COL[design], splicing=TRUE)
  plotDEXSeq(eDat, "90522", expression = FALSE, 
             norCounts = TRUE, cex.axis = 1.2, cex = 1.3, 
             lwd = 2, legend = TRUE,
             color=COL, 
             color.samples=COL[design], splicing=TRUE)
  dev.off()
}

fn <- "Edge-plot.png"
plotExp(fn, eData.Edge)
fn <- "Ex-plot.png"
plotExp(fn, eData.Ex)


###################################################
### code chunk number 33: SessionInfo
###################################################
sessionInfo()


