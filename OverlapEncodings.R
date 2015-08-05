### R code from vignette source 'OverlapEncodings.Rnw'

###################################################
### code chunk number 1: options
###################################################
options(width=100)


###################################################
### code chunk number 2: untreated1_chr4
###################################################
library(pasillaBamSubset)
untreated1_chr4()


###################################################
### code chunk number 3: readGappedAlignments
###################################################
library(GenomicRanges)
library(Rsamtools)
flag0 <- scanBamFlag(isDuplicate=FALSE, isValidVendorRead=TRUE)
param0 <- ScanBamParam(flag=flag0)
gal14 <- readGappedAlignments(untreated1_chr4(), use.names=TRUE, param=param0)


###################################################
### code chunk number 4: ngap
###################################################
head(gal14)
table(ngap(gal14))


###################################################
### code chunk number 5: makeTranscriptDbFromUCSC
###################################################
library(GenomicFeatures)
dm3_refGene_txdb <- makeTranscriptDbFromUCSC(genome="dm3", tablename="refGene")
exbytx <- exonsBy(dm3_refGene_txdb, by="tx")


###################################################
### code chunk number 6: check-for-trans-splicing
###################################################
table(elementLengths(runLength(seqnames(exbytx))))
table(elementLengths(runLength(strand(exbytx))))


###################################################
### code chunk number 7: exbytx_strand
###################################################
exbytx_strand <- unlist(runValue(strand(exbytx)), use.names=FALSE)


###################################################
### code chunk number 8: ov14
###################################################
ov14 <- findOverlaps(gal14, exbytx, ignore.strand=TRUE)


###################################################
### code chunk number 9: grl14of
###################################################
grl14o <- grglist(gal14, order.as.in.query=TRUE)
grl14f <- flipQuery(grl14o)


###################################################
### code chunk number 10: ovenc14of
###################################################
ovenc14o <- encodeOverlaps(grl14o, exbytx, hits=ov14)
ovenc14f <- encodeOverlaps(grl14f, exbytx, hits=ov14)


###################################################
### code chunk number 11: ovenc14
###################################################
grl14o_strand <- unlist(runValue(strand(grl14o)), use.names=FALSE)
ovenc14 <- selectEncodingWithCompatibleStrand(ovenc14o, ovenc14f,
                                              grl14o_strand, exbytx_strand,
                                              hits=ov14)
ovenc14


###################################################
### code chunk number 12: ovenc14_again
###################################################
ovenc14_again <- encodeOverlaps(grl14o, exbytx, hits=ov14, flip.query.if.wrong.strand=TRUE)
stopifnot(identical(ovenc14_again, ovenc14))


###################################################
### code chunk number 13: ovenc14_table
###################################################
unique_ovenc14 <- levels(encoding(ovenc14))
length(unique_ovenc14)
head(unique_ovenc14)
ovenc14_table <- table(encoding(ovenc14))
tail(sort(ovenc14_table))


###################################################
### code chunk number 14: compatible_unique_ovenc14
###################################################
sort(ovenc14_table[isCompatibleWithSplicing(unique_ovenc14)])


###################################################
### code chunk number 15: ov14_is_compat
###################################################
ov14_is_compat <- isCompatibleWithSplicing(ovenc14)
table(ov14_is_compat)  # 476124 "compatible" overlaps


###################################################
### code chunk number 16: gal14_ncompat
###################################################
gal14_ncompat <- tabulate(queryHits(ov14)[ov14_is_compat], nbins=length(gal14))
elementMetadata(gal14)$ncompat <- gal14_ncompat
head(gal14)
table(gal14_ncompat)


###################################################
### code chunk number 17: nb_non_zero_gal14_ncompat
###################################################
sum(gal14_ncompat != 0)


###################################################
### code chunk number 18: exbytx_ncompat14
###################################################
exbytx_ncompat14 <- tabulate(subjectHits(ov14)[ov14_is_compat], nbins=length(exbytx))
names(exbytx_ncompat14) <- names(exbytx)
tail(table(exbytx_ncompat14))


###################################################
### code chunk number 19: almost_compatible_unique_ovenc14
###################################################
sort(ovenc14_table[isCompatibleWithSkippedExons(unique_ovenc14)])


###################################################
### code chunk number 20: ov14_is_almostcompat
###################################################
ov14_is_almostcompat <- isCompatibleWithSkippedExons(ovenc14)
table(ov14_is_almostcompat)  # 837 "almost compatible" overlaps


###################################################
### code chunk number 21: gal14_nalmostcompat
###################################################
gal14_nalmostcompat <- tabulate(queryHits(ov14)[ov14_is_almostcompat], nbins=length(gal14))
elementMetadata(gal14)$nalmostcompat <- gal14_nalmostcompat
head(gal14)
table(gal14_nalmostcompat)


###################################################
### code chunk number 22: nb_non_zero_gal14_nalmostcompat
###################################################
sum(gal14_nalmostcompat != 0) 


###################################################
### code chunk number 23: exbytx_nalmostcompat14
###################################################
exbytx_nalmostcompat14 <- tabulate(subjectHits(ov14)[ov14_is_almostcompat], nbins=length(exbytx))
names(exbytx_nalmostcompat14) <- names(exbytx)
table(exbytx_nalmostcompat14)


###################################################
### code chunk number 24: aln_shows_nov_splice_jct
###################################################
aln_shows_nov_splice_jct <- gal14_nalmostcompat != 0L &
                            gal14_ncompat == 0L
head(which(aln_shows_nov_splice_jct))


###################################################
### code chunk number 25: is_nov_splice_jct
###################################################
is_nov_splice_jct <- queryHits(ov14) %in% which(aln_shows_nov_splice_jct)


###################################################
### code chunk number 26: narrow_is_nov_splice_jct
###################################################
is_nov_splice_jct <- is_nov_splice_jct & ov14_is_almostcompat


###################################################
### code chunk number 27: extractSkippedExonRanks
###################################################
skpexrk <- extractSkippedExonRanks(ovenc14)[is_nov_splice_jct]
table(elementLengths(skpexrk))


###################################################
### code chunk number 28: tx2skpexrk
###################################################
names(skpexrk) <- queryHits(ov14)[is_nov_splice_jct]
f <- names(exbytx)[subjectHits(ov14)[is_nov_splice_jct]]
tx2skpexrk <- split(skpexrk, f)


###################################################
### code chunk number 29: tx2skpexrk_names
###################################################
head(names(tx2skpexrk))  # transcript TxDb internal ids


###################################################
### code chunk number 30: tx10_details
###################################################
tx2skpexrk[["10"]]


###################################################
### code chunk number 31: tx58_details
###################################################
tx2skpexrk[["58"]]


###################################################
### code chunk number 32: untreated3_chr4
###################################################
untreated3_chr4()


###################################################
### code chunk number 33: readGappedAlignmentPairs
###################################################
galp34 <- readGappedAlignmentPairs(untreated3_chr4(), use.names=TRUE, param=param0)
head(galp34)


###################################################
### code chunk number 34: first_last
###################################################
head(first(galp34))
head(last(galp34))


###################################################
### code chunk number 35: isProperPair
###################################################
table(isProperPair(galp34))


###################################################
### code chunk number 36: keep_only_proper_pairs
###################################################
galp34 <- galp34[isProperPair(galp34)]


###################################################
### code chunk number 37: ov34
###################################################
ov34 <- findOverlaps(galp34, exbytx, ignore.strand=TRUE)


###################################################
### code chunk number 38: ovenc34
###################################################
grl34 <- grglist(galp34, order.as.in.query=TRUE)
ovenc34 <- encodeOverlaps(grl34, exbytx, hits=ov34, flip.query.if.wrong.strand=TRUE)
ovenc34


###################################################
### code chunk number 39: ovenc34_table
###################################################
unique_ovenc34 <- levels(encoding(ovenc34))
length(unique_ovenc34)
head(unique_ovenc34)
ovenc34_table <- table(encoding(ovenc34))
tail(sort(ovenc34_table))


###################################################
### code chunk number 40: compatible_unique_ovenc34
###################################################
sort(ovenc34_table[isCompatibleWithSplicing(unique_ovenc34)])


###################################################
### code chunk number 41: ov34_is_compat
###################################################
ov34_is_compat <- isCompatibleWithSplicing(ovenc34)
table(ov34_is_compat)  # 95801 "compatible" overlaps


###################################################
### code chunk number 42: galp34_ncompat
###################################################
galp34_ncompat <- tabulate(queryHits(ov34)[ov34_is_compat], nbins=length(galp34))
elementMetadata(galp34)$ncompat <- galp34_ncompat
head(galp34)
table(galp34_ncompat)


###################################################
### code chunk number 43: nb_non_zero_galp34_ncompat
###################################################
sum(galp34_ncompat != 0)


###################################################
### code chunk number 44: exbytx_ncompat34
###################################################
exbytx_ncompat34 <- tabulate(subjectHits(ov34)[ov34_is_compat], nbins=length(exbytx))
names(exbytx_ncompat34) <- names(exbytx)
tail(table(exbytx_ncompat34))


###################################################
### code chunk number 45: sessionInfo
###################################################
sessionInfo()


