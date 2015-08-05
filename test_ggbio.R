library(ggbio)
options(width=72)
## theme_set(theme_grey(base_size = 9))


#genomic range objecy
library(ggbio)
set.seed(1)
N <- 200
library(GenomicRanges)
gr <- GRanges(seqnames = sample(c("chr1", "chr2", "chr3"),
                                size = N, replace = TRUE),
              IRanges(start = sample(1:300, size = N, replace = TRUE),
                      width = sample(70:75, size = N,replace = TRUE)),
              strand = sample(c("+", "-", "*"), size = N, 
                              replace = TRUE),
              value = rnorm(N, 10, 3), score = rnorm(N, 100, 30),
              group = sample(c("Normal", "Tumor"), 
                             size = N, replace = TRUE),
              pair = sample(letters, size = N, 
                            replace = TRUE))

# a  bug in the cbind!
methods:::bind_activation(FALSE)
gr
autoplot(gr)
autoplot(gr, stat = "coverage", geom = "area")



# test bed file into gRanges format
library(rtracklayer)
bedfilepath <- system.file("data", "input_bernsteinNA.bed", package = "rtracklayer")


bedfilepath
import(bedfilepath, asRangedData = FALSE)


library(rtracklayer)
ex_bed <- "data/input_bernsteinNA.bed"


lit = import("data/input_bernsteinNA_annote.bed")
head(lit)
autoplot(lit, layout = "karyogram", aes(color = name, fill = name))


test_path <- system.file("tests", package = "rtracklayer")
test_bed <- file.path(test_path, "test.bed")
head(test_bed)
autoplot(test_bed, aes(ﬁll = name))


test <- import(test_bed)
head(test)
autoplot(test)
qplot(test)

autoplot(ex_bed, coord="genome")


lit
autoplot(lit, aes(ﬁll = name))

head(lit)
autoplot(lit)


p_lit <- qplot(lit)
print(p_lit)
autoplot(lit)


library(rtracklayer)

autoplot("data/input_bernsteinNA.bed", aes(fill = name))


#############################
# plot with karyogram

data(hg19IdeogramCyto, package = "biovizBase")
head(hg19IdeogramCyto)
getOption("biovizBase")$cytobandColor
autoplot(hg19IdeogramCyto, layout = "karyogram", cytoband = gieStain)


#karyogram
library(biovizBase)
data(hg19IdeogramCyto)
p1 <- autoplot(hg19IdeogramCyto, layout = "karyogram", aes(fill  =
  gieStain))
p1

library(GenomicRanges)
hg19 <- keepSeqlevels(hg19IdeogramCyto, paste0("chr", c(1:22, "X", "Y")))
hg19
autoplot(hg19, layout = "karyogram", cytoband=FALSE)
biovizBase::isIdeogram(hg19)
autoplot(hg19, layout = "karyogram", cytoband = FALSE, aes(fill = gieStain)) +  scale_fill_giemsa()

data(darned_hg19_subset500, package = "biovizBase")
dn <- darned_hg19_subset500
head(dn)

data(hg19Ideogram, package = "biovizBase")

#get the sequence length
seqlengths(dn) <- seqlengths(hg19Ideogram)[names(seqlengths(dn))]
seqlengths(dn)
head(dn)
seqlengths(hg19Ideogram)[names(seqlengths(dn))]



autoplot("data/input_bernsteinNA.bed", layout = "karyogram", )

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ggbio)
data(genesymbol, package = "biovizBase")
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
model <- exonsBy(txdb, by = "tx")
model17 <- subsetByOverlaps(model, genesymbol["RBM17"])
exons <- exons(txdb)
exon17 <- subsetByOverlaps(exons, genesymbol["RBM17"])
## reduce to make sure there is no overlap just for example
exon.new <- reduce(exon17)
## suppose
set.seed(1)
values(exon.new)$sample1 <- rnorm(length(exon.new), 10, 3)
values(exon.new)$sample2 <- rnorm(length(exon.new), 10, 10)
values(exon.new)$significant <- c(TRUE, rep(FALSE, length(exon.new) - 1))


\section{Plot}
<<>>=
  require(ggbio)

#karyotype
bernstein <- autoplot("data/input_bernsteinNA.bed", layout = "karyogram", aes(color = name, fill = name)) 
bernstein

?read.bed
dn <- import("data/simpleRepeat.bed")
?import

# coerced into grange data
dn <- as(dn, "GRanges")
class(dn)
head(dn)


p <- bernstein + layout_karyogram(dn, aes(x = start, y = end), ylim = c(10, 30), geom = "line",
                                  color = "red")

p


methods:::bind_activation(FALSE)

myIdeo <- GRanges(seqnames = names(chrom.length), 
                  ranges = IRanges(start = 1, width = chrom.length))

seqlevels(myIdeo) = names(chrom.length)
seqlengths(myIdeo) = (chrom.length)




#############################################################################

crc1 <- system.file("extdata", "crc1-missense.csv", package = "biovizBase")
crc1 <- read.csv(crc1)
library(GenomicRanges)
mut.gr <- with(crc1, GRanges(Chromosome, IRanges(Start_position, End_position),
                             strand = Strand))

values(mut.gr) <- subset(crc1, select = -c(Start_position, End_position,
                                           Chromosome))
data("hg19Ideogram", package = "biovizBase")
seqs <- seqlengths(hg19Ideogram)
## subset_chr
chr.sub <- paste("chr", 1:22, sep = "")
chr.sub
## levels tweak
seqlevels(mut.gr) <- c(chr.sub, "chrX")
mut.gr <- keepSeqlevels(mut.gr, chr.sub)
seqs.sub <- seqs[chr.sub]
## remove wrong position
bidx <- end(mut.gr) <= seqs.sub[match(as.character(seqnames(mut.gr)), names(seqs.sub))]
mut.gr <- mut.gr[which(bidx)]
## assign_seqlengths
seqlengths(mut.gr) <- seqs.sub
## reanme to shorter names
new.names <- as.character(1:22)
names(new.names) <- paste("chr", new.names, sep = "")
new.names
mut.gr


hg19Ideo <- hg19Ideogram
hg19Ideo <- keepSeqlevels(hg19Ideogram, chr.sub)
hg19Ideo <- renameSeqlevels(hg19Ideo, new.names)
head(hg19Ideo)
rearr <- read.csv(system.file("extdata", "crc-rearrangment.csv", package = "biovizBase"))
rearr
hg19Ideo
theme_grey()
library(ggbio)
library("grid")
blankground <- function() {
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.margin = unit(0,"null"),
        plot.margin = rep(unit(0,"null"),4),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank()
        )
}

rm('p')
p <- ggplot()

p <- p+layout_circle(myIdeo, geom = "ideo", fill = "gray70", radius = 30,
                              trackWidth = 4)
p
