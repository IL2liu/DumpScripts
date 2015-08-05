########################################
# draw chromatin states on chromosomes #
########################################

library(ggbio)
## requrie 
p <- plotIdeogram()




p <- plotIdeogram(genome = "hg19")
p
attr(p, "ideogram.data")
plotIdeogram(hg19IdeogramCyto, "chr1")

library(biovizBase)
data(hg19IdeogramCyto)
plotIdeogram(hg19IdeogramCyto, "chr1", aspect.ratio = 1/20)
plotIdeogram(hg19IdeogramCyto, "chr1", aspect.ratio = 1/20, zoom.region = c(1e+07,5e+07))

p <- plotIdeogram(hg19IdeogramCyto, "chr1", aspect.ratio = 1/20)
p + xlim(1e+07, 5e+07)
p <- plotIdeogram(hg19IdeogramCyto, "chr1")
df <- data.frame(x = seq(from = 5e+07, to = 9e+07, length = 100), y = rnorm(100))
p2 <- qplot(data = df, x = x, y = y, geom = "line") + ylab("")
tracks(p, p2 = p2, heights = c(1.2, 5))
p1 <- plotIdeogram(mm9, "chr1")
