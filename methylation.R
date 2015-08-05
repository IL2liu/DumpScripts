ideoDMC <- function(methylDiff.obj, chrom.length, difference = 25, 
                    qvalue = 0.01, circos = FALSE, title = "test", hyper.col = "magenta", 
                    hypo.col = "green") {
  require(methylKit)
  require(GenomicRanges)
  require(ggbio)
  
  # chrom.length
  myIdeo <- GRanges(seqnames = names(chrom.length), ranges = IRanges(start = 1, 
                                                                     width = chrom.length))
  seqlevels(myIdeo) = names(chrom.length)
  seqlengths(myIdeo) = (chrom.length)
  
  
  hypo = get.methylDiff(methylDiff.obj, difference = difference, qvalue = qvalue, 
                        type = "hypo")
  hyper = get.methylDiff(methylDiff.obj, difference = difference, qvalue = qvalue, 
                         type = "hyper")
  
  g.per = as(hyper, "GRanges")
  g.per = keepSeqlevels(g.per, names(chrom.length))
  
  g.po = as(hypo, "GRanges")
  g.po = keepSeqlevels(g.po, names(chrom.length))
  values(g.po)$id = "hypo"
  values(g.per)$id = "hyper"
  
  if (circos) {
    
    p <- ggplot() + layout_circle(myIdeo, geom = "ideo", fill = "gray70", 
                                  radius = 39, trackWidth = 2)
    
    
    p <- p + layout_circle(c(g.po, g.per), geom = "point", 
                           size = 1, aes(x = midpoint, 
                                         y = meth.diff, color = id), radius = 25, trackWidth = 30) +              
                                           scale_colour_manual(values = c(hyper.col, hypo.col))
    p + layout_circle(myIdeo, geom = "text", aes(label = seqnames), 
                      vjust = 0, radius = 55, trackWidth = 7) + opts(title = title)
    
  } else {
    
    p <- ggplot() + layout_karyogram(myIdeo, cytoband = FALSE)
    p + layout_karyogram(c(g.po, g.per), geom = "point", size = 1, 
                         aes(x = midpoint, 
                             y = meth.diff, color = id)) + scale_colour_manual(values = c(hyper.col, 
                                                                                          hypo.col)) + opts(title = title)
    
  }
}