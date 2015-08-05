###########################################
# Annotate your genome in circos's style  #
###########################################

circosHuman_GenomeAnnotation <- function(in.chrAnnot, out.chrAnnot, 
                                         title, in.lab, out.lab)
{ 
  library(ggbio)
  library(rtracklayer) #for import purpose
  
  # get the hg19 ideogram
  library("BSgenome.Hsapiens.UCSC.hg19")
  chr.len = seqlengths(Hsapiens)  # get chromosome lengths
  # remove X,Y,M and random chromosomes
  chr.len = chr.len[grep("_|M", names(chr.len), invert = T)] 
  names(chr.len)
  hg19Ideo <- GRanges(seqnames = names(chr.len), 
                      ranges = IRanges(start = 1,width = chr.len))
  seqlevels(hg19Ideo) = names(chr.len)
  seqlengths(hg19Ideo) = (chr.len)
  
  shortChr.names <- as.character(1:22)
  names(shortChr.names) <- paste("chr", shortChr.names, sep = "")
  shortChr.names  
  hg19Ideo <- renameSeqlevels(hg19Ideo, shortChr.names)
  
  # load dataset to annotate inside chromosome
  in.chr <- import(in.chrAnnot)
  in.chr <- as(in.chr,"GRanges")
  
  # adjust the chr size
  in.chr = keepSeqlevels(in.chr, names(chr.len))
  seqlevels(in.chr) = names(chr.len)
  seqlengths(in.chr) = (chr.len)
  
  # load dataset to annotate outside chromosome
  out.chr <- import(out.chrAnnot)
  out.chr <- as(out.chr, "GRanges")
  
  # adjust the chr size
  out.chr = keepSeqlevels(out.chr, names(chr.len))
  seqlevels(out.chr) = names(chr.len)
  seqlengths(out.chr) = (chr.len)
  
  # circo in action!
  p <- ggplot() + layout_circle(hg19Ideo, geom = "ideo", fill = "gray70", 
                                radius = 30,trackWidth = 4)
  p <- p + layout_circle(hg19Ideo, geom = "scale", size = 2, radius = 35, 
                         trackWidth = 2)
  p <- p + layout_circle(hg19Ideo, geom = "text", aes(label = seqnames), 
                         vjust = 0,radius = 38, trackWidth = 7)
  
  # annotate inside chromosomes
  p <- p + layout_circle(in.chr, geom = "rect", color = "black", 
                         radius = 30,trackWidth = 4)
  # annotate outside chromosomes
  p <- p + layout_circle(out.chr, geom = "rect", color = "steelblue",
                         radius = 26,trackWidth = 4)
  
  # annotate with text
  p <- p + labs(title = title) + 
      annotate("text",x=50,y=50,label=in.lab, color= "black", size=3.0, hjust=1)+
      annotate("text",x=50,y=47,label=out.lab, color= "steelblue", size=3.0, hjust=1)
    
  return (p)
}

###################
# working example #
###################
# set the working directory
parent_dir = '/TEMP_DDN/users/gfilion/rlim/Desktop/Projects/CHIPdatasets/HumanAlignmentData/fasq-cells:H1-hESC'
r_dir = 'Ranalysis/check_gemMappingNoRepeat'
work_dir = file.path(parent_dir, r_dir)
setwd(work_dir)

# draw circos
input_bernstein <- circosHuman_GenomeAnnotation(in.chrAnnot="data/randomly1000PickedSimpleRepeat.bed",
                             out.chrAnnot="data/randomly1000PickedBernsteinNA.bed",
                            in.lab="Simple Repeats", out.lab="Gem Repeats",
                            title="Input-H1hESC-Bernstein")

input_synder <- circosHuman_GenomeAnnotation(in.chrAnnot="data/randomly1000PickedSimpleRepeat.bed",
                                             out.chrAnnot="data/randomly1000PickedSynderRepeat.bed",
                                             in.lab="Simple Repeats", out.lab="Gem Repeats",
                                             title="Input-H1hESC-Synder")

par(mfrow=c(2,2), oma=c(0,0,2,0))
input_bernstein
input_synder
