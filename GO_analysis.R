

# set global chunk options
# for figures
opts_chunk$set(fig.path='figs/', fig.align='center', fig.show='hold',dev='CairoPDF', out.width='.4\\linewidth')
# replacing "=" into "->" to make it R thing
options(replace.assign=TRUE,width=90)
# caching chunks
opts_chunk$set(cache.extra = R.version,cache.path='cache/')
opts_chunk$set(cache.extra = rand_seed)



## #library for mapping using ensembl database
## library(biomaRt)
## ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
## #input data containing all the gene_ids from TSSs that are associated with a color domain
## #this data set were obtained by
## colorAll_Id <- read.delim('input/geneID', header=F, sep='\ ')
## all_id <- colorAll_Id$V1
## all_goIDs = getBM(attributes=c('ensembl_gene_id', 'go_id'), filters='ensembl_gene_id', values=all_id, mart=ensembl)
## head(all_goIDs)
## write.table(all_goIDs, 'output/all_goIDs.txt',quote=FALSE, row.names=F, col.names=F, sep='\t')



## ColorGO_slim <- read.delim('output/Color_ExpressedGO_MajorityRule_GF20130505.txt', header=F)
## head(ColorGO_slim)
## rownames(ColorGO_slim) <- ColorGO_slim$V1
## ColorGO_slim <- ColorGO_slim[,2:5]
## colnames(ColorGO_slim) <- c('Red', 'Pink', 'Yellow', 'Black')
## 
## # add the backgound genes associated with colors
## backgroundGO_Color <- c(5799, 990, 1164, 1043)
## ColorGO_slim['Totals', ] <- backgroundGO_Color
## tail(ColorGO_slim)
## ColorGO_slim$row.sum <- rowSums(ColorGO_slim)
## 



## source('GO-compute-hyper-pvals.R')
## 
## # select the GO-terms < 0.001 FDR
## get_TOPGO <- function(GO_table, threshold){
##   #bilateral: - for depletion and + for enrichment
##   testHyperGO <- GO.pvals(GO_table, "bilateral")
##   adj_testHyperGO <- apply(testHyperGO, 2,
##                            function(x) p.adjust(abs(x), method="BH"))
##   row_index <- unique(which(adj_testHyperGO < threshold, arr.ind=TRUE)[,1])
##   TOP_p <- adj_testHyperGO[row_index,]
##   TOP_go <- GO_table[row_index,]
##   list_o <- list("TOP_p"= TOP_p,
##                  "TOP_go"= TOP_go)
##   return (list_o)
## }
## 
## # get the selected GO-terms with its associated adj_p values
## TOP_GO <- get_TOPGO(ColorGO_slim, 0.001)
## 



## TOP_GOgene <- TOP_GO$TOP_go
## head(TOP_GOgene)
## TOP_GOgene
## 
## # obs/exp
## ratTOP_GOgene <- apply(TOP_GOgene[, 1:4],1,  function(x) x/backgroundGO_Color)
## ratTOP_GOgene <- t(ratTOP_GOgene)
## # normalize to one
## NormTOP_GOgene <- sweep(ratTOP_GOgene, 1, rowSums(ratTOP_GOgene), FUN="/")



## 
## # loading's PCA plots for gene-go associations
## GeneGoBinary <- read.table('input/GeneGOBinary.txt',
##                            header=TRUE)
## GeneGoBinary <- as.matrix(t(GeneGoBinary))
## head(GeneGoBinary)
## colnames(GeneGoBinary) <- GeneGoBinary["GO_ID",]
## GeneGoBinary <- GeneGoBinary[-nrow(GeneGoBinary),]
## # make it a numeric matrix
## class(GeneGoBinary) <- "numeric"
## 
## # general PCA representation for gene-go associations
## pca_GeneGo <- prcomp(GeneGoBinary, scale.=T)
## plot(pca_GeneGo$rotation[,1:2],
##      cex=log(colSums(GeneGoBinary) ))
## dev.off()
## 
## # PCA for only enriched GO
## pca_GeneGORotation <- pca_GeneGo$rotation
## head(pca_GeneGORotation)
## dim(pca_GeneGORotation)
## rownames(pca_GeneGORotation)
## 
## NormTOP_GOgene <- sweep(ratTOP_GOgene, 1, rowSums(ratTOP_GOgene), FUN="/")
## head(NormTOP_GOgene)
## rownames(NormTOP_GOgene)
## 
## Enriched_pcaGeneGO <- pca_GeneGORotation[rownames(NormTOP_GOgene),1:2]
## dim(Enriched_pcaGeneGO)
## head(Enriched_pcaGeneGO)
## rownames(NormTOP_GOgene)
## 
## get_enrichGOColor <- function(enrichTable){
##   ratios <- vector()
##   colors_max <- vector()
##   for (i in 1: nrow(enrichTable)){
##     col_max <- which.max(enrichTable[i,])
##     colors_max[i] <- colnames(enrichTable)[col_max]
##     ratios[i] <- enrichTable[i, col_max]
##   }
##   list_o <- list("ratios"=ratios,
##                  "colors"=colors_max)
##   return (list_o)
## }
## 
## max_GOCOLOR <- get_enrichGOColor(NormTOP_GOgene)
## Enriched_pcaGeneGO <- as.data.frame(Enriched_pcaGeneGO)
## Enriched_pcaGeneGO$color <- max_GOCOLOR$colors
## Enriched_pcaGeneGO$ratio <- (max_GOCOLOR$ratios)*10
## Enriched_pcaGeneGO$no <- 1:nrow(Enriched_pcaGeneGO)
## size_enrichGO <- log(colSums(GeneGoBinary[,rownames(NormTOP_GOgene)]))
## 
## # get the mapping table
## goTermMap <- read.delim('input/GOtermMapTable.txt',
##                         header=FALSE)
## head(goTermMap)
## goTermMap[goTermMap$go_term == "mRNA processing",]
## colnames(goTermMap) <- c('go_id', 'go_term')
## Enriched_pcaGeneGOTerm <- merge(Enriched_pcaGeneGO,goTermMap,
##                                 by.x = "row.names",
##                                 by.y= "go_id")
## head(Enriched_pcaGeneGOTerm)
## factor(Enriched_pcaGeneGOTerm$color)
## 
## # for go_term annotation on plot
## selected_no <- c(7, 14, 12, 5)
## annotation_enrich <- Enriched_pcaGeneGOTerm[which(Enriched_pcaGeneGOTerm$no %in% selected_no),]
## 
## library(ggplot2)
## p <- ggplot(Enriched_pcaGeneGOTerm, aes(PC1, PC2))
## p <- p + geom_point(aes(colour=factor(color),
##                    size = size_enrichGO,
##                    ))   +
##         geom_text(data=annotation_enrich,
##                    aes(PC1,PC2,label=go_term,
##                        family="mono"),
##                    angle=45, vjust=-0.2,
##                    colour="azure4")
## 
## pdf('figs/PCA-GO.pdf', useDingbats=FALSE)
## 
## p + scale_colour_manual(values=c('black', 'deeppink2',
##                                  'red', 'gold2'))+
## scale_size(range = c(3, 15)) +
## theme(legend.position="none") +
## xlab('')+ylab('')
## 
## dev.off()
## 



## library("phenotypicForest")
## library(plyr)
## head(Enriched_pcaGeneGOTerm)
## head(NormTOP_GOgene)
## dim(NormTOP_GOgene)
## 
## library(reshape2)
## long.NormTop_GOgene <- melt(NormTOP_GOgene, variable.name="score")
## head(long.NormTop_GOgene)
## head(goTermMap)
## colnames(long.NormTop_GOgene) <- c("go_id", "color", "value")
## long.NormTop_GOgene <- merge(long.NormTop_GOgene, goTermMap,
##                              by.x="go_id", by.y="go_id" )
## head(long.NormTop_GOgene)
## long.NormTop_GOgene <- long.NormTop_GOgene[, c(4,4,2,3)]
## colnames(long.NormTop_GOgene) <- c("family", "item", "score", "value")
## attributes(long.NormTop_GOgene$family)$levels[1:20]
## 
## functions <-  as.vector(attributes(long.NormTop_GOgene$family)$levels[21:59])
## functions
## 
## core_funcs <- c("DNA metabolic proces","aging",
##                 "biosynthetic process","carbohydrate metabolic process",
##                 "catabolic process" ,"cell cycle",
##                 "cell death","cell division",
##                 "cell proliferation","cellular amino acid metabolic process" ,
##                 "cellular component assembly" , "cellular nitrogen compound metabolic process" ,
##                 "cellular protein modification process" ,"chromosome organization" ,
##                 "chromosome segregation",  "cofactor metabolic process",
##                 "cytoskeleton organization" , "cytoskeleton-dependent intracellular transport",
##                 "generation of precursor metabolites and energy", "growth" ,"lipid metabolic process",
##                 "locomotion", "mRNA processing" , "membrane organization",
##                 "mitochondrion organization", "mitosis" ,
##                 "nucleocytoplasmic transport" , "pigmentation",
##                 "plasma membrane organization", "protein complex assembly" ,"protein folding",
##                 "protein maturation", "protein targeting" , "reproduction",
##                 "response to stress", "ribonucleoprotein complex assembly" ,
##                 "ribosome biogenesis", "small molecule metabolic process",
##                 "sulfur compound metabolic process", "tRNA metabolic process" ,
##                 "translation" , "transport" ,
##                 "vacuolar transport" , "vesicle-mediated transport"
##                 )
## 
## specific_func <- c("anatomical structure development","anatomical structure formation involved in morphogenesis",
##                    "cell adhesion","cell differentiation","transmembrane transport",
##                    "cell morphogenesis", "cell motility",
##                    "cell-cell signaling", "developmental maturation",
##                    "embryo development", "extracellular matrix organization",
##                    "neurological system process", "signal transduction", "homeostatic process", "immune system process"
##                    )
## 
## family_GO <- factor(ifelse(long.NormTop_GOgene$family%in%specific_func, "specific", "core"))
## long.NormTop_GOgene$family <- family_GO
## head(long.NormTop_GOgene)
## factor(long.NormTop_GOgene$score)
## 
## 
## pdf("figs/polarGO.pdf", height=10, width=10, useDingbats=FALSE)
## p <- polarHistogram(long.NormTop_GOgene,
##                     guides=c(0,0,0,0))+
##      scale_fill_manual(values=c("black", "deeppink2",
##                                 "red", "gold2"))+
##     theme(legend.position="none")
## 
## dev.off()



## 
## head(Enriched_pcaGeneGOTerm)
## head(TOP_GOgene)
## head(size_enrichGO)
## 
## get_PC_EnrichGOTerm <- Enriched_pcaGeneGOTerm[, c(1,7,2,3)]
## # combine with the color-associated GOS
## spieDataSets <- merge(get_PC_EnrichGOTerm, TOP_GOgene,
##       by.x="Row.names", by.y="row.names")[, 1:8]
## # combine with the size
## spieDataSets <- merge(test, size_enrichGO,
##       by.x="Row.names", by.y="row.names")
## colnames(spieDataSets) <- c("GO_ID", "GO_Term",
##                             "X", "Y",
##                             "Red", "Pink", "Yellow", "Black",
##                             "Size")
## head(spieDataSets)
## spieDataSets$No <- 1:nrow(spieDataSets)
## 
## pdf("figs/PCAspie.pdf", useDingbats=FALSE,
##     width=8, height=8)
## library(scales)
## # plot the layout
## plot(c(min(spieDataSets$X)- 0.05, max(spieDataSets$X)+0.1),
##      c(min(spieDataSets$Y)-0.05, max(spieDataSets$Y)+0.1),
##      type='n', xlab="", ylab="",frame=FALSE)
## rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col =
##      "grey")
## # add the spiechart
## for (i in 1:nrow(spieDataSets)){
##   addspie(as.vector(backgroundGO_Color), as.numeric(spieDataSets[i, 5:8]),
##         x = spieDataSets$X[i], y=spieDataSets$Y[i],
##         R = sqrt(spieDataSets$Size[i])/100,
##         col=alpha(c("red", "deeppink2",
##               "gold2", "black"), 0.7),
##         xlab="", ylab="", xaxt="n", yaxt="n")
## }
## text(spieDataSets$X, spieDataSets$Y, labels = (spieDataSets$No), font=2,
##      pos=3, cex=1, col="white", adj=1, lwd=2)
## dev.off()
## 
## 



spieDataSets <- spieDataSets[, c(10,1,2)]
library(xtable)
head(spieDataSets)
print(xtable(spieDataSets), include.rownames=FALSE)



sessionInfo()



library(knitr)
knit("GO_analysis.Rnw" ) # compile to tex
purl("GO_analysis.Rnw", documentation = 0) # extract R code only
knit2pdf("GO_analysis.Rnw")


