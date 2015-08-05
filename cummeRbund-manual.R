### R code from vignette source 'cummeRbund-manual.Rnw'

###################################################
### code chunk number 1: init
###################################################
options(width=65)


###################################################
### code chunk number 2: loadLib
###################################################
library(cummeRbund)


###################################################
### code chunk number 3: read
###################################################
cuff <- readCufflinks(dir=system.file("extdata", package="cummeRbund"))
cuff


###################################################
### code chunk number 4: add_features
###################################################
#annot<-read.table("gene_annotation.tab",sep="\t",header=T,na.string="-")
#addFeatures(cuff,annot,level="genes")


###################################################
### code chunk number 5: global_plots_1
###################################################
dens<-csDensity(genes(cuff))
dens


###################################################
### code chunk number 6: global_plots_dens
###################################################
dens<-csDensity(genes(cuff))
dens
print(dens)


###################################################
### code chunk number 7: global_plots_2
###################################################
b<-csBoxplot(genes(cuff))
b


###################################################
### code chunk number 8: global_plots_box
###################################################
b<-csBoxplot(genes(cuff))
b
print(b)


###################################################
### code chunk number 9: global_plots_3
###################################################
s<-csScatter(genes(cuff),"hESC","Fibroblasts",smooth=T)
s


###################################################
### code chunk number 10: global_plots_scatter
###################################################
s<-csScatter(genes(cuff),"hESC","Fibroblasts",smooth=T)
s
print(s)


###################################################
### code chunk number 11: global_plots_4
###################################################
m<-MAplot(genes(cuff),"hESC","Fibroblasts")
m


###################################################
### code chunk number 12: global_plots_MA
###################################################
m<-MAplot(genes(cuff),"hESC","Fibroblasts")
m
print(m)


###################################################
### code chunk number 13: global_plots_5
###################################################
v<-csVolcano(genes(cuff),"hESC","Fibroblasts")
v


###################################################
### code chunk number 14: global_plots_volcano
###################################################
v<-csVolcano(genes(cuff),"hESC","Fibroblasts")
v
print(v)


###################################################
### code chunk number 15: data_access_1
###################################################
gene.features<-features(genes(cuff))
head(gene.features)
gene.fpkm<-fpkm(genes(cuff))
head(gene.fpkm)
isoform.fpkm<-fpkm(isoforms(cuff))
head(isoform.fpkm)
gene.diff<-diffData(genes(cuff))
head(gene.diff)


###################################################
### code chunk number 16: data_access_2
###################################################
sample.names<-samples(genes(cuff))
head(sample.names)
gene.featurenames<-featureNames(genes(cuff))
head(gene.featurenames)


###################################################
### code chunk number 17: data_access_3
###################################################
gene.matrix<-fpkmMatrix(genes(cuff))
head(gene.matrix)


###################################################
### code chunk number 18: create_geneset_1
###################################################
data(sampleData)
myGeneIds<-sampleIDs
myGeneIds
myGenes<-getGenes(cuff,myGeneIds)
myGenes


###################################################
### code chunk number 19: geneset_plots_1
###################################################
h<-csHeatmap(myGenes,cluster='both')
h


###################################################
### code chunk number 20: geneset_plots_heatmap
###################################################
h<-csHeatmap(myGenes,cluster='both')
h
print(h)


###################################################
### code chunk number 21: geneset_plots_1.5
###################################################
b<-expressionBarplot(myGenes)
b


###################################################
### code chunk number 22: geneset_plots_barplot
###################################################
b<-expressionBarplot(myGenes)
b
print(b)


###################################################
### code chunk number 23: geneset_plots_2
###################################################
s<-csScatter(myGenes,"Fibroblasts","hESC",smooth=T)
s


###################################################
### code chunk number 24: geneset_plots_scatter
###################################################
s<-csScatter(myGenes,"Fibroblasts","hESC",smooth=T)
s
print(s)


###################################################
### code chunk number 25: geneset_plots_3
###################################################
v<-csVolcano(myGenes,"Fibroblasts","hESC")
v


###################################################
### code chunk number 26: geneset_plots_volcano
###################################################
v<-csVolcano(myGenes,"Fibroblasts","hESC")
v
print(v)


###################################################
### code chunk number 27: geneset_plots_4
###################################################
ih<-csHeatmap(isoforms(myGenes),cluster='both',labRow=F)
ih


###################################################
### code chunk number 28: geneset_plots_isoform_heatmap
###################################################
ih<-csHeatmap(isoforms(myGenes),cluster='both',labRow=F)
ih
print(ih)


###################################################
### code chunk number 29: geneset_plots_5
###################################################
den<-csDendro(myGenes)


###################################################
### code chunk number 30: geneset_plots_dendro
###################################################
den<-csDendro(myGenes)
plot(den)


###################################################
### code chunk number 31: gene_level_1
###################################################
myGeneId<-"PINK1"
myGene<-getGene(cuff,myGeneId)
myGene
head(fpkm(myGene))
head(fpkm(isoforms(myGene)))


###################################################
### code chunk number 32: gene_plots_1
###################################################
gl<-expressionPlot(myGene)
gl


###################################################
### code chunk number 33: gene_plots_line
###################################################
gl<-expressionPlot(myGene)
gl
print(gl)


###################################################
### code chunk number 34: gene_plots_2
###################################################
gb<-expressionBarplot(myGene)
gb


###################################################
### code chunk number 35: gene_plots_bar
###################################################
gb<-expressionBarplot(myGene)
gb
print(gb)


###################################################
### code chunk number 36: gene_plots_3
###################################################
igb<-expressionBarplot(isoforms(myGene))
igb


###################################################
### code chunk number 37: gene_plots_bar_isoforms
###################################################
igb<-expressionBarplot(isoforms(myGene))
igb
print(igb)


###################################################
### code chunk number 38: get_sig_1
###################################################
mySigGenes<-getSig(cuff,alpha=0.05,level='genes')
head(mySigGenes)
length(mySigGenes)


###################################################
### code chunk number 39: get_sig_2
###################################################
hESC_vs_iPS.sigIsoforms<-getSig(cuff,x='hESC',y='iPS',alpha=0.05,level='isoforms')
head(hESC_vs_iPS.sigIsoforms)
length(hESC_vs_iPS.sigIsoforms)


###################################################
### code chunk number 40: get_sig_3
###################################################
mySigTable<-getSigTable(cuff,alpha=0.01,level='genes')
head(mySigTable,20)


###################################################
### code chunk number 41: geneset_plots_5
###################################################
ic<-csCluster(myGenes,k=4)
head(ic$cluster)
icp<-csClusterPlot(ic)
icp


###################################################
### code chunk number 42: geneset_plots_cluster
###################################################
ic<-csCluster(myGenes,k=4)
head(ic$cluster)
icp<-csClusterPlot(ic)
icp
print(icp)


###################################################
### code chunk number 43: specificity_1
###################################################
myGenes.spec<-csSpecificity(myGenes)
head(myGenes.spec)


###################################################
### code chunk number 44: similar_1
###################################################
mySimilar<-findSimilar(cuff,"PINK1",n=20)
mySimilar.expression<-expressionPlot(mySimilar,logMode=T,showErrorbars=F)


###################################################
### code chunk number 45: similar_plots_1
###################################################
mySimilar<-findSimilar(cuff,"PINK1",n=20)
mySimilar.expression<-expressionPlot(mySimilar,logMode=T,showErrorbars=F)
print(mySimilar.expression)


###################################################
### code chunk number 46: similar_2
###################################################
myProfile<-c(500,0,400)
mySimilar2<-findSimilar(cuff,myProfile,n=10)
mySimilar2.expression<-expressionPlot(mySimilar2,logMode=T,showErrorbars=F)


###################################################
### code chunk number 47: similar_plots_2
###################################################
myProfile<-c(500,0,400)
mySimilar2<-findSimilar(cuff,myProfile,n=10)
mySimilar2.expression<-expressionPlot(mySimilar2,logMode=T,showErrorbars=F)
print(mySimilar2.expression)


###################################################
### code chunk number 48: close_connection
###################################################
end<-sqliteCloseConnection(cuff@DB)


###################################################
### code chunk number 49: session
###################################################
sessionInfo()


