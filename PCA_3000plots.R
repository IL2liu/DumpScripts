#######################################################################
# Plots PCAs of the chip profiles that have been filtered by wavelets #
#######################################################################

# datasets
# big table contains the combined CHIP profiles with ca.200,000 randomly picked sequences
raw_3000bin <- read.delim("~/Desktop/Projects/CHIPdatasets/HumanAlignmentData/fasq-cells:H1-hESC/3000_rawSelected/bigTable.bed")

# annotation file
#1. load the annotation file for the raw fastq files downloaded from Encode
annot_rawfastq <- read.csv(comm='#', "~/Desktop/Projects/CHIPdatasets/HumanAlignmentData/Fasq-H1-hESC_annotation.csv", header=F)
#2. load the the file for the sample names
sample_annot_rawfastq <- read.table("~/HumanAlignmentData/sample_fastqH1-hESC_annotation.txt", quote="\"")
#3. grep the file numbers
fastqFile_numbers = sub('.*H1-(\\d{3}).*', '\\1', annot_rawfastq$V11)
#4. grep only the files that were included in these datasets, 
#note there are three files were excluded as these files were failed with their md5sum check (during dowloading from Encode)
annot_Dataset_raw_3000 = paste('X', fastqFile_numbers, sep='') %in% colnames(raw_3000bin)
sample_Annot_raw_3000 = sample_annot_rawfastq$V2[annot_Dataset_raw_3000]
sampleNo_Annot_raw_3000 = fastqFile_numbers[annot_Dataset_raw_3000]

#data transformation
mat_raw_3000 = as.matrix(raw_3000bin[,ncol(raw_3000bin)])
mat_raw_log_3000 = log(mat_raw_3000+1)
pca_log_raw_3000 = prcomp(mat_raw_log_3000, scale.=T)

#PCA
gradient <- colorRampPalette(c('seagreen3', 'black', 'red'))(256)


for (i in 1:length(sampleNo_Annot_raw_3000)){
   file_no = sampleNo_Annot_raw_3000[i]
   pic_f = paste('pics', file_no, sep='/')
   pdf(paste(pic_f, '.pdf', sep=''), width=12)

   par(mfrow=c(1,2) ,oma=c(0,0,2,0))
   #PC1 and PC3
   plot(pca_log_raw_3000$x[,1], pca_log_raw_3000$x[,3], pch='.',
        col=gradient[256*rank(mat_raw_log_3000[,i])/nrow(mat_raw_log_3000)],
        xlab='PC1', ylab='PC3', frame=FALSE)
   #PC2 and PC5
   plot(pca_log_raw_3000$x[,2], pca_log_raw_3000$x[,5], pch='.',
        col=gradient[256*rank(mat_raw_log_3000[,i])/nrow(mat_raw_log_3000)],
        xlab='PC2', ylab='PC5', frame=FALSE)
   title(sample_Annot_raw_3000[i], outer=TRUE)
   dev.off()
}
