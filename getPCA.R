###################################
# PCA on the randomly picked loci #
###################################


getwd()

get_dataSets <- function(dataset,inputs=FALSE, Colmat_dim, selected=FALSE, selection){
  
  # dataset : input table
  # inputs as a vector, e.g : c("X008","X195")
  # Colmat_dim = dimension of matrix substracted with the input samples
  # e.g dim of the dataset: 176861 97 with 2 inputs, it becomes 5:95
  # as the first 4 columns are chr, start, end ,and state
  
  # output:
  # a list of two datasets
  # 1. clean_dataset: annotated dataset (chr,start,end,state)
  # 2. mat_clean_dataset: dataset without annotation
  
  if (selected){
    dataset = dataset[, which(colnames(dataset) %in% selection)]
  } else{
  #remove the input
  dataset = dataset[, -which(colnames(dataset) %in% inputs)]
  }
  #remove the NA
  clean_dataset = dataset[complete.cases(dataset),]
  # remove rows that contain only zeros
  mat_clean_dataset = clean_dataset[,Colmat_dim]
  ind_nonZero = rowSums(mat_clean_dataset == 0)!= ncol(mat_clean_dataset) 
  mat_clean_dataset = mat_clean_dataset[ind_nonZero,]
  clean_dataset = clean_dataset[ind_nonZero,]
  returnList = list("clean_dataset"=clean_dataset, "mat_dataset"=mat_clean_dataset) 
  return (returnList)
  
}



# load datasets
dataset_91raw <- read.delim("na/randomlypicked_bigtable.raw")
dataset_91bin <- read.delim("randomlypicked_bigtable.bin")
selected_91raw <- read.delim("na/randomlypickedSelected_bigtable.raw")
selected_91bin <- read.delim("randomlypickedSelected_bigtable.bin")

head(dataset_91raw)
head(dataset_91bin)
head(selected_91raw)
head(selected_91bin)

dim(dataset_91raw)
dim(dataset_91bin)
dim(selected_91raw)
dim(selected_91bin)

raw91_Profiles <- get_dataSets(dataset_91raw, inputs=c("X008","X195"), 
                               Colmat_dim=5:95)
bin91_Profiles <-  get_dataSets(dataset_91bin, inputs=c("X008","X195"), 
                                Colmat_dim=5:95)
selected_raw91_Profiles <- get_dataSets(selected_91raw, Colmat_dim=5:14,
                                        selected=TRUE, 
                                        selection=c("chr","start", "end","state",
                                                    "X009", "X036", "X002", "X029", 
                                                    "X033", "X037", "X028", "X007", 
                                                    "X034", "X032"))
selected_bin91_Profiles <- get_dataSets(selected_91bin, Colmat_dim=5:14,
                                        selected=TRUE,
                                        selection=c("chr","start", "end","state",
                                                    "X009", "X036", "X002", "X029", 
                                                    "X033", "X037", "X028", "X007", 
                                                    "X034", "X032"),5:14)

head(raw91_Profiles$clean_dataset)
head(bin91_Profiles$clean_dataset)
head(selected_raw91_Profiles$clean_dataset)
head(selected_bin91_Profiles$clean_dataset)
dim(selected_raw91_Profiles$clean_dataset)

