### R code from vignette source 'vignette.Rnw'

###################################################
### code chunk number 1: vignette.Rnw:27-57
###################################################
  library("phenotypicForest")
  library(reshape2)
  
  randomName<-function(n=1,syllables=3){
    vowels<-c("a","e","i","o","u","y")
    consonants<-letters[!(letters %in% vowels)]
    replicate(n,
              paste(
                rbind(sample(consonants,syllables,replace=TRUE),
                      sample(vowels,syllables,replace=TRUE)),
                sep='',collapse='')
              )
  }

  toyData<-function(nPhenotype,nSNP,nPhenotypeGroups){
    df<-data.frame(
      phenotype=rep(randomName(nPhenotype),1,each=nSNP),
      value=rep(1:nSNP,nPhenotype)+runif(nSNP*nPhenotype,min=-0.1,max=0.1),
      lowerBound=runif(nPhenotype*nSNP,min=0.0,max=0.1),
      upperBound=runif(nPhenotype*nSNP,min=0.0,max=0.1),
      phenotypeGroup=rep(sample(toupper(randomName(nPhenotypeGroups)),nPhenotype,replace=TRUE),1,each=nSNP),
      SNP=paste("rs",rep(sample(100000,nSNP),nPhenotype),sep='')
      )
    
    df<-within(df,{
      lowerBound<-value-lowerBound
      upperBound<-value+upperBound}
               )
    df  
  }


###################################################
### code chunk number 2: vignette.Rnw:60-66
###################################################
# creating a simple phorest
  set.seed(42)
  df<-toyData(17,7,4)
  print(head(df))
  p<-phorest(df,connectingLines=TRUE)
  print(p)


###################################################
### code chunk number 3: vignette.Rnw:73-74
###################################################
print(p)


###################################################
### code chunk number 4: vignette.Rnw:89-105
###################################################
 set.seed(42)
 df<-toyData(10,1,3)
 print(df)
 p1<-phorest(df)
 
 df$phenotypeGroup<-NULL # delete phenotypeGroup column
  
 userDefined<-list(
 "JONOSE"=c("xuxoho","xyhyma","suvaku","bogixi","byzylo","besero"),
 "RYBYGU"=c("sudiru","ziqiju"),
 "GYKUXU"=c("mapyxy","cyniky"))

  p2<-phorest(df,phenotypeGroups=userDefined)
  
  print(p1)
  print(p2)


###################################################
### code chunk number 5: vignette.Rnw:116-132
###################################################
# creating some random SNP groups
set.seed(42)
nSNP<-200
nSNPGroup<-4
df<-toyData(17,nSNP,4)
df$value<-rnorm(nrow(df))
tmp<-data.frame(
  SNP=unique(df$SNP),
  SNPGroup=sample(paste("SNP Group#",1:nSNPGroup,sep=''),nSNP,replace=TRUE))
df<-merge(df,tmp,by="SNP")

p<-phorest(
  df,
  largeSNPSet=TRUE,
  title='default plot for large SNP sets')
print(p)


###################################################
### code chunk number 6: vignette.Rnw:138-139
###################################################
print(p)


###################################################
### code chunk number 7: vignette.Rnw:144-153
###################################################
# specifying the aggregating functions
p<-phorest(
  df,
  largeSNPSet=TRUE,
  aggregatingFunction=function(x) mean(x,na.rm=TRUE),
  aggregatingLowerBound=function(x) mean(x,na.rm=TRUE)-sd(x,na.rm=TRUE),
  aggregatingUpperBound=function(x) mean(x,na.rm=TRUE)+sd(x,na.rm=TRUE),
  title='mean and standard deviation')
print(p)


###################################################
### code chunk number 8: vignette.Rnw:158-176
###################################################
library(reshape2)
set.seed(42)
nFamily<-20
nItemPerFamily<-sample(1:6,nFamily,replace=TRUE)
nValues<-3
randomWord<-function(n,nLetters=5) 
   replicate(n,paste(sample(letters,nLetters,replace=TRUE),sep='',collapse='')) 
 
 df<-data.frame(
   family=rep(randomWord(nFamily),nItemPerFamily),
   item=randomWord(sum(nItemPerFamily),3))
 
 df<-cbind(df,as.data.frame(matrix(runif(nrow(df)*nValues),nrow=nrow(df),ncol=nValues)))
 
 df<-melt(df,c("family","item"),variable.name="score") # from wide to long
 print(head(df))
 p<-polarHistogram(df,familyLabel=FALSE)
 print(p)


###################################################
### code chunk number 9: vignette.Rnw:181-182
###################################################
print(p)


###################################################
### code chunk number 10: vignette.Rnw:191-193
###################################################
  # not run
  # p<-phorest(df,columnNames=c("SNP"="snpid","phenotype"="assays))


###################################################
### code chunk number 11: vignette.Rnw:202-204
###################################################
  p<-p+opts(title="put title here")+xlab("label for x-axis")+ylab("label for y-axis")  
  print(p)


