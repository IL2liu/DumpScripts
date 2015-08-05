# to run with source(...,chdir=TRUE)

filename<-"test-phorest.html"

# helpers
addHTML<-function(filename,html) cat(html,file=filename,append=TRUE)  

makeTable<-function(testID,ggplotObject,title='',scale=1,dpi=50){
  ggsave(paste("images/test_",testID,".png",sep=''),scale=scale,dpi=dpi)
  compare<-paste("compare images/truth_",testID,".png images/test_",testID,".png images/diff_",testID,".png",sep='')
  system(compare,ignore.stdout=FALSE,ignore.stderr=FALSE)
  
  convert<-paste("convert -delay 50 images/truth_",testID,".png images/test_",testID,".png images/diff_",testID,".png -loop 0 images/flicker_",testID,".gif",sep='')
  system(convert,ignore.stdout=FALSE,ignore.stderr=FALSE)
  
  html<-"<TABLE BORDER=1>"
  html<-paste(html,"<TR><TD COLSPAN=3>Test:",testID," ",title,"</TD></TR>",sep='')
  html<-paste(html,"<TR>")
  html<-paste(html,"<TD><IMG SRC='images/test_",testID,".png' width='300px' ></TD>",sep='')
  html<-paste(html,"<TD><IMG SRC='images/truth_",testID,".png' width='300px'></TD>",sep='')
  html<-paste(html,"<TD><A HREF='images/flicker_",testID,".gif'><IMG SRC='images/diff_",testID,".png' width='300px'></A></TD>",sep='')
  html<-paste(html,"</TR></TABLE>")
  html  
}
# random Name
# gives a random combination of syllables (consonant+vowel) to make up a random name. 
# Useful for generating readable fake IDs.
# C. Ladroue

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
    phenotypeGroup=rep(sample(randomName(nPhenotypeGroups),nPhenotype,replace=TRUE),1,each=nSNP),
    SNP=paste("rs",rep(sample(100000,nSNP),nPhenotype),sep='')
    )
  
  df<-within(df,{
    lowerBound<-value-lowerBound
    upperBound<-value+upperBound}
             )
  df  
}

# create file and header
cat("",file=filename,append=FALSE)
addHTML(filename,"<HTML>")
addHTML(filename,"<TITLE>Tests for phorest()</TITLE>")
addHTML(filename,"<BODY>")

# small SNP set
  addHTML(filename,"<H1>small SNP set</H1>")
  set.seed(42)
  nPhenotype<-15
  nPhenotypeGroups<-3
  nSNP<-10
  df<-toyData(nPhenotype,nSNP,nPhenotypeGroups)
  
  # missing values and shuffle
  df<-df[sample(nrow(df),floor(0.66*nPhenotype*nSNP)),]
  
  testID<-1
  p<-phorest(df)
  addHTML(filename,makeTable(testID,p,title=""))

  testID<-2
  p<-phorest(df,connectingLines=TRUE)
  addHTML(filename,makeTable(testID,p,title="connectingLines=TRUE"))

  # SNP Groups
  df<-toyData(nPhenotype,nSNP,nPhenotypeGroups)
  nSNPGroup<-3
  df$SNPGroup<-rep(paste("Group #",sort(1:nSNP%%nSNPGroups),sep=''),nPhenotype)

  df<-df[sample(nrow(df),floor(0.66*nPhenotype*nSNP)),]

  testID<-3
  p<-phorest(df,connectingLines=TRUE)
  addHTML(filename,makeTable(testID,p,title="with $SNPGroups"))

addHTML(filename,"<H1>Large SNP Set</H1>")

  set.seed(42)
  nPhenotype<-15
  nPhenotypeGroups<-3
  nSNP<-50
  df<-toyData(nPhenotype,nSNP,nPhenotypeGroups)
  nSNPGroup<-3
  df$SNPGroup<-rep(paste("Group #",sort(1:nSNP%%nSNPGroup),sep=''),nPhenotype)

  df<-df[sample(nrow(df),floor(0.66*nPhenotype*nSNP)),]

  testID<-4
  p<-phorest(df,largeSNPSet=TRUE)
  addHTML(filename,makeTable(testID,p,title=""))

  p<-phorest(df,
             largeSNPSet=TRUE,
             aggregatingFunction=function(x) mean(x,na.rm=TRUE),
             aggregatingUpperBound=function(x) max(x,na.rm=TRUE),
             aggregatingLowerBound=function(x) min(x,na.rm=TRUE)
             )
  testID<-5
  p<-phorest(df,largeSNPSet=TRUE)
  addHTML(filename,makeTable(testID,p,title="aggregating functions: min,mean,max"))

# closing file
addHTML(filename,"</BODY>")
addHTML(filename,"</HTML>")

# library("testthat")
# context("phorest, small SNP set")

# Hashes of ggplot objects is not stable enough to be used for tests
# library("digest")
# md5<-function(p) digest(serialize(p,NULL),algo='md5')
# 
# test_that("phorest, small SNP set ",{
#   set.seed(100)
#   nPhenotype<-15
#   nPhenotypeGroup<-3
#   nSNP<-6
#   df<-data.frame(phenotype=rep(LETTERS[1:nPhenotype],1,each=nSNP),
#                  SNP=paste("rs",rep(sample(100000,nSNP),nPhenotype,each=1),sep=''),
#                  value=rep(1:nSNP,nPhenotype,each=1)+runif(nPhenotype*nSNP)*0.8)
#   
#   df$lowerBound<-df$value-runif(nrow(df))*0.1
#   df$upperBound<-df$value+runif(nrow(df))*0.1
#   
#   
#   df<-df[sample(nrow(df),floor(0.66*nPhenotype*nSNP)),]
#   
#   p<-phorest(df,connectingLines=TRUE)
#   cat("\n",md5(p))
#   expect_equal(md5(p),"50f1dc56cdf2c134bd12734249a4d222")
#   
#   ## Grouping phenotypes
#   permutation<-sample(LETTERS[1:nPhenotype])
#   index<-c(1,sort(sample(2:nPhenotype,nPhenotypeGroup-1)),nPhenotype+1);
#   phenotypeGroup<-lapply(1:nPhenotypeGroup,function(k) permutation[index[k]:(index[k+1]-1)])
#   groupNames<-paste("group#",1:nPhenotypeGroup,sep='')
#   
#   p<-phorest(df,connectingLines=TRUE, groups=phenotypeGroup, groupNames=groupNames)
#   cat("\n",md5(p))
#   expect_equal(md5(p),"3ab07b66053b0427537cef0eca72140b")
# })
# 
# context("phorest, large SNP set")
# 
# test_that("phorest, large SNP set ",{
#   set.seed(100)
#   nPhenotype<-15
#   nPhenotypeGroup<-3
#   nSNP<-30
#   nSNPGroup<-4
#   
#   df<-data.frame(SNPGroup=rep(1:nSNPGroup,1,each=nSNP*nPhenotype))
#   SNPNames<-paste("rs",sample(100000,nSNP*nSNPGroup),sep='')
#   df<-within(df, {
#     SNP<-rep(SNPNames,1, each=nPhenotype)
#     value<-rnorm(nrow(df),mean=SNPGroup,sd=1/SNPGroup)
#     phenotype<-rep(LETTERS[1:nPhenotype],nSNP*nSNPGroup,each=1)
#     SNPGroup<-factor(paste("Group#",SNPGroup,sep=''))
#   })
#     
#   p<-phorest(df,largeSNPSet=TRUE)
#   cat("\n",md5(p))
#   
#   expect_equal(md5(p),"16e75f0305e088c70d5b04eead3702fc")
#   
#   ## Grouping phenotypes
#   set.seed(42)
#   permutation<-sample(LETTERS[1:nPhenotype])
#   index<-c(1,sort(sample(2:nPhenotype,nPhenotypeGroup-1)),nPhenotype+1);
#   phenotypeGroup<-lapply(1:nPhenotypeGroup,function(k) permutation[index[k]:(index[k+1]-1)])
#   groupNames<-paste("group#",1:nPhenotypeGroup,sep='')
#   
#   p<-phorest(df,connectingLines=TRUE, groups=phenotypeGroup, groupNames=groupNames,largeSNPSet=TRUE)
#  cat("\n",md5(p))
#   expect_equal(md5(p),"8e97586d9a7649fa8f1d9fbf40450b92")
#   })
