#date created: 131212

###########################
# Discretization pipeline #
###########################   

####################################################################
# !before run make the output directories , i.e, /pics and /tabs!  #
####################################################################

for (fname in list.files(pattern="*.bed")){
   # load the dataset
   data = read.delim(fname)
   HMMdata = data[,c(1,4,5)] # grab only chromosomes, CHIP, and input samples

   # the engine of discretizer ... here we go...
   require(HummingBee)
   fit <- NA
   try(fit <- BaumWelch.NB(HMMdata, m=2))
   if (is.na(fit)) {
      err_m = paste(fname,"failed", sep=":")
      write(err_m, stderr())
      next
   }

   # get the unique file-number for output's name
   file_no = sub('.*H1-(\\d{3}).*', '\\1', fname)
   pic_f = paste('pics', file_no, sep='/')
   pdf(paste(pic_f, '.pdf', sep=''))
   par(mfrow=c(2,1))
   plot(HMMdata[1:10000,2], col=fit$vPath, type='h', 
         ylab="read", xlab="chr1:1-10000", main=paste("CHIP", file_no, sep="-"))
   plot(HMMdata[1:10000,3], col=fit$vPath, type='h', 
         ylab="read", xlab="chr1:1-10000", main=paste("Input", file_no, sep="-"))
   dev.off()
   #write the table
   label = data[,1:3] #chr,start,end
   HMMoutput = cbind(label,fit$vPath-1)
   tab_f = paste('tabs', file_no, sep='/')
   write.table(HMMoutput, file=paste(tab_f,'.bed', sep=''), col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
}
