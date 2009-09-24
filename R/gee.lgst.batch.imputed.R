#dichotomous traits only
#phenfile: phenotype file name in quotation marks,must provide
#genfile" genotype file name in quotation marks, must provide
#outfile: output file name in quotation marks,must provide
#pedfile: famid id fa mo sex 
#famid is cluster id
gee.lgst.batch.imputed=function(genfile,phenfile,pedfile,outfile,phen,covars=NULL){
  print(paste("phenotype data = ", phenfile))
  print(paste("genotype data = ", genfile))
  print(paste("pedigree data = ", pedfile))
  print(paste("Result of GEE analyses =",outfile))
  if(is.null(covars)){
    print("Covariates = NONE")
  }else{
    print(paste("Covariate(s) =",covars,collapse=", "))
  }
  print("Running GEE")

  read.in.data <- function(phenfile,genfile,pedfile) {
  print("Reading in Data")
  ped.dat <- read.csv(gzfile(genfile),header=TRUE,na.string="")
  snp.names <- names(ped.dat)[-1]
  pedigree <- read.csv(pedfile,header=TRUE)
  gntp.all <- merge(pedigree,ped.dat,by="id")

#read in phenotype data
  phen.dat=read.csv(phenfile,header=TRUE)
  phen.name=colnames(phen.dat)[-1]
  n.snp=length(names(gntp.all))

  if(length(grep("^sex$",colnames(phen.dat)))==0) {
  phensnp.dat<-merge(gntp.all,phen.dat,by=c("id"))
  } else {
## sex is one of the columns in the phenotype file
  phensnp.dat<-merge(gntp.all,phen.dat,by=c("id","sex"))
  }
  print("Done reading in data")
  return(list(data=phensnp.dat,snps=snp.names,phen.name=phen.name))
}

  phensnp.dat <- read.in.data(phenfile,genfile,pedfile)
  snplist<-phensnp.dat$snps
  test.dat<-phensnp.dat$data
  if (!is.null(covars) & sum(snplist %in% covars)>=1) {
     names(test.dat)[which(names(test.dat)==paste(snplist[snplist %in% covars],".x",sep=""))] <- snplist[snplist %in% covars]
     covars[covars %in% snplist] <- paste(covars[covars %in% snplist],".y",sep="")
  }
  test.dat<-test.dat[order(test.dat$famid),]
  library(gee)

  final1<-c()
  if (is.null(covars)){
     print(paste("No Covariates, Running:",phen))
     if (length(snplist)<2) { 
        temp.out <- c(phen,snplist,gee.lgst.imputed(snp=test.dat[,snplist],phen=phen,test.dat=test.dat))
     } else {
        temp.out <-as.data.frame(apply(test.dat[,phensnp.dat$snps],2,gee.lgst.imputed,phen=phen,test.dat=test.dat))
  	 temp.out <-cbind(rep(phen,ncol(temp.out)),colnames(temp.out),t(temp.out))
       }
     final1 <- rbind(final1,temp.out)
  } else {
     print(paste("Covariates, Running:",phen))
     if (length(snplist)<2) { 
        temp.out <- c(phen,snplist,gee.lgst.imputed(snp=test.dat[,snplist],phen=phen,test.dat=test.dat,covar=covars))
     } else {
        temp.out <-as.data.frame(apply(test.dat[,phensnp.dat$snps],2,gee.lgst.imputed,phen=phen,test.dat=test.dat,covar=covars))
  	 temp.out <-cbind(rep(phen,ncol(temp.out)),colnames(temp.out),t(temp.out))
       }
     final1 <- rbind(final1,temp.out)
  }
	
  colnames(final1)=c("phen","snp","N","Nd","AF","AFd","beta","se","remark","pval")
  write.table(as.matrix(final1),outfile,col.names=T, row.names=F,quote=F,sep=",",na="")
  warnings()
}