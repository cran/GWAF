
auto <- function(genopath,phenfile,pedfile,outfile,phen,covars=NULL,analysis,lib.loc,model=NULL,kinmat=NULL,col.names=F,sep.ped=",",sep.phe=",",sep.gen=","){

gfiles <- list.files(path=genopath,full.names=T)
gfs <- list.files(path=genopath,full.names=F)

if (length(gfs)==0) stop(paste("No genotype files in ",genopath,"!!",sep=""))
if (!analysis %in% c("lme","lme.imputed","gee","gee.imputed")) stop("Please choose one appropriate option from lme, lme.imputed, gee, and gee.imputed for analysis")
if (!file.exists(phenfile)) stop(paste(phenfile," does not exist!",sep=""))
if (!file.exists(pedfile)) stop(paste(pedfile," does not exist!",sep=""))
if (file.exists(outfile)) stop(paste(outfile," already exists!",sep=""))
if (!is.null(kinmat)) {
   trykin <- try(load(kinmat))
   if (inherits(trykin,"try-error")) stop(paste('kinship matrix does not exist at ',kinmat))
}

pheno <- read.table(phenfile,as.is=T,header=T,sep=sep.phe)
if (!phen %in% names(pheno)) stop(paste(phen," does not exist in ",phenfile,"!!",sep=""))
if (!is.null(covars) & sum(covars %in% names(pheno))!=length(covars)) stop(paste("Not all covariates exist in ",phenfile,"!!",sep=""))
when <- format(Sys.time(), "%Y-%b-%d-%H-%M")

if (analysis=="lme") {
   cmds <- character(8)
   if (is.null(model)) model <- "a"
   for (j in 1:length(gfs)){
       cmds[1] <- paste("library(GWAF,lib.loc='",lib.loc,"')",sep="")
       cmds[2] <- paste("lme.batch(phenfile='",phenfile,"',genfile='",gfiles[j],"',",sep="")
       cmds[3] <- paste("pedfile='",pedfile,"',",sep="")
       cmds[4] <- paste("outfile='",outfile,"',",sep="")
       if (is.null(covars)) cmds[5] <- paste("kinmat='",kinmat,"',",sep="") else 
           cmds[5] <- paste("covars=c('",noquote(paste(covars,collapse="','")),"'),kinmat='",kinmat,"',",sep="")                   
       cmds[6] <- paste("phen='",phen,"',model='",model,"',sep.ped='",sep.ped,"',sep.phe='",sep.phe,"',sep.gen='",sep.gen,"',col.names=F)",sep="")
       cmds[7] <- paste("write(paste('",analysis," analysis is done for ",phen," with ",gfiles[j],"!',sep=''),paste('",phen,".lme.",when,".",j,"_",gfs[j],".log',sep=''))",sep="")
       cmds[8] <- paste("system(paste('rm ",phen,".lme.",when,"*.R ",phen,".lme.",when,"*.sh',sep=''))",sep="")
       write(cmds,paste(phen,".lme.",when,".",j,".R",sep=""),ncol=1)
       write("#$ -o /dev/null",file=paste(phen,".lme.",when,".",j,".sh",sep=""))
       write("#$ -e /dev/null",file=paste(phen,".lme.",when,".",j,".sh",sep=""),append=T)
       write(paste("R < ",phen,".lme.",when,".",j,".R --no-save",sep=""),file=paste(phen,".lme.",when,".",j,".sh",sep=""),append=T)
       write(paste("qsub -cwd ",phen,".lme.",when,".",j,".sh",sep=""),paste(phen,".lme.",when,".lst",sep=""),ncol=1,append=T)
   }
   if (model %in% c("a","d","r")){
  	write(c("phen","snp","n0","n1","n2","h2q","beta","se","chisq","df","model","pval"),outfile,sep=",",ncol=12)
   } else write(c("phen","snp","n0","n1","n2","h2q","beta10","beta20","beta21","se10","se20","se21","chisq","df","model","pval"),outfile,ncol=16,sep=",")
   print(paste("Quit R and submit all jobs by using ksh ",phen,".lme.",when,".lst",sep=""))

}  else 
if (analysis=="lme.imputed") {
   cmds <- character(8)
   if (!is.null(model)) stop(paste("No model option for imputed analysis!!"))
   for (j in 1:length(gfs)){
       cmds[1] <- paste("library(GWAF,lib.loc='",lib.loc,"')",sep="")
       cmds[2] <- paste("lme.batch.imputed(phenfile='",phenfile,"',genfile='",gfiles[j],"',",sep="")
       cmds[3] <- paste("pedfile='",pedfile,"',",sep="")
       cmds[4] <- paste("outfile='",outfile,"',",sep="")
       if (is.null(covars)) cmds[5] <- paste("kinmat='",kinmat,"',",sep="") else 
          cmds[5] <- paste("covars=c('",noquote(paste(covars,collapse="','")),"'),kinmat='",kinmat,"',",sep="")                   
       cmds[6] <- paste("phen='",phen,"',sep.ped='",sep.ped,"',sep.phe='",sep.phe,"',sep.gen='",sep.gen,"',col.names=F)",sep="")
       cmds[7] <- paste("write(paste('",analysis," analysis is done for ",phen," with ",gfiles[j],"!',sep=''),paste('",phen,".lme.imputed.",when,".",j,"_",gfs[j],".log',sep=''))",sep="")
       cmds[8] <- paste("system(paste('rm ",phen,".lme.imputed.",when,"*.R ",phen,".lme.imputed.",when,"*.sh',sep=''))",sep="")
       write(cmds,paste(phen,".lme.imputed.",when,".",j,".R",sep=""),ncol=1)
       write("#$ -o /dev/null",file=paste(phen,".lme.imputed.",when,".",j,".sh",sep=""))
       write("#$ -e /dev/null",file=paste(phen,".lme.imputed.",when,".",j,".sh",sep=""),append=T)
       write(paste("R < ",phen,".lme.imputed.",when,".",j,".R --no-save",sep=""),file=paste(phen,".lme.imputed.",when,".",j,".sh",sep=""),append=T)
       write(paste("qsub -cwd ",phen,".lme.imputed.",when,".",j,".sh",sep=""),paste(phen,".lme.imputed.",when,".lst",sep=""),ncol=1,append=T)
   }
   write(c("phen","snp","N","AF","h2q","beta","se","pval"),outfile,ncol=8,sep=",")
   print(paste("Quit R and submit all jobs by using ksh ",phen,".lme.imputed.",when,".lst",sep=""))

} else
if (analysis=="gee") {       
   cmds <- character(8)
   if (is.null(model)) model <- "a"
   for (j in 1:length(gfs)){
       cmds[1] <- paste("library(GWAF,lib.loc='",lib.loc,"')",sep="")
       cmds[2] <- paste("gee.lgst.batch(phenfile='",phenfile,"',genfile='",gfiles[j],"',",sep="")
       cmds[3] <- paste("pedfile='",pedfile,"',",sep="")
       cmds[4] <- paste("outfile='",outfile,"',",sep="")
       if (is.null(covars)) cmds[5] <- paste("',model='",model,"',",sep="") else 
           cmds[5] <- paste("covars=c('",noquote(paste(covars,collapse="','")),"'),model='",model,"',",sep="")                   
       cmds[6] <- paste("phen='",phen,"',sep.ped='",sep.ped,"',sep.phe='",sep.phe,"',sep.gen='",sep.gen,"',col.names=F)",sep="")
       cmds[7] <- paste("write(paste('",analysis," analysis is done for ",phen," with ",gfiles[j],"!',sep=''),paste('",phen,".gee.",when,".",j,"_",gfs[j],".log',sep=''))",sep="")
       cmds[8] <- paste("system(paste('rm ",phen,".gee.",when,"*.R ",phen,".gee.",when,"*.sh',sep=''))",sep="")
       write(cmds,paste(phen,".gee.",when,".",j,".R",sep=""),ncol=1)
       write("#$ -o /dev/null",file=paste(phen,".gee.",when,".",j,".sh",sep=""))
       write("#$ -e /dev/null",file=paste(phen,".gee.",when,".",j,".sh",sep=""),append=T)
       write(paste("R < ",phen,".gee.",when,".",j,".R --no-save",sep=""),file=paste(phen,".gee.",when,".",j,".sh",sep=""),append=T)
       write(paste("qsub -cwd ",phen,".gee.",when,".",j,".sh",sep=""),paste(phen,".gee.",when,".lst",sep=""),ncol=1,append=T)
   }
  if (model %in% c("a","d","r")) {  	
	write(c("phen","snp","n0","n1","n2","nd0","nd1","nd2","miss.0","miss.1","miss.diff.p","beta","se","chisq","df","model","remark","pval"),
             outfile,sep=",",ncol=18)
  } else
	write(c("phen","snp","n0","n1","n2","nd0","nd1","nd2","miss.0","miss.1","miss.diff.p","beta10","beta20","beta21",
			"se10","se20","se21","chisq","df","model","remark","pval"),outfile,sep=",",ncol=22)
   print(paste("Quit R and submit all jobs by using ksh ",phen,".gee.",when,".lst",sep=""))

}  else 
if (analysis=="gee.imputed") {       
   cmds <- character(7)
   if (!is.null(model)) stop(paste("No model option for imputed analysis!!"))
   for (j in 1:length(gfs)){
       cmds[1] <- paste("library(GWAF,lib.loc='",lib.loc,"')",sep="")
       cmds[2] <- paste("gee.lgst.batch.imputed(phenfile='",phenfile,"',genfile='",gfiles[j],"',",sep="")
       cmds[3] <- paste("pedfile='",pedfile,"',",sep="")
       cmds[4] <- paste("outfile='",outfile,"',sep.ped='",sep.ped,"',sep.phe='",sep.phe,"',sep.gen='",sep.gen,"',col.names=F,",sep="")
       if (is.null(covars)) cmds[5] <- paste("phen='",phen,"')",sep="") else 
           cmds[5] <- paste("covars=c('",noquote(paste(covars,collapse="','")),"'),phen='",phen,"')",sep="")                   
       cmds[6] <- paste("write(paste('",analysis," analysis is done for ",phen," with ",gfiles[j],"!',sep=''),paste('",phen,".gee.imputed.",when,".",j,"_",gfs[j],".log',sep=''))",sep="")
       cmds[7] <- paste("system(paste('rm ",phen,".gee.imputed.",when,"*.R ",phen,".gee.imputed.",when,"*.sh',sep=''))",sep="") 
       write(cmds,paste(phen,".gee.imputed.",when,".",j,".R",sep=""),ncol=1)
       write("#$ -o /dev/null",file=paste(phen,".gee.imputed.",when,".",j,".sh",sep=""))
       write("#$ -e /dev/null",file=paste(phen,".gee.imputed.",when,".",j,".sh",sep=""),append=T)
       write(paste("R < ",phen,".gee.imputed.",when,".",j,".R --no-save",sep=""),file=paste(phen,".gee.imputed.",when,".",j,".sh",sep=""),append=T)
       write(paste("qsub -cwd ",phen,".gee.imputed.",when,".",j,".sh",sep=""),paste(phen,".gee.imputed.",when,".lst",sep=""),ncol=1,append=T)
   }
   write(c("phen","snp","N","Nd","AF","AFd","beta","se","remark","pval"),outfile,ncol=10,sep=",")
   print(paste("Quit R and submit all jobs by using ksh ",phen,".gee.imputed.",when,".lst",sep=""))
}  

}

