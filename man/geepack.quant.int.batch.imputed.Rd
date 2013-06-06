\name{geepack.quant.int.batch.imputed}
\alias{geepack.quant.int.batch.imputed}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ function to test gene-environment or gene-gene interactions between a continuous trait and a batch of imputed SNPs in families using 
Generalized Estimation Equation model }
\description{
   Fit Generalized Estimation Equation (GEE) model to test gene-environment or gene-gene interactions for a continuous phenotype 
and all imputed SNPs in a genotype file in family data under additive genetic model. The interaction term is the product of SNP genotype (allelic dosage) and a covariate for interaction (\code{cov.int}). 
The covariate for interaction (\code{cov.int}) can be SNP genotype (gene-gene interaction) or an environmental factor (gene-environment interaction). Only one
interaction term is allowed. When \code{cov.int} is dichotomous, stratified analyses can be requested by specifying \code{sub}="Y". The covariance between the main 
effect (SNP) and the interaction effect is provided in the output when stratified analysis is not requested. Each pedigree is treated as 
a cluster, with independence working correlation matrix used in the robust variance estimator.
This function applies the same interaction test to all imputed SNPs in the genotype data. 
In each test for interaction, the \code{geese} function from \code{geepack} package is used.
}
\usage{
geepack.quant.int.batch.imputed(phenfile,genfile,pedfile,phen,covars,cov.int,sub="N",
outfile,col.names=T,sep.ped=",",sep.phe=",",sep.gen=",")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phenfile}{ a character string naming the phenotype file for reading (see format requirement in details) }
  \item{genfile}{ a character string naming the (imputed) genotype file for reading (see format requirement in details) }
  \item{pedfile}{ a character string naming the pedigree file for reading (see format requirement in details)}
  \item{outfile}{ a character string naming the result file for writing }
  \item{phen}{ a character string for a phenotype name in \code{phenfile} }
  \item{covars}{ a character vector for covariates in \code{phenfile} }
  \item{cov.int}{ a character string naming the covariate for interaction, the covariate has to be included in \code{covars} }
  \item{sub}{ "N" (default) for no stratified analysis, and "Y" for requesting stratified analyses (only when \code{cov.int} is dichotomous) }
  \item{col.names}{ a logical value indicating whether the output file should contain column names }
  \item{sep.ped}{ the field separator character for pedigree file }
  \item{sep.phe}{ the field separator character for phenotype file }
  \item{sep.gen}{ the field separator character for genotype file }
}
\details{
  Similar to the details for \code{geepack.quant.int.batch} function but here the SNP data contains imputed genotypes (allele dosages) 
that are continuous and range from 0 to 2. 
}
\value{
   Please see value in \code{geepack.quant.int.batch} function.
 
}
\references{   Liang, K.Y. and Zeger, S.L. (1986)
  Longitudinal data analysis using generalized linear models.
  \emph{Biometrika}, \bold{73} 13--22. 

  Zeger, S.L. and Liang, K.Y. (1986)
  Longitudinal data analysis for discrete and continuous outcomes.
  \emph{Biometrics}, \bold{42} 121--130.

  Yan, J and Fine, J. (2004) Estimating equations for association structures. \emph{Stat Med}, \bold{23} 859--874.  
}
\author{ Qiong Yang <qyang@bu.edu> and Ming-Huei Chen <mhchen@bu.edu> }

\examples{
\dontrun{
geepack.quant.int.batch.imputed(phenfile="simphen.csv",genfile="simgen.csv",
pedfile="simped.csv",phen="SIMQT",outfile="simout.csv",col.names=T,covars=c("sex",age"),
cov.int="sex",sub="Y",sep.ped=",",sep.phe=",",sep.gen=",")
}
}
