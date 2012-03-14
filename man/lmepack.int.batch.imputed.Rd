\name{lmepack.int.batch.imputed}
\alias{lmepack.int.batch.imputed}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ function to test gene-environment or gene-gen interaction and associations between a continuous trait and a batch of imputed SNPs in families using 
Linear Mixed Effects model }
\description{
   Fit linear mixed effect model to test gene-environment or gene-gen interaction and genetic association for a continuous phenotype 
and all imputed SNPs in a genotype file under additive genetic model. The interaction term is the product of SNP genotype (allelic dosage) and a covariate for interaction (\code{cov.int}). 
The covariate for interaction (\code{cov.int}) can be SNP genotype (gene-gene interaction) or an environmental factor (gene-environment interaction). Only one
interaction term is allowed. When (\code{cov.int}) is dichotomous, stratified analyses can be requested by specifying \code{sub}="Y". The covariance between the main 
effect (SNP) and the interaction effect is provided in the output when stratified analysis is not requested. The SNP genotype and the interaction are treated
as fixed effectd, and a random effect correlated according to degree of relatedness within a family is also fitted. 
In each test for trait-SNP assocaition or interaction, the \code{lmekin()} function which is modified from
the same named function in package \code{kinship} is used. As \code{kinship} package has been archived, \code{GWAF} package now
includes functions from \code{kinship} package to support the modified \code{lmekin} function. Archived \code{kinship} package can be found 
at http://cran.r-project.org/web/packages/kinship/index.html.
}
\usage{
lmepack.int.batch.imputed(phenfile,genfile,pedfile,phen,kinmat,covars,cov.int,sub="N",
outfile,col.names=T,sep.ped=",",sep.phe=",",sep.gen=",")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phenfile}{ a character string naming the phenotype file for reading (see format requirement in details) }
  \item{genfile}{ a character string naming the (imputed) genotype file for reading (see format requirement in details) }
  \item{pedfile}{ a character string naming the pedigree file for reading (see format requirement in details)}
  \item{outfile}{ a character string naming the result file for writing }
  \item{phen}{ a character string for a phenotype name in phenfile }
  \item{covars}{ a character vector for covariates in phenfile }
  \item{cov.int}{ a character string naming the covariate for interaction, the covariate has to be included in \code{covars} }
  \item{sub}{ "N" (default) for no stratified analysis, and "Y" for requesting stratified analyses (only when cov.int is dichotomous) }
  \item{kinmat}{ a character string naming the file where kinship coefficient matrix is kept }
  \item{col.names}{ a logical value indicating whether the output file should contain column names }
  \item{sep.ped}{ the field separator character for pedigree file }
  \item{sep.phe}{ the field separator character for phenotype file }
  \item{sep.gen}{ the field separator character for genotype file }
}
\details{
  Similar to the details for 'lmepack.int.batch' function but here the SNP data contains imputed genotypes (allele dosages) 
that are continuous and range from 0 to 2. 
}
\value{
   Please see value in 'lmepack.int.batch' function.
 
}
\references{ 
 kinship package: mixed-effects Cox models, sparse matrices, and modeling data from large pedigrees.
Beth Atkinson (atkinson@mayo.edu) for pedigree functions.Terry Therneau (therneau@mayo.edu) for all other functions.
2007. Ref Type: Computer Program http://cran.r-project.org/. 

 Abecasis, G. R., Cardon, L. R., Cookson, W. O., Sham, P. C., & Cherny, S. S. Association analysis in 
a variance components framework. \emph{Genet Epidemiol}, \bold{21} Suppl 1, S341-S346 (2001).
}
\author{ Qiong Yang <qyang@bu.edu> and Ming-Huei Chen <mhchen@bu.edu> }