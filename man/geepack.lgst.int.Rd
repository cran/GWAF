\name{geepack.lgst.int}
\alias{geepack.lgst.int}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ function for testing gene-environment or gene-gene interaction between a dichotomous trait and a genotyped SNP in family data using Generalized Estimation Equation model}
\description{
  Fit logistic regression via Generalized Estimation Equation (GEE) to test gene-environment or gene-gene interaction between a dichotomous phenotype 
and one genotyped SNP in a genotype file in family data under additive genetic model. The interaction term is the product of SNP genotype and the covariate for interaction (\code{cov.int}). 
The covariate for interaction (\code{cov.int}) can be SNP genotype (gene-gene interaction) or an environmental factor (gene-environment interaction). Only one
interaction term is allowed. When \code{cov.int} is dichotomous, stratified analyses can be requested by specifying \code{sub}="Y". The covariance between the main 
effect (SNP) and the interaction effect is provided in the output when stratified analysis is not requested. Each family is treated as 
a cluster, with independence working correlation matrix used in the robust variance estimator. The interaction 
test is carried out by the \code{geese} function from package \code{geepack}. This function is called in \code{geepack.lgst.int.batch} function to apply interaction test to all SNPs in a 
genotype file.
}
\usage{
geepack.lgst.int(snp,phen,test.dat,covar,cov.int,sub="N")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{snp}{ genotype data of a SNP }
  \item{phen}{ a character string for a phenotype name in \code{test.dat} }
  \item{test.dat}{ the product of merging phenotype, genotype and pedigree data, should be ordered by "famid" }
  \item{covar}{ a character vector for covariates in \code{test.dat} }
  \item{cov.int}{ a character string naming the covariate for interaction, the covariate has to be included in \code{covar} }
  \item{sub}{ "N" (default) for no stratified analysis, and "Y" for requesting stratified analyses (only when cov.int is dichotomous) }
}
\details{
     The \code{geepack.lgst.int} function tests gene-environment or gene-genn interaction between a dichtomous trait and a SNP  
from a dataset that contains phenotype, genotype and pedigree data (\code{test.dat}), where the dataset needs to be ordered by famid. Please also see details in 
details for \code{geepack.lgst.int.batch} function.
}
\value{
  Please see value in \code{geepack.lgst.int.batch} function.
}
\references{   Liang, K.Y. and Zeger, S.L. (1986)
  Longitudinal data analysis using generalized linear models.
  \emph{Biometrika}, \bold{73} 13--22. 

  Zeger, S.L. and Liang, K.Y. (1986)
  Longitudinal data analysis for discrete and continuous outcomes.
  \emph{Biometrics}, \bold{42} 121--130.

  Yan, J and Fine, J. (2004) Estimating equations for association structures. \emph{Stat Med}, \bold{23} 859--874.  
}
\author{ Qiong Yang <qyang@bu.edu> and Ming-Huei Chen <mhchen@bu.edu>  }
\seealso{ \code{geese} function from package \code{geepack} }
