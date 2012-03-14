\name{geepack.lgst.int.imputed}
\alias{geepack.lgst.int.imputed}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ function for testing gene-environment or gene-gene interaction between a dichotomous trait and an imputed SNP in family data using Generalized Estimation Equation model}
\description{
  Fit logistic regression via Generalized Estimation Equation (GEE) to test gene-environment or gene-gene interaction between a dichotomous phenotype 
and one imputed SNP in a genotype file under additive genetic model. The interaction term is the product of SNP genotype and a covariate for interaction (\code{cov.int}). 
The covariate for interaction (\code{cov.int}) can be SNP genotype (gene-gene interaction) or an environmental factor (gene-environment interaction). Only one
interaction term is allowed. When \code{cov.int} is dichotomous, stratified analyses can be requested by specifying \code{sub}="Y". The covariance between the main 
effect (SNP) and the interaction effect is provided in the output when stratified analysis is not requested. Each family is treated as 
a cluster, with independence working correlation matrix used in the robust variance estimator.
This function is called in \code{geepack.lgst.int.batch.imputed} function to apply interaction test to all imputed SNPs in a 
genotype file. The interaction test is carried out by the \code{geese} function from package \code{geepack}.
}
\usage{
geepack.lgst.int.imputed(snp,phen,test.dat,covar,cov.int,sub="N")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{snp}{ genotype data of a SNP }
  \item{phen}{ a character string for a phenotype name in \code{test.dat} }
  \item{test.dat}{ the product of merging phenotype, genotype and pedigree data, should be ordered by "famid" }
  \item{covar}{ a character vector for covariates in \code{test.dat} }
  \item{cov.int}{ a character string naming the covariate for interaction, the covariate has to be included in \code{covar} }
  \item{sub}{ "N" (default) for no stratified analysis, and "Y" for requesting stratified analyses (only when \code{cov.int} is dichotomous) }
}
\details{
  Similar to the details for \code{geepack.lgst.int.batch} function but here the SNP data contains imputed genotypes (allele dosages) 
that are continuous and range from 0 to 2. 
}
\value{
  Please see value in \code{geepack.lgst.int.batch.imputed} function.
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