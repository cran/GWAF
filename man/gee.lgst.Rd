\name{gee.lgst}
\alias{gee.lgst}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ function for testing association between a dichotomous trait and a genotyped SNP in family data using GEE}
\description{
  Fit logistic regression via GEE to test association between a dichotomous phenotype 
and one genotyped SNP in a genotype file with user specified genetic model. Each family is treated as 
a cluster, with independence working correlation matrix used in the robust variance estimator.
This function is called in gee.lgst.batch function to apply association test to all SNPs in the 
genotype data.
}
\usage{
gee.lgst(snp, phen, test.dat, covar=NULL, model="a")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{snp}{ genotype data of a SNP }
  \item{phen}{ a character string for a phenotype name in phenfile }
  \item{test.dat}{ the product of merging phenotype, genotype and pedigree data, should be ordered by "famid" }
  \item{covar}{ a character vector for covariates in phenfile }
  \item{model}{ a single character of 'a','d','g', or 'r', with 'a'=additive, 'd'=dominant, 'g'=general and 'r'=recessive models }
}
\details{
     The \code{gee.lgst} function tests association between a dichtomous trait and a SNP from a dataset that contains phenotype, genotype and 
pedigree data, where the dataset needs to be ordered by famid.
}
\value{
  Please see output in \code{gee.lgst.batch}.
}
\references{   Liang, K.Y. and Zeger, S.L. (1986)
  Longitudinal data analysis using generalized linear models.
  \emph{Biometrika}, \bold{73} 13--22. 

  Zeger, S.L. and Liang, K.Y. (1986)
  Longitudinal data analysis for discrete and continuous outcomes.
  \emph{Biometrics}, \bold{42} 121--130.

  Vincent J Carey.Ported to R by Thomas Lumley (versions 3.13 and 4.4) and Brian Ripley (version 4.13). gee: Generalized Estimation Equation solver. 
  [4.13]. 2007. Ref Type: Computer Program, http://cran.r-project.org/
}
\author{ Qiong Yang <qyang@bu.edu> and Ming-Huei Chen <mhchen@bu.edu>  }
\seealso{ \code{gee()} function from package \code{gee} }
