\name{geepack.lgst}
\alias{geepack.lgst}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ function for testing association between a dichotomous trait and a genotyped SNP in family data using Generalized Estimation Equation model}
\description{
  Fit logistic regression via GEE to test association between a dichotomous phenotype 
and one genotyped SNP in a genotype file with user specified genetic model. Each family is treated as 
a cluster, with independence working correlation matrix used in the robust variance estimator. The trait-SNP association 
test is carried out by the \code{geese} function from package \code{geepack}. 
This function is called in \code{geepack.lgst.batch} function to apply association test to all SNPs in the 
genotype data.
}
\usage{
geepack.lgst(snp, phen, test.dat, covar = NULL, model = "a")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{snp}{ genotype data of a SNP }
  \item{phen}{ a character string for a phenotype name in \code{test.dat} }
  \item{test.dat}{ the product of merging phenotype, genotype and pedigree data, should be ordered by "famid" }
  \item{covar}{ a character vector for covariates in \code{test.dat} }
  \item{model}{ a single character of 'a','d','g', or 'r', with 'a'=additive, 'd'=dominant, 'g'=general and 'r'=recessive models }
}
\details{
     The \code{geepack.lgst} function tests association between a dichtomous trait and a SNP from a dataset that contains phenotype, genotype and 
pedigree data (test.dat), where the dataset needs to be ordered by famid.
}
\value{
  Please see output in \code{geepack.lgst.batch}. 
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
\examples{
\dontrun{
geepack.lgst(snp=data[,"rs123"],phen="CVD",test.dat=data,model="a",covar=c("age","sex"))
}
}