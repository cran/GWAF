\name{geepack.lgst.imputed}
\alias{geepack.lgst.imputed}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ function for testing association between a dichotomous trait and an imputed SNP in family data using Generalized Estimation Equation model}
\description{
  Fit logistic regression via Generalized Estimation Equation (GEE) to test association between a dichotomous phenotype 
and one imputed SNP in a genotype file in family data under additive genetic model. Each family is treated as 
a cluster, with independence working correlation matrix used in the robust variance estimator. The trait-SNP association 
test is carried out by the \code{geese} function from package \code{geepack}. 
This function is called in \code{geepack.lgst.batch.imputed} function to apply association test to all imputed SNPs in a 
genotype file.
}
\usage{
geepack.lgst.imputed(snp, phen, test.dat, covar = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{snp}{ imputed genotype data of a SNP}
  \item{phen}{ a character string for a phenotype name in \code{test.dat}  }
  \item{test.dat}{ the product of merging phenotype, genotype and pedigree data, should be ordered by "famid" }
  \item{covar}{ a character vector for covariates in \code{test.dat} }
}
\details{
  Similar to the details for \code{geepack.lgst} function but here the SNP data contains imputed genotypes (allele dosages) 
that are continuous and range from 0 to 2. In addition, the user-specified genetic model argument is not available.
}
\value{
  Please see output in \code{geepack.lgst.batch.imputed}.
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