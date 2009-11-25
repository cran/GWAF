\name{gee.lgst.imputed}
\alias{gee.lgst.imputed}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ function for testing association between a dichotomous trait and an imputed SNP in family data using GEE}
\description{
  Fit logistic regression via GEE to test association between a dichotomous phenotype 
and one imputed SNP in a genotype file. Each family is treated as 
a cluster, with independence working correlation matrix used in the robust variance estimator.
This function is called in \code{gee.lgst.batch.imputed} function to apply association test to all imputed SNPs in the 
genotype data.
}
\usage{
gee.lgst.imputed(snp, phen, test.dat, covar = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{snp}{ imputed genotype data of a SNP}
  \item{phen}{ a character string for a phenotype name in phenfile }
  \item{test.dat}{ the product of merging phenotype, genotype and pedigree data, should be ordered by "famid" }
  \item{covar}{ a character vector for covariates in phenfile }
}
\details{
  Similar to the details for \code{gee.lgst} function but here the SNP data contains imputed genotypes (allele dosages) 
that are continuous and range from 0 to 2. In addition, the user specified genetic model argument is not available.
}
\value{
  Please see output in \code{gee.lgst.batch.imputed}.
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