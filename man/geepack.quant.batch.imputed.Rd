\name{geepack.quant.batch.imputed}
\alias{geepack.quant.batch.imputed}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ function to test associations between a continuous trait and a batch of imputed SNPs in families using 
Generalized Estimation Equation model }
\description{
  Fit Generalized Estimation Equation (GEE) model to test association between a continuous phenotype 
and all imputed SNPs in a genotype file in family data under additive genetic model. Each family is treated as a cluster, with independence working correlation matrix used in the 
robust variance estimator. The proportion of phenotype variation explained by the tested SNP is not provided. 
This function applies the same trait-SNP association test to all imputed SNPs in the genotype data. 
The trait-SNP association test is carried out by using the \code{geese} function from package \code{geepack}.
}
\usage{
geepack.quant.batch.imputed(phenfile, genfile, pedfile, phen, covars = NULL, 
outfile, col.names = T, sep.ped = ",", sep.phe = ",", sep.gen = ",")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phenfile}{ a character string naming the phenotype file for reading (see format requirement in details) }
  \item{genfile}{ a character string naming the genotype file for reading (see format requirement in details) }
  \item{pedfile}{ a character string naming the pedigree file for reading (see format requirement in details)}
  \item{phen}{ a character string for a phenotype name in \code{phenfile}  }
  \item{covars}{ a character vector for covariates in \code{phenfile} }
  \item{outfile}{ a character string naming the result file for writing }
  \item{col.names}{ a logical value indicating whether the output file should contain column names }
  \item{sep.ped}{ the field separator character for pedigree file }
  \item{sep.phe}{ the field separator character for phenotype file }
  \item{sep.gen}{ the field separator character for genotype file }
}
\details{
  Similar to the details for \code{geepack.quant.batch} function but here the SNP data contains imputed genotypes (allele dosages) 
that are continuous and range from 0 to 2. In addition, the user specified genetic model argument is not available.
}
\value{
   No value is returned. Instead, results are written to \code{outfile}.

  \item{phen }{phenotype name}
  \item{snp }{SNP name}
  \item{N }{the number of individuals in analysis}
  \item{AF }{imputed allele frequency of coded allele}
  \item{beta }{regression coefficient of SNP covariate}
  \item{se }{standard error of \code{beta}}
  \item{pval }{p-value of testing \code{beta} not equal to zero}
 
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