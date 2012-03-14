\name{lmepack.batch.imputed}
\alias{lmepack.batch.imputed}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ function to test associations between a continuous trait and a batch of imputed SNPs in families using 
Linear Mixed Effects model }
\description{
  Fit linear mixed effects (LME) model to test associations between a continuous phenotype 
and all imputed SNPs in a genotype file in family data under additive genetic model. The SNP genotype is treated
as a fixed effect, and a random effect correlated according to degree of relatedness within a 
family is also fitted. In each trait-SNP assocaition test, the \code{lmekin} function which is modified from
the same named function in package \code{kinship} is used. As \code{kinship} package has been archived, \code{GWAF} package now
includes functions from \code{kinship} package to support the modified \code{lmekin} function. Archived \code{kinship} package can be found 
at http://cran.r-project.org/web/packages/kinship/index.html.
}
\usage{
lmepack.batch.imputed(phenfile, genfile, pedfile, phen, kinmat, covars = NULL, 
outfile, col.names = T, sep.ped = ",", sep.phe = ",", sep.gen = ",")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phenfile}{ a character string naming the phenotype file for reading (see format requirement in details) }
  \item{genfile}{ a character string naming the genotype file for reading (see format requirement in details) }
  \item{pedfile}{ a character string naming the pedigree file for reading (see format requirement in details)}
  \item{phen}{ a character string for a phenotype name in \code{phenfile}  }
  \item{kinmat}{ a character string naming the file where kinship coefficient matrix is kept }
  \item{covars}{ a character vector for covariates in \code{phenfile} }
  \item{outfile}{ a character string naming the result file for writing }
  \item{col.names}{ a logical value indicating whether the output file should contain column names }
  \item{sep.ped}{ the field separator character for pedigree file }
  \item{sep.phe}{ the field separator character for phenotype file }
  \item{sep.gen}{ the field separator character for genotype file }
}
\details{
  Similar to the details for \code{lmepack.batch} function but here the SNP data contains imputed genotypes (allele dosages) 
that are continuous and range from 0 to 2. In addition, the user specified genetic model argument is not available.
}
\value{
   No value is returned. Instead, results are written to \code{outfile}.

  \item{phen }{phenotype name}
  \item{snp }{SNP name}
  \item{N }{the number of individuals in analysis}
  \item{AF }{imputed allele frequency of coded allele}
  \item{h2q }{the portion of phenotypic variation explained by the SNP}
  \item{beta }{regression coefficient of SNP covariate}
  \item{se }{standard error of \code{beta}}
  \item{pval }{p-value of testing \code{beta} not equal to zero}
 
}
\references{ 
 kinship package: mixed-effects Cox models, sparse matrices, and modeling data from large pedigrees.
Beth Atkinson (atkinson@mayo.edu) for pedigree functions.Terry Therneau (therneau@mayo.edu) for all other functions.
2007. Ref Type: Computer Program http://cran.r-project.org/. 

 Abecasis, G. R., Cardon, L. R., Cookson, W. O., Sham, P. C., & Cherny, S. S. Association analysis in 
a variance components framework. \emph{Genet Epidemiol}, \bold{21} Suppl 1, S341-S346 (2001).
}
\author{ Qiong Yang <qyang@bu.edu> and Ming-Huei Chen <mhchen@bu.edu> }