\name{lmeVpack.batch.imputed}
\alias{lmeVpack.batch.imputed}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ function to efficiently test associations between a continuous trait and a batch of imputed SNPs in families using 
Linear Mixed Effects model }
\description{
  A faster version of function \code{lmepack.batch.imputed}. Unlike \code{lmepack.batch.imputed} that estimates polygenic variation
for every SNP in a batch of imputed SNPs, \code{lmeVpack.batch.imputed} only estimates once for a batch of imputed SNPs. Particularly 
recommended for analyzing 1000G imputed genotype data.
}
\usage{
lmeVpack.batch.imputed(phenfile, genfile, pedfile, phen, kinmat, covars = NULL, 
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
  Similar to the details for \code{lmepack.batch.imputed}.
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
 coxme package: mixed-effects Cox models, sparse matrices, and modeling data from large pedigrees.
Beth Atkinson (atkinson@mayo.edu) for pedigree functions.Terry Therneau (therneau@mayo.edu) for all other functions.
2007. Ref Type: Computer Program http://cran.r-project.org/. 

 Abecasis, G. R., Cardon, L. R., Cookson, W. O., Sham, P. C., & Cherny, S. S. Association analysis in 
a variance components framework. \emph{Genet Epidemiol}, \bold{21} Suppl 1, S341-S346 (2001).
}
\author{ Qiong Yang <qyang@bu.edu> and Ming-Huei Chen <mhchen@bu.edu> }

\examples{
\dontrun{
lmeVpack.batch.imputed(phenfile="simphen.csv",genfile="simgen.csv",pedfile="simped.csv",
phen="SIMQT",kinmat="simkmat.Rdata",outfile="simout.csv",covars=c("age","sex"),
sep.ped=",",sep.phe=",",sep.gen=",")
}
}

