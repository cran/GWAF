\name{lmepack.batch}
\alias{lmepack.batch}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ function to test genetic associations between a continuous trait and a batch of genotyped SNPs in families using 
Linear Mixed Effects model  }
\description{
   Fit linear mixed effects (LME) model to test associations between a continuous phenotype 
and all SNPs in a genotype file in family data under user specified genetic model. The SNP genotype is treated
as a fixed effect, and a random effect correlated according to degree of relatedness within a 
family is also fitted. In each trait-SNP assocaition test, the \code{lmekin} function from
package \code{coxme} is used. 
}
\usage{
lmepack.batch(phenfile, genfile, pedfile, phen, kinmat, model = "a", covars = NULL, 
outfile, col.names = T, sep.ped = ",", sep.phe = ",", sep.gen = ",")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genfile}{ a character string naming the genotype file for reading (see format requirement in details) }
  \item{phenfile}{ a character string naming the phenotype file for reading (see format requirement in details) }
  \item{pedfile}{ a character string naming the pedigree file for reading (see format requirement in details) }
  \item{outfile}{ a character string naming the result file for writing }
  \item{phen}{ a character string for a phenotype name in \code{phenfile} }
  \item{covars}{ a character vector for covariates in \code{phenfile} }
  \item{model}{ a single character of 'a','d','g', or 'r', with 'a'=additive, 'd'=dominant, 'g'=general and 'r'=recessive models }
  \item{kinmat}{ a character string naming the file where kinship coefficient matrix is kept }
  \item{col.names}{ a logical value indicating whether the output file should contain column names }
  \item{sep.ped}{ the field separator character for pedigree file }
  \item{sep.phe}{ the field separator character for phenotype file }
  \item{sep.gen}{ the field separator character for genotype file }
}
\details{
   The \code{lmepack.batch} function first reads in and merges phenotype-covariates, genotype 
and pedigree files, then tests the association of \code{phen} against all SNPs in \code{genfile}.
\code{genfile} contains unique individual id and genotype data, with the column names being "id" and SNP names.
For each SNP, the genotype data should be coded as 0, 1, 2 indicating the numbers of the coded alleles. The SNP name in genotype file should not have any 
dash, '-' and other special characters(dots and underscores are OK). \code{phenfile} contains unique individual id, 
phenotype and covariates data, with the column names being "id" and phenotype and 
covaraite names. \code{pedfile} contains pedigree informaion, with the column names being 
"famid","id","fa","mo","sex". In all files, missing value should be an empty space, except missing parental id in \code{pedfile}.
 SNPs with low genotype counts (especially minor allele homozygote) may be omitted
or analyzed with dominant model. The \code{lmepack.batch} function fits LME model using a modified \code{lmekin} function from \code{coxme} package. 
}
\value{
  No value is returned. Instead, results are written to \code{outfile}.
When the genetic model is 'a', 'd' or 'r', the result includes the following columns.
When the genetic model is 'g', \code{beta} and \code{se} are replaced with \code{beta10},
\code{beta20}, \code{beta21}, \code{se10}, \code{se20}, \code{se21}.

  \item{phen }{phenotype name}
  \item{snp }{SNP name}
  \item{n0 }{the number of individuals with 0 copy of coded alleles}
  \item{n1 }{the number of individuals with 1 copy of coded alleles}
  \item{n2 }{the number of individuals with 2 copies of coded alleles}
  \item{h2q }{the portion of phenotypic variation explained by the SNP}
  \item{beta }{regression coefficient of SNP covariate}
  \item{se }{standard error of \code{beta}}
  \item{chisq }{Chi-square statistic for testing \code{beta} not equal to zero}
  \item{df }{degree of freedom of the Chi-square statistic}
  \item{model }{model actually used in the analysis}
  \item{pval }{p-value of the chi-square statistic}
  \item{  }{  }
  \item{beta10 }{regression coefficient of genotype with 1 copy of coded allele vs. that with 0 copy}
  \item{beta20 }{regression coefficient of genotype with 2 copy of coded allele vs. that with 0 copy}
  \item{beta21 }{regression coefficient of genotype with 2 copy of coded allele vs. that with 1 copy}
  \item{se10 }{standard error of \code{beta10}}
  \item{se20 }{standard error of \code{beta20}}
  \item{se21 }{standard error of \code{beta21}}
}
\references{ 
 coxme package: mixed-effects Cox models, sparse matrices, and modeling data from large pedigrees.
Beth Atkinson (atkinson@mayo.edu) for pedigree functions.Terry Therneau (therneau@mayo.edu) for all other functions.
2007. Ref Type: Computer Programhttp://cran.r-project.org/web/packages/coxme/. 

 Abecasis, G. R., Cardon, L. R., Cookson, W. O., Sham, P. C., & Cherny, S. S. Association analysis in 
a variance components framework. \emph{Genet Epidemiol}, \bold{21} Suppl 1, S341-S346 (2001).
}

\author{ Qiong Yang <qyang@bu.edu> and Ming-Huei Chen <mhchen@bu.edu> }

\examples{
\dontrun{
lmepack.batch(phenfile="simphen.csv",genfile="simgen.csv",pedfile="simped.csv",
phen="SIMQT",kinmat="simkmat.Rdata",model="a",outfile="simout.csv",col.names=T,
sep.ped=",",sep.phe=",",sep.gen=",")
}
}
