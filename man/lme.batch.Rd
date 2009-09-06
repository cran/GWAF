\name{lme.batch}
\alias{lme.batch}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ function to test genetic association for a continuous trait in families using 
Linear Mixed Effects model  }
\description{
   Fit linear mixed effect model to test association between a continuous phenotype 
and all SNPs in a genotype file with user specified genetic model. The SNP genotype is treated
as fixed effects, and a random effect correlated according to degree of relatedness within a 
family is also fitted.
In each trait-SNP assocaition test, the \code{lmekin()} function which is modified from
the same named function in package \code{kinship} is used. 
}
\usage{
lme.batch(phenfile, genfile, pedfile, phen, kinmat, model = "a", 
covars = NULL, outfile)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genfile}{ a character string naming the genotype file for reading(see format requirement in details) }
  \item{phenfile}{ a character string naming the phenotype file for reading (see format requirement in details) }
  \item{pedfile}{ a character string naming the pedigree file for reading(see format requirement in details) }
  \item{outfile}{ a character string naming the result file for writing }
  \item{phen}{ a character string for a phenotype name in phenfile }
  \item{covars}{ a character vector for covariates in phenfile }
  \item{model}{ a single character of 'a','d','g', or 'r', with 'a'=additive, 'd'=dominant, 'g'=general and 'r'=recessive models }
  \item{kinmat}{ a character string naming the file where kinship coefficient matrix is kept }
}                                                                                              
\details{
   The \code{lme.batch} function first reads in and merges comma-delimited phenotype-covariates, genotype 
and pedigree files, then tests the association of \code{phen} against all SNPs in \code{genfile}.
\code{genfile} is a comma delimited file, with the column names being "id" and SNP names separated by comma.
For each SNP, the genotype data should be coded as 0, 1, 2 indicating the numbers of less frequent 
alleles. The SNP name in genotype file should not have any 
dash, '-' and other special characters(dots and underscores are OK). \code{phenfile} is a comma delimited file, with the column names being "id" and phenotype and 
covaraite names separated by comma. \code{pedfile} is a comma delimited file, with the column names being 
"famid","id","fa","mo","sex". In all files, missing value should be an empty space. SNPs with low genotype counts (especially minor allele homozygote) may be omitted
or analyzed with dominant model. The 'lme.batch' function fits linear mixed effects (LME) model using  
with 'lme.batch' function from GWAF package and a modified 'lmekin' function from \code{kinship} package. 
}
\value{
  No value is returned. Instead, results are written to \code{outfile}.
When the genetic model is 'a', 'd' or 'r', the result includes the following columns.
When the genetic model is 'g', \code{beta} and \code{se} are replaced with \code{beta10},
\code{beta20},\code{beta21},\code{se10},\code{se20},\code{se21} .

  \item{phen }{phenotype name}
  \item{snp }{SNP name}
  \item{n0 }{the number of individuals with 0 copy of minor alleles}
  \item{n1 }{the number of individuals with 1 copy of minor alleles}
  \item{n2 }{the number of individuals with 2 copies of minor alleles}
  \item{h2q }{the portion of phenotypic variation explained by the SNP}
  \item{beta }{regression coefficient of SNP covariate}
  \item{se }{standard error of \code{beta}}
  \item{chisq }{Chi-square statistic for testing \code{beta} not equal to zero}
  \item{df }{degree of freedom of the Chi-square statistic}
  \item{model }{model actually used in the analysis}
  \item{pval }{p-value of the chi-square statistic}
  \item{  }{  }
  \item{beta10 }{regression coefficient of genotype with 1 copy of minor allele vs. that with 0 copy}
  \item{beta20 }{regression coefficient of genotype with 2 copy of minor allele vs. that with 0 copy}
  \item{beta21 }{regression coefficient of genotype with 2 copy of minor allele vs. that with 1 copy}
  \item{se10 }{standard error of \code{beta10}}
  \item{se20 }{standard error of \code{beta20}}
  \item{se21 }{standard error of \code{beta21}}
}
\references{ 
 kinship package: mixed-effects Cox models, sparse matrices, and modeling data from large pedigrees.
Beth Atkinson (atkinson@mayo.edu) for pedigree functions.Terry Therneau (therneau@mayo.edu) for all other functions.
2007. Ref Type: Computer Program http://cran.r-project.org/. 

 Abecasis, G. R., Cardon, L. R., Cookson, W. O., Sham, P. C., & Cherny, S. S. Association analysis in 
a variance components framework. \emph{Genet Epidemiol}, \bold{21} Suppl 1, S341-S346 (2001).
}
\author{ Qiong Yang <qyang@bu.edu> and Ming-Huei Chen <mhchen@bu.edu> }
