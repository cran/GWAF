\name{geepack.lgst.batch}
\alias{geepack.lgst.batch}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ function to test genetic associations between a dichotomous trait and a batch of genotyped SNPs in families using 
Generalized Estimation Equation model }
\description{
   Fit logistic regression via Generalized Estimation Equation (GEE) to test associations between a dichotomous phenotype 
and all genotyped SNPs in a genotype file in family data with user specified genetic model. Each pedigree is treated as 
a cluster, with independence working correlation matrix used in the robust variance estimator.
This function applies the same trait-SNP association test to all SNPs in the genotype data. 
The trait-SNP association test is carried out by \code{geepack.lgst} function where the 
the \code{geese} function from package \code{geepack} is used.
}
\usage{
geepack.lgst.batch(genfile, phenfile, pedfile, outfile, phen, covars = NULL, 
model = "a", col.names = T, sep.ped = ",", sep.phe = ",", sep.gen = ",")
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
  \item{col.names}{ a logical value indicating whether the output file should contain column names }
  \item{sep.ped}{ the field separator character for pedigree file }
  \item{sep.phe}{ the field separator character for phenotype file }
  \item{sep.gen}{ the field separator character for genotype file }
}
\details{
   The \code{geepack.lgst.batch} function first reads in and merges phenotype-covariates, genotype 
and pedigree files, then tests the association of \code{phen} against all SNPs in \code{genfile}.
\code{genfile} contains unique individual id and genotype data, with the column names being "id" and SNP names.
For each genotyped SNP, the genotype data should be coded as 0, 1, 2 indicating the numbers of the coded alleles. The SNP names in genotype file should not have any 
dash, '-' and other special characters(dots and underscores are OK). \code{phenfile} contains unique individual id, 
phenotype and covariates data, with the column names being "id" and phenotype and 
covaraite names. \code{pedfile} contains pedigree informaion, with the column names being 
"famid","id","fa","mo","sex". In all files, missing value should be an empty space, except missing parental id in \code{pedfile}.
 Only phenotypes with two categories are analyzed. A phenotype should be coded as 
0 and 1, with 1 denoting affected and 0 unaffected. SNPs with low genotype counts 
(especially minor allele homozygote) may be omitted or analyzed with dominant model or
analyzed with logistic regression. 
The \code{geepack.lgst.batch} function fits GEE model using each pedigree as a cluster 
 with \code{geepack.lgst} function from GWAF package and \code{geese} function from \code{geepack} package. 
}
\value{
  No value is returned. Instead, results are written to \code{outfile}.
When the genetic model is 'a', 'd' or 'r', the result includes the following columns.
When the genetic model is 'g', \code{beta} and \code{se} are replaced with \code{beta10},
\code{beta20}, \code{beta21}, \code{se10}, \code{se20}, and \code{se21} .

  \item{phen }{phenotype name}
  \item{snp }{SNP name}
  \item{n0 }{the number of individuals with 0 copy of coded alleles}
  \item{n1 }{the number of individuals with 1 copy of coded alleles}
  \item{n2 }{the number of individuals with 2 copies of coded alleles}
  \item{nd0 }{the number of individuals with 0 copy of coded alleles in affected sample}
  \item{nd1 }{the number of individuals with 1 copy of coded alleles in affected sample}
  \item{nd2 }{the number of individuals with 2 copies of coded alleles in affected sample}
  \item{miss.0 }{Genotype missing rate in unaffected sample }
  \item{miss.1 }{Genotype missing rate in affected sample }
  \item{miss.diff.p }{P-value of differential missingness test between unaffected and affected samples }  
  \item{beta }{regression coefficient of SNP covariate}
  \item{se }{standard error of \code{beta}}
  \item{chisq }{Chi-square statistic for testing \code{beta} not equal to zero}
  \item{df }{degree of freedom of the Chi-square statistic}
  \item{model }{model actually used in the analysis}
  \item{remark }{warning or additional information for the analysis, 'not converged' indicates the 
GEE analysis did not converge; 'logistic reg' indicates GEE model is replaced by logistic regression;
'exp count<5' indicates any expected count is less than 5 in phenotype-genotype table; 'not converged 
and exp count<5', 'logistic reg & exp count<5' are noted similarly; 'collinearity' indicates collinearity
exists between SNP and some covariates}
  \item{pval }{p-value of the chi-square statistic}
  \item{  }{  }

  \item{beta10 }{regression coefficient of genotype with 1 copy of coded allele vs. that with 0 copy}
  \item{beta20 }{regression coefficient of genotype with 2 copy of coded allele vs. that with 0 copy}
  \item{beta21 }{regression coefficient of genotype with 2 copy of coded allele vs. that with 1 copy}
  \item{se10 }{standard error of \code{beta10}}
  \item{se20 }{standard error of \code{beta20}}
  \item{se21 }{standard error of \code{beta21}}
}
\author{ Qiong Yang <qyang@bu.edu> and Ming-Huei Chen <mhchen@bu.edu>  }

\examples{
\dontrun{
geepack.lgst.batch(phenfile="simphen.csv",genfile="simgen.csv",pedfile="simped.csv",
phen="SIMQT",model="a",outfile="simout.csv",sep.ped=",",sep.phe=",",sep.gen=",")
}
}
