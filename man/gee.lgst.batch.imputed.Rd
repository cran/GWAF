\name{gee.lgst.batch.imputed}
\alias{gee.lgst.batch.imputed}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ function to test genetic association between a dichotomous trait and a batch of imputed SNPs in families using 
GEE }
\description{
  Fit logistic regression via GEE to test association between a dichotomous phenotype 
and all imputed SNPs in a genotype file. Each family is treated as 
a cluster, with independence working correlation matrix used in the robust variance estimator.
This function applies the same trait-SNP association test to all SNPs in the imputed genotype data. 
The trait-SNP association test is carried out by gee.lgst function where the 
the \code{gee()} function from package \code{gee} is used.
}
\usage{
gee.lgst.batch.imputed(genfile, phenfile, pedfile, outfile, phen, covars = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genfile}{ a character string naming the genotype file for reading (see format requirement in details) }
  \item{phenfile}{ a character string naming the phenotype file for reading (see format requirement in details) }
  \item{pedfile}{ a character string naming the pedigree file for reading (see format requirement in details) }
  \item{outfile}{ a character string naming the result file for writing }
  \item{phen}{ a character string for a phenotype name in phenfile }
  \item{covars}{ a character vector for covariates in phenfile }
}
\details{
  Similar to the details for 'gee.lgst.batch' but here the SNP data contains imputed genotypes (allele dosages) 
that are continuous and range from 0 to 2. In addition, the user 
specified genetic model argument is not available.
}
\value{
  No value is returned. Instead, results are written to \code{outfile}.

  \item{phen }{phenotype name}
  \item{snp }{SNP name}
  \item{N }{the number of individuals in analysis}
  \item{Nd }{the number of individuals in affected sample in analysis}
  \item{AF }{imputed allele frequency of coded allele}
  \item{AFd }{imputed allele frequency of coded allele in affected sample}
  \item{beta }{regression coefficient of SNP covariate}
  \item{se }{standard error of \code{beta}}
  \item{remark }{warning or additional information for the analysis, note that the genotype counts are based
on rounded imputed genotypes; 'not converged' indicates the 
GEE analysis did not converge; 'logistic reg' indicates GEE model is replaced by logistic regression;
'exp count<5' indicates any expected count is less than 5 in phenotype-genotype table; 'not converged 
and exp count<5', 'logistic reg & exp count<5' are noted similarly; 'collinearity' indicates collinearity
exists between SNP and some covariates}
  \item{pval }{p-value of the chi-square statistic}
}

\author{ Qiong Yang <qyang@bu.edu> and Ming-Huei Chen <mhchen@bu.edu>  }