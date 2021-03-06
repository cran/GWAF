\name{geepack.lgst.int.batch}
\alias{geepack.lgst.int.batch}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ function to test gene-environment or gene-gene interactions between a dichotomous trait and a batch of genotyped SNPs in families using 
Generalized Estimation Equation model }
\description{
   Fit logistic regression via Generalized Estimation Equation (GEE) to test gene-environment or gene-gene interactions for a dichotomous phenotype 
and all genotyped SNPs in a genotype file in family data under additive genetic model. The interaction term is the product of SNP genotype and a covariate for interaction (\code{cov.int}). 
The covariate for interaction (\code{cov.int}) can be SNP genotype (gene-gene interaction) or an environmental factor (gene-environment interaction). Only one
interaction term is allowed. When \code{cov.int} is dichotomous, stratified analyses can be requested by specifying \code{sub}="Y". The covariance between the main 
effect (SNP) and the interaction effect is provided in the output when stratified analysis is not requested. Each pedigree is treated as 
a cluster with independence working correlation matrix used in the robust variance estimator.
This function applies the same interaction test to all SNPs in a genotype file. 
The interaction test is carried out by \code{geepack.lgst.int} function from \code{GWAF} where the 
the \code{geese} function from package \code{geepack} is used.
}
\usage{
geepack.lgst.int.batch(genfile,phenfile,pedfile,outfile,phen,covars,cov.int,sub="N",
col.names=T,sep.ped=",",sep.phe=",",sep.gen=",")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genfile}{ a character string naming the genotype file for reading (see format requirement in details) }
  \item{phenfile}{ a character string naming the phenotype file for reading (see format requirement in details) }
  \item{pedfile}{ a character string naming the pedigree file for reading (see format requirement in details) }
  \item{outfile}{ a character string naming the result file for writing }
  \item{phen}{ a character string for a phenotype name in \code{phenfile} }
  \item{covars}{ a character vector for covariates in \code{phenfile} }
  \item{cov.int}{ a character string naming the covariate for interaction, the covariate has to be included in \code{covars} }
  \item{sub}{ "N" (default) for no stratified analysis, and "Y" for requesting stratified analyses (only when \code{cov.int} is dichotomous) }
  \item{col.names}{ a logical value indicating whether the output file should contain column names }
  \item{sep.ped}{ the field separator character for pedigree file }
  \item{sep.phe}{ the field separator character for phenotype file }
  \item{sep.gen}{ the field separator character for genotype file }
}
\details{
   The \code{geepack.lgst.int.batch} function first reads in and merges phenotype-covariates, genotype 
and pedigree files, then tests gene-environment or gene-gen interaction for \code{phen} against all SNPs in \code{genfile}.
Only one interaction term is allowed, so is the covariate for interaction (\code{cov.int}). When \code{cov.int} is dichotomous, stratified analyses can be 
requested by specifying \code{sub}="Y". The covariance between the main effect (SNP) and the interaction effect is 
provided in the output when stratified analysis is not requested.
\code{genfile} contains unique individual id and genotype data, with the column names being "id" and SNP names.
For each genotyped SNP, the genotype data should be coded as 0, 1, 2 indicating the numbers of the coded alleles. The SNP names in genotype file should not have any 
dash, '-' and other special characters(dots and underscores are OK). \code{phenfile} contains unique individual id, 
phenotype and covariates data, with the column names being "id" and phenotype and 
covaraite names. \code{pedfile} contains pedigree informaion, with the column names being 
"famid","id","fa","mo","sex". In all files, missing value should be an empty space, except missing parental id in \code{pedfile}.
 Only phenotypes with two categories are analyzed. A phenotype should be coded as 
0 and 1, with 1 denoting affected and 0 unaffected. SNPs with low genotype counts 
(especially minor allele homozygote) may be omitted or analyzed with logistic regression. 
The \code{geepack.lgst.int.batch} function fits GEE model using each pedigree as a cluster 
 with \code{geepack.lgst.int} function from \code{GWAF} package and \code{geese} function from \code{geepack} package. 
}
\value{
  No value is returned. Instead, results are written to \code{outfile}.
If stratified analyses are requested, the result file will include the following columns. Otherwise, \code{cov_beta_snp_beta_int} will be included instead of 
the results from stratified analyses, that is, \code{beta_snp_cov0}, \code{se_snp_cov0}, \code{pval_snp_cov0}, \code{beta_snp_cov1}, \code{se_snp_cov1}, 
and \code{pval_snp_cov1}.

  \item{phen }{phenotype name}
  \item{snp }{SNP name}
  \item{covar_int }{the covariate for interaction }
  \item{n }{sample size used in analysis}
  \item{AF }{allele frequency of the coded allele}
  \item{nd }{the number of individuals in affected sample}
  \item{AFd }{allele frequency of the coded allele in affected sample}
  \item{model }{genetic model used in analysis, additive model only}
  \item{beta_snp }{regression coefficient of SNP covariate}
  \item{se_snp }{standard error of \code{beta_snp}}
  \item{pval_snp }{p-value of testing \code{beta_snp} not equal to zero}
  \item{beta_snp_cov0 }{regression coefficient of SNP covariate in stratified analysis using the subset where \code{cov.int} level is 0}
  \item{se_snp_cov0 }{standard error of \code{beta_snp_cov0}}
  \item{pval_snp_cov0 }{p-value of testing \code{beta_snp_cov0} not equal to zero}
  \item{beta_snp_cov1 }{regression coefficient of SNP covariate in stratified analysis using the subset where \code{cov.int} level is 1}
  \item{se_snp_cov1 }{standard error of \code{beta_snp_cov1}}
  \item{pval_snp_cov1 }{p-value of testing \code{beta_snp_cov1} not equal to zero}
  \item{beta_int }{regression coefficient of the interaction term}
  \item{se_int }{standard error of \code{beta_int}}
  \item{pval_int }{p-value of testing \code{beta_int} not equal to zero}
  \item{remark }{warning or additional information for the analysis, 'not converged' indicates the 
GEE analysis did not converge; 'logistic reg' indicates GEE model is replaced by logistic regression;
'exp count<5' indicates any expected count is less than 5 in phenotype-genotype table; 'not converged 
and exp count<5', 'logistic reg & exp count<5' are noted similarly; 'collinearity' indicates collinearity
exists between SNP and some covariates}

}
\author{ Qiong Yang <qyang@bu.edu> and Ming-Huei Chen <mhchen@bu.edu>  }

\examples{
\dontrun{
geepack.lgst.int.batch(phenfile="simphen.csv",genfile="simgen.csv",pedfile="simped.csv",
phen="CVD",outfile="simout.csv",covars=c("sex","age"),cov.int="age",
sep.ped=",",sep.phe=",",sep.gen=",")
}
}
