\name{glmm.lgst}
\alias{glmm.lgst}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ function for testing association between a dichotomous trait and a genotyped SNP in family data using Generalized Linear Mixed Effects model}
\description{
  Fit a Generalized Linear Mixed Effects model (GLMM) with logistic link and a normal distributed random intercept for each cluster  to test association between a dichotomous phenotype 
and one genotyped SNP in a genotype file with user specified genetic model. Each family is treated as 
a cluster. When analyze rare variants for dichotomous traits, GLMM,  as implemeted by this function, is recommended over other methods such as GEE. The trait-SNP association 
test is carried out by the \code{lmer} function from package \code{lme4}. 
This function is called in \code{glmm.lgst.batch} function to apply association test to all SNPs in the 
genotype data.
}
\usage{
glmm.lgst(snp, phen, test.dat, covar = NULL, model = "a")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{snp}{ genotype data of a SNP }
  \item{phen}{ a character string for a phenotype name in \code{test.dat} }
  \item{test.dat}{ the product of merging phenotype, genotype and pedigree data, should be ordered by "famid" }
  \item{covar}{ a character vector for covariates in \code{test.dat} }
  \item{model}{ a single character of 'a','d','g', or 'r', with 'a'=additive, 'd'=dominant, 'g'=general and 'r'=recessive models }
}
\details{
     The \code{glmm.lgst} function tests association between a dichtomous trait and a SNP from a dataset that contains phenotype, genotype and 
pedigree data (test.dat), where the dataset needs to be ordered by famid.
}
\value{
  Please see output in \code{glmm.lgst.batch}. 
}
\references{   
Bates D, Maechler M and Bolker B. (2010) lme4: Linear mixed-effects models using S4 classes
Longitudinal data analysis using generalized linear models. http://cran.r-project.org/web/packages/lme4/. 

Chen MH, Liu X, Wei F, Larson MG, Fox CS, Vasan RS and Yang Q. (2011) A comparison of strategies for analyzing 
dichotomous outcomes in genome-wide association studies with general pedigrees. Genetic Epidemiology 35:650-657. 

}
\author{ Qiong Yang <qyang@bu.edu> and Ming-Huei Chen <mhchen@bu.edu>  }
\seealso{ \code{lmer} function from package \code{lme4} }
\examples{
\dontrun{
glmm.lgst(snp=data[,"rs123"],phen="CVD",test.dat=data,model="a",covar=c("age","sex"))
}
}