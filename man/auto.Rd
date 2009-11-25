\name{auto}
\alias{auto}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ function to generate scripts for genome-wide association analysis using GWAF}
\description{
  Given a phenotype file, a pedigree file, a phenotype of interest, covariates, analysis of interest
(can be 'lme', 'lme.imputed', 'gee' or 'gee.imputed'), a path or directory (\code{genopath}) that keeps geotype files, 
and other arguments, \code{auto} function generates one R script and one shell script that excutes the
R script for each genotype file, and one list file that excutes all shell scripts in batch mode. 
}
\usage{
auto(genopath, phenfile, pedfile, outfile, phen, covars = NULL, analysis, lib.loc, model = NULL, 
kinmat = NULL, col.names = F, sep.ped = ",", sep.phe = ",", sep.gen = ",")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genopath}{ a character string indicating the path or directory that keeps genotype files to be analyzed }
  \item{phenfile}{ a character string naming the phenotype file for reading (see format requirement in details) }
  \item{pedfile}{ a character string naming the pedigree file for reading (see format requirement in details) }
  \item{outfile}{ a character string naming the result file for writing }
  \item{phen}{ a character string for a phenotype name in phenfile }
  \item{covars}{ a character vector for covariates in phenfile }
  \item{analysis}{ a character string indicating the analysis of interest available in \code{GWAF} package, can be 'lme', 'lme.imputed', 'gee' or 'gee.imputed' }
  \item{lib.loc}{ a character string indicating the location of \code{GWAF} package }
  \item{model}{ a single character of 'a','d','g', or 'r', with 'a'=additive, 'd'=dominant, 'g'=general and 'r'=recessive models; Not appropriate/needed for analyzing imputed SNPs }
  \item{kinmat}{ a character string naming the file where kinship coefficient matrix is kept; needed for LME analyses }
  \item{col.names}{ a logical value indicating whether the output file should contain column names }
  \item{sep.ped}{ the field separator character for pedigree file }
  \item{sep.phe}{ the field separator character for phenotype file }
  \item{sep.gen}{ the field separator character for genotype file }
}
\details{
  \code{auto} function generates one R script and one shell script that excutes the R script for each genotype file, and one list file 
that can excute all shell scripts in batch mode to analyze all genotype data in \code{genopath}. These scripts are named based on the phenotype of interest,
the analysis of interest and the time these scripts are generated. After creating these scripts, \code{auto} 
function gives a message telling the user how to submit ALL the jobs (using ksh XXXX.lst). When a submitted job
is completed, a log file indicating which genotype file was analyzed will be generated and the R script and the 
shell script will be removed. The number of log files should equal to the number of genotype files in \code{genopath}, if all jobs 
are completed. All the results will be written and appended to the user specified single output file. Different
\code{outfile} should be assigned for different \code{genopath} to avoid over-writting.
}
\value{
    No value is returned. Instead, results are written to \code{outfile}.
}
\author{ Ming-Huei Chen <mhchen@bu.edu> and Qiong Yang <qyang@bu.edu>}