\name{GWplot}
\alias{GWplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ function for making genome-wide p-values plot }
\description{
  GWplot function plots -log_{10} p-value based on SNP's chromosomal position in bitmap format.  
}
\usage{
GWplot(data, pval, pos, chr, chr.plot=c(1:22,"X"), title.text="", 
ylim=Inf, outfile, cutoff1=5e-08, cutoff2=4e-07)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ a dataframe that contains p-values, chromosome number and physical position of SNPs }
  \item{pval}{ a character string correspong to the name of the p-value column }
  \item{pos}{ a character string correspong to the name of column with SNP physical positions  }
  \item{chr}{ a character string correspong to the name of column with SNP chromosome number }
  \item{chr.plot}{ the chsomosomes of interest for GWplot; either 1:22 or c(1:22,"X"), default chr.plot=c(1:22,"X"), "X" for X chromosome }
  \item{title.text}{ the title of the genome-wide p-value plot }
  \item{ylim}{ the maximum of -log_{10} p-value to be plotted, useful when not want to plot extremely small p-values }
  \item{outfile}{ the file name (xxxx.bmp) for output plot in bitmap format }
  \item{cutoff1}{ genome-wide significance; default is 5E-8 ; p-values below this threshold will be highlighted in red}
  \item{cutoff2}{ suggestive genome-wide significance; default is 4E-7; p-values below this threshold but above cutoff1 will be highlighted in blue }
}
\details{
  When the dataset has 0 p-value, \code{GWplot} will generate pvalzero.csv that contain the results with 0 p-value and make the genome-wide p-values
plot by replacing 0 p-value with 5E-324. P-values that reach genome-wide significance are displayed in red color; P-values that reach 
suggestive genome-wide significance but not genome-wide significance are displayed in blue color. 
}
\author{ Qiong Yang <qyang@bu.edu> and Ming-Huei Chen <mhchen@bu.edu> }
