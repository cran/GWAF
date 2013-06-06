\name{qq}
\alias{qq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ function to make Qantile-Qantile (QQ) plot for p-values}
\description{
  qq function makes the QQ plot of p-values against a uniform (0,1) distribution. 
The genomic control parameter for one degree freedom chi-square statistics corresponding to the p-values is also 
plotted.
}
\usage{
qq(pvalue, outfile)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pvalue}{ P-values of interest. }
  \item{outfile}{ the file name (xxxx.bmp) for output QQ plot in bitmap format }
}
\author{ Qiong Yang <qyang@bu.edu> and Ming-Huei Chen <mhchen@bu.edu> }
\examples{
\dontrun{
qq(pvalue=result[,"pval"],outfile="QQ_test.bmp")
}
}