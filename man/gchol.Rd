\name{gchol}
\alias{gchol}
\alias{gchol.bdsmatrix}
\title{
  Generalized Cholesky decompostion
}
\description{
Perform the generalized Cholesky decompostion of a real symmetric matrix. The function and R documentation are included from the 
archived \code{kinship} package that can be found at http://cran.r-project.org/web/packages/kinship/index.html.
}
\usage{
gchol(x, tolerance=1e-10)
}
\arguments{
\item{x}{
the symmetric matrix to be factored
}
\item{tolerance}{
the numeric tolerance for detection of singular columns in x.
}
}
\value{
an object of class \code{gchol} containing the
generalized Cholesky decompostion.
It has the appearance of a lower triangular matrix.
}
\details{
A symmetric matrix A can be decomposed as LDL', where L is a lower
triangular matrix with 1's on the diagonal, L' is the transpose of
L, and D is diagonal.
The inverse of L is also lower-triangular, with 1's on the diagonal.  
If all elements of D are positive, then A must be symmetric
positive definite (SPD), and the solution can be reduced the usual
Cholesky decomposition U'U where U is upper triangular and
U = sqrt(D) L'.

The main advantage of the generalized form is that it admits
of matrices that are not of full rank: D will contain zeros 
marking the redundant columns, and the rank of A is the
number of non-zero columns.  If all elements of D are zero or
positive, then A is a non-negative definite (NND) matrix.
The generalized form also has the (quite minor) numerical advantage
of not requiring square roots during its calculation.
To extract the components of the decompostion, use the 
\code{diag} and \code{as.matrix}
functions.

The \code{solve} has a method for gchol decompostions,
and there are gchol methods for block diagonal symmetric
(\code{bdsmatrix}) matrices as well.
}
\keyword{array}
% docclass is function
% Converted by Sd2Rd version 1.21.
