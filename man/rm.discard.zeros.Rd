\name{rm.discard.zeros}
\alias{rm.discard.zeros}

\title{ Discard rows and columns from a matrix that exclusively contain zero-valued off-diagonal matrix elements}

\description{ The function removes those rows and columns from an input matrix that exclusively contain zero-valued off-diagonal elements.}


\usage{ rm.discard.zeros(mat, tol = 0, silent = FALSE)}
%- maybe also 'usage' for other objects documented here.


\arguments{
  \item{mat}{ Input matrix (typically after using \code{\link{rm.denoise.mat}})}
  \item{tol}{ A (small) real number specifying a thresholf for removal of matrix elements (see 'Details').}
  \item{silent}{ A logical variable that determines if the number of removed rows and columns is printed by the function or not.}
}


\details{
The diagonal of the matrix is not included when counting the zeros in a row/column, i.e. a row/column is actually removed if the diagonal element is the only non-zero element in that row/column. The tolerance \code{tol} specifies a threshold. Matrix elements below this threshold will be treated as if they were zero.
}


\value{ A matrix with zero-valued rows/colums removed.}

\author{ Uwe Menzel <uwemenzel@gmail.com>}

\seealso{ Remove noise from a random matrix by applying a threshold: \code{\link{rm.denoise.mat}} }

\examples{

set.seed(777)
random.matrix <- create.rand.mat(size = 1000, distrib = "norm")$rand.matr
dim(random.matrix)

## After identification of a proper threshold:
cleaned.matrix <- rm.denoise.mat(random.matrix, threshold = 3.2, keep.diag = TRUE)	
cl2.matrix = rm.discard.zeros(cleaned.matrix)				
 
}



