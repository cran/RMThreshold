\name{rm.matrix.validation}
\alias{rm.matrix.validation}
\title{Validate input matrix prior to threshold computation}

\description{
The function checks if the input matrix is well-conditioned for the algorithm used by \code{RMThreshold}. The matrix must be real-valued, symmetric, and large. Rank and sparseness of the matrix are calculated. Diagnostic plots are created. 
}

\usage{ rm.matrix.validation(rand.mat, unfold.method = "gaussian", 
		bandwidth = "nrd0", nr.fit.points = 51, discard.outliers = TRUE)
}

\arguments{
  \item{rand.mat}{ A random, real-valued, symmetric input matrix.}
  \item{unfold.method}{ A string variable that determines which type of unfolding algorithm is used. Must be one of 'gaussian' (Gaussian kernel density) or 'spline' (cubic spline interpolation on the cumulative distribution function).}
  \item{bandwidth}{ Bandwidth used to calculate the Gaussian kernel density. Only active if \code{unfold.method = 'gaussian'} is used. See the description of the \code{density} function.}   
  \item{nr.fit.points}{ Number of evenly spaced supporting points used for the cubic spline to the empirical cumulative distribution function.}  
  \item{discard.outliers}{ A logical variable that determines if outliers are to be discarded from the spectrum of eigenvalues.}     
}


\details{
The input matrix must be real-valued and symmetric (a correlation or mutual information matrix self-evidently is). The matrix must not be too sparse (if so, you are probably done without thresholding). The rank of the matrix must not be too low in order to obtain a sufficient number of non-zero eigenvalues. Furthermore, the matrix must not be too small because Random Matrix Theory applies for large (theoretically infinite) matrices only. The function creates a diagnostic plot, showing the empirical eigenvalue distribution and the distribution of the spacings between them. The eigenvalue distribution of the input matrix should approximately resemble the Wigner semi-circle, while the spacings should resemble the Wigner-Dyson distribution (Wigner surmise).   
}

\value{
A list containing the following entries:
  \item{sparseness}{ The sparseness of the input matrix.}
  \item{rank}{ The rank of the input matrix.}  
  \item{validation.plot}{ The name of the valdation plot.} 
  \item{unfold.plot}{ The name of the plot which can be used to check if eigenvalue unfolding worked correctly.}
  \item{nr.outliers.removed}{ The number of eigenvalue outliers that have been removed. Only if \code{discard.outliers = TRUE} was used.}  
}


\references{
  \url{https://en.wikipedia.org/wiki/Random_matrix} 
  
  Wigner, E. P. , \emph{Characteristic vectors of bordered matrices with infinite dimensions}, Ann. Math. 62, 548-564, 1955. 
  
  Mehta, M., \emph{Random Matrices}, 3nd edition. Academic Press, 2004.
    
  Furht, B. and Escalante, A. (eds.), \emph{Handbook of Data Intensive Computing}, Springer Science and Business Media, 2011.  
}

\author{ Uwe Menzel <uwemenzel@gmail.com> }

\note{It is highly recommended to check the input matrix using this function.}


\seealso{ Creating a random matrix: \code{\link{create.rand.mat}}}


\examples{

## Run with pre-defined random matrix:
load(system.file("data", "random.matrix.RData", package="RMThreshold"))
dim(random.matrix)		# 1000 1000   should be big enough

\dontrun{
res <- rm.matrix.validation(random.matrix)
res <- rm.matrix.validation(random.matrix, discard.outliers = FALSE)	
res <- rm.matrix.validation(random.matrix, unfold.method = "spline")
res <- rm.matrix.validation(random.matrix, unfold.method = "spline", discard.outliers = FALSE)
}

## Run with self-created random matrix:
\dontrun{
rmat <- create.rand.mat(size = 1000, distrib = "norm", mean = 0, stddev = 1)$rand.matrix	
res <- rm.matrix.validation(rmat)
res <- rm.matrix.validation(rmat, discard.outliers = FALSE)
res <- rm.matrix.validation(rmat, unfold.method = "spline")
res <- rm.matrix.validation(rmat, unfold.method = "spline", discard.outliers = FALSE)
}

\dontrun{
  library(igraph)

  ## Create noisy matrix and validate:
  g <- erdos.renyi.game(1000, 0.1)	
  adj = as.matrix(get.adjacency(g))
  rm.matrix.validation(adj)	# Wigner-Dyson case, unstructured matrix, noise

  ## Create modular (block-diagonal) matrix and validate:
  matlist = list()
  for (i in 1:4) matlist[[i]] = get.adjacency(erdos.renyi.game(250, 0.1))	
  mat <- bdiag(matlist)	# block-diagonal matrix 		 
  rm.matrix.validation(as.matrix(mat))	# Exponential case, modular matrix

}

}

% Add one or more standard keywords, see file 'KEYWORDS' in the R documentation directory.
\keyword{ models }
\keyword{ Random Matrix Theory }
\keyword{ Correlation matrix }
\keyword{ Mutual information }
\keyword{ Signal noise separation }



