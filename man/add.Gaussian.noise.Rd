\name{add.Gaussian.noise}
\alias{add.Gaussian.noise}
\title{ Add Gaussian noise to a matrix}

\description{ The function adds Gaussian (i.e. normally distributed) noise to a matrix.}

\usage{add.Gaussian.noise(mat, mean = 0, stddev = 1, symm = TRUE)}

\arguments{
  \item{mat}{ Input matrix.}
  \item{mean}{ Mean of the Gaussian noise to be added.}
  \item{stddev}{ Standard deviation of the Gaussian noise to be added.}
  \item{symm}{ A logical variable that determines if the matrix is to be symmetrized after adding the noise.}
}


\details{
The function uses the \code{rnorm} function to create the normally distributed noise and adds it to the input matrix. Optionally, the matrix is symmetrized by adding it's transpose and dividing by \eqn{\sqrt 2}.}


\value{ The input matrix with noise added, optionally symmetrized.}

\references{
\url{https://en.wikipedia.org/wiki/Gaussian_noise}
}

\author{ Uwe Menzel <uwemenzel@gmail.com>}

\note{ The matrix can not be symmetrized if it is not quadratic.}


\seealso{ Random generation for the normal distribution: \code{\link[stats]{rnorm}} }

\examples{

\dontrun{
N = 500
some.mat = matrix(rep(1, N*N), nrow = N)	
some.mat[1:3, 1:10]
res <- rm.matrix.validation(some.mat)		# not really a proper matrix for this approach.
}

## It can help to add Gaussian noise to an improper matrix
\dontrun{
noisy.matrix <- add.Gaussian.noise(some.mat, mean = 0, stddev = 1, symm = TRUE)
noisy.matrix[1:3, 1:10]
res <- rm.matrix.validation(noisy.matrix)	# better!
res <- rm.get.threshold(noisy.matrix)		# 4.2	
}

}

\keyword{ models }
\keyword{ Gaussian noise }
\keyword{ Random Matrix Theory }
