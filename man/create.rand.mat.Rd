\name{create.rand.mat}
\alias{create.rand.mat}
\title{ Create a real-valued, symmetric random matrix}

\description{ The function creates a real-valued, symmetric random matrix of desired dimension. Two alternatives for the probability distribution of the matrix elements are provided.} 

\usage{
create.rand.mat(size = 1000, distrib = c("norm", "unif"), mean = 0, stddev = 1)
}

\arguments{
  \item{size}{ Dimension of the (quadratic) matrix.}
  \item{distrib}{ Desired probability distribution of the matrix elements. Can be \code{norm} or \code{unif}.}
  \item{mean}{ Desired mean of the normal distribution. Only active if \code{distrib = 'norm'} was chosen.}
  \item{stddev}{ Desired standard deviation of the normal distribution. Only active if \code{distrib = 'norm'} was chosen.}
}

\details{
The function creates a real-valued, symmetrical random matrix of desired dimension. Two alternatives for the probability distribution of the matrix elements are provided: normal and uniform. If \code{distrib = 'norm'}, the mean and the standard deviation can additionally be chosen. If \code{distrib = 'unif'}, the matrix elements will be uniformly distributed in the interval (-1,1). 
}


\value{
A list containing the following components:
  \item{mean.diag}{ The mean of the diagonal elements of the matrix.}
  \item{stddev.diag}{ The standard deviation of the diagonal elements of the matrix.}
  \item{mean.triangle}{ The mean of the upper triangle of the matrix (diagonal excluded).}
  \item{stddev.triangle}{ The standard deviation of the upper triangle of the matrix (diagonal excluded).}
  \item{rand.matrix}{ The random matrix created.}
}

\author{Uwe Menzel <uwemenzel@gmail.com>}

\examples{
  set.seed(777)
  random.matrix <- create.rand.mat(size = 1000, distrib = "norm")$rand.matrix
  dim(random.matrix)
}

\keyword{ models }
\keyword{ Random matrix }



