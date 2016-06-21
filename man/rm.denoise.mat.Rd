\name{rm.denoise.mat}
\alias{rm.denoise.mat}
\title{ Remove noise from a random matrix by applying a threshold}

\description{ Matrix elements with an absolute value below the given threshold are set to zero. }

\usage{ rm.denoise.mat(mat, threshold, keep.diag = TRUE)}

\arguments{
  \item{mat}{ The noisy input matrix.}
  \item{threshold}{ Numerical value of the threshold.}
  \item{keep.diag}{ A logical variable that determines if the diagonal of the matrix is thresholded or not. The default is \code{keep.diag = T}. In that case, diagonal matrix elements are not touched.}
}

\details{ The function outputs the number of non-zero matrix elements before and after thresholding.}

\value{ The thresholded matrix.}

\author{ Uwe Menzel <uwemenzel@gmail.com>}


\seealso{ Estimate an objective threshold for signal-noise separation in random matrices: \code{\link{rm.matrix.validation}} }


\examples{

load(system.file("data", "random.matrix.RData", package="RMThreshold"))
dim(random.matrix)


## After identification of a proper candidate threshold:
cleaned.matrix <- rm.denoise.mat(random.matrix, threshold = 3.2, keep.diag = TRUE)	  
  
}

