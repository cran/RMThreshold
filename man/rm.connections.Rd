\name{rm.connections}
\alias{rm.connections}
\title{ Create ordered list of largest matrix elements}

\description{ The function creates a data frame which is sorted according to the (absolute) magnitude of the matrix elements.}

\usage{ rm.connections(mat, nr.list = 30, abs.val = TRUE, fn = NULL)}


\arguments{
  \item{mat}{ Input matrix.}
  \item{nr.list}{ Number of matrix elements to show.}
  \item{abs.val}{ Logical variable determining if absolute values should be used for sorting.}
  \item{fn}{ A file name. If not \code{NULL}, the data frame is saved to that file.}
}

\details{
This function can for instance be useful if pairs of samples with the largest correlation/mutual information are to be identified. By default, the matrix elements are sorted according to their absolute values. The list will not be saved if no filename is invoked, otherwise it will be saved to a tab-separated text file. 
}

\value{
A data frame containing the values of the largest \code{nr.list} matrix elements, together with the respective row- and column numbers. If present, the referring row- and column names are also included.
}

\author{ Uwe Menzel <uwemenzel@gmail.com>}

\examples{

load(system.file("data", "random.matrix.RData", package="RMThreshold"))
dim(random.matrix)

## After identification of a proper threshold:
cleaned.matrix <- rm.denoise.mat(random.matrix, threshold = 3.2, keep.diag = TRUE)	
cl2.matrix = rm.discard.zeros(cleaned.matrix)				 
  
}

\keyword{ models }

