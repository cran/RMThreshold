\name{RMThreshold-package}
\alias{RMThreshold-package}
\alias{RMThreshold}
\docType{package}

\title{ Signal-Noise Separation in Correlation Matrices by using Eigenvalue Spectrum Analysis}

\description{
The package provides an algorithm that can be used to determine an objective threshold for signal-noise separation in large random matrices.    
}



\details{
The package provides an algorithm which can be used to determine an objective threshold for signal-noise separation in large random matrices (correlation matrices, mutual information matrices, network adjacency matrices). The package makes use of the results of Random Matrix Theory (RMT). The algorithm increments a suppositional threshold monotonically, thereby recording the eigenvalue spacing distribution of the matrix. According to RMT, that distribution undergoes a characteristic change when the threshold properly separates signal from noise. The modular structure of the matrix (or of the corresponding network) can be unraveled if such a threshold is found. 
}


\author{
Uwe Menzel <uwemenzel@gmail.com>
}


\references{

\url{https://en.wikipedia.org/wiki/Random_matrix} 

Wigner, E. P. , \emph{Characteristic vectors of bordered matrices with infinite dimensions}, Ann. Math. 62, 548-564, 1955.

Mehta, M., \emph{Random Matrices}, 3nd edition. Academic Press, 2004.

Furht, B. and Escalante, A. (eds.), \emph{Handbook of Data Intensive Computing}, Springer Science and Business Media, 2011.

Luo, F. et al., \emph{Constructing gene co-expression networks and predicting functions of unknown genes by random matrix theory.} BMC Bioinformatics, 2007.

}

\keyword{models}
\keyword{ Random Matrix Theory }

\examples{
\dontrun{
  set.seed(777)
  random.mat <- create.rand.mat(size = 1000, distrib = "norm")$rand.matrix
  res <- rm.matrix.validation(random.mat)	# ok
  res <- rm.get.threshold(random.mat)		# threshold about 3.19	 			
  rm.show.plots(res$comparison.plots) 							
  cleaned.matrix <- rm.denoise.mat(random.mat, threshold = 3.2)	
  cleaned.matrix <- rm.discard.zeros(cleaned.matrix)
}					
}
