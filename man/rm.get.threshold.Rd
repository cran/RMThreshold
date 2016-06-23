\name{rm.get.threshold}
\alias{rm.get.threshold}
\title{Estimate an objective threshold for signal-noise separation in random matrices}

\description{
This is the main function of the package. A suppositional threshold is incremented monotonically, thereby recording the eigenvalue spacing distribution (nearest neighbor spacing distribution, NNSD) of the thresholded matrix. According to Random Matrix Theory, that distribution undergoes a characteristic change when the threshold properly separates signal from noise. By subsequent removal of the below-threshold matrix elements, the modular structure of the matrix - or of a network characterized by this matrix - can be unraveled. The function works for real symmetric matrices.
}
		

\usage{ rm.get.threshold(rand.mat, nr.thresholds = 51, 
    unfold.method = "gaussian", bandwidth = "nrd0", nr.fit.points = 51,
    dist.method = "LL", nr.breaks = 51, discard.outliers = TRUE, 
    discard.zeros = FALSE, min.mat.dim = 40, max.ev.spacing = 3, interval = NULL,
    interactive = TRUE, smooth.par = 0.5, plot.comp = TRUE,
    save.fit = FALSE, plot.spacing = FALSE, wait.seconds = 0)
}

		
\arguments{
  \item{rand.mat}{ A random, real-valued, symmetric matrix.}
  \item{nr.thresholds}{ Number of equidistant thresholds being probed and recorded.}
  \item{unfold.method}{ A string variable that determines which type of unfolding algorithm is used. Must be one of 'gaussian' (Gaussian kernel density) or 'spline' (cubic spline interpolation on the cumulative distribution function).}
  \item{bandwidth}{ Bandwidth used to calculate the Gaussian kernel density. Only active if \code{unfold.method = 'gaussian'} is used. See description of the \code{density} function.} 
  \item{dist.method}{ A string variable that determines which type of distance estimation to the limiting distributions is used. Must be one of 'LL' (Log Likelihood) or 'KLD' (Kullback-Leibler Distance).}
  \item{nr.fit.points}{ Number of supporting points used for the cubic spline to the empirical cumulative distribution function. Only active if \code{unfold.method = 'spline'} is used.}  
  \item{nr.breaks}{ Number of bins used in the analysis to subdivide the empirical eigenvalue spacing distribution.}     
  \item{smooth.par}{ Parameter controlling the degree of smoothing in the loess regression curve presented in the final plot (distance vs. threshold).} 
  \item{discard.outliers}{ A logical variable that determines if outliers are to be discarded from the spectrum of eigenvalues (see 'Details').}   
  \item{discard.zeros}{ A logical variable specifying if rows and columns exclusively containing zeros outside the main diagonal are to be removed from the matrix after applying each threshold (see 'Details').}   
  \item{min.mat.dim}{ By thresholding a matrix, rows and columns exclusively containing zeros in the off-diagonal matrix elements likely emerge. The parameter \code{min.mat.dim} determines the minimum number of non-zero rows and columns of the probed matrix during the thresholding. The thresholding loop is stopped if the number of non-zero rows and columns is getting below \code{min.mat.dim}.}  
  \item{max.ev.spacing}{ A number determining the maximum eigenvalue spacing that is considered when plotting and analyzing the empirical eigenvalue spacing distribution (see 'Details'). }   
  \item{interval}{ Interval of thresholds that is searched through. A numeric vector with two components (minimum and maximum threshold). The default is \code{interval = NULL} which means that the search interval ranges from the minimum absolute value of all matrix elements to the maximum absolute value of all matrix elements.} 
   \item{interactive}{ A logical variable that determines if the user wants to choose the candidate thresholds interactively, by mouse clicks.}     
  \item{plot.comp}{ A logical variable that determines if the empirical distribution of the eigenvalue spacings is displayed in a plot window during function execution.}      
  \item{save.fit}{ A logical variable that determines if a plot of the spline fit to the empirical cumulative distribution function is saved for each threshold. Can be used to check if fitting works well.} 
  \item{plot.spacing}{ A logical variable that determines if a scatterplot showing the eigenvalue spacing is saved at each suppositional threshold. Can be used to check if unfolding of the eigenvalue spacing works correctly (see 'Details').}  
  \item{wait.seconds}{ A numerical variable that, if set to non-zero values, enables viewing the plots with the eigenvalue spacing distribution during function execution for a predetermined duration. Useful on fast computers. Setting the variable to a specific number will cause the function to show the actual eigenvalue spacing distribution at least that number of seconds.}     
}

\details{
The function \code{rm.get.threshold} is the main function of the package. It takes a random matrix as input and applies a sequence of suppositional thresholds on it by setting all matrix elements to zero whose absolute value is below the actual threshold. The eigenvalues of the modified matrix are calculated at each threshold by using the \code{\link[base]{eigen}} function of the R \code{base} package. The eigenvalue spectrum is then unfolded, i.e. it is calibrated in such a way that the average global spacing between the eigenvalues is constant over the whole spectrum. The latter can be tracked by setting \code{plot.spacing = T}. Two methods are provided for unfolding: one method is based on calculation of the Gaussian kernel density of the eigenvalue spectrum; another method is based on fitting a cubic spline function to the cumulative empirical eigenvalue distribution. The latter is determined by the parameter \code{unfold.method}. See the references for details of the unfolding procedure. For each threshold, a distance between the empirical eigenvalue spacing distribution (NNSD) and both limiting distributions (Wigner-Dyson distribution and Exponential distribution, respectively) is estimated. Two methods of distance computation are implemented: a method based on computation of the log likelihood of the empirical eigenvalue spacing presupposing each of the limiting distributions, and a method based on calculation of the Kullback-Leibler divergence between empirical eigenvalue spacing distribution and these limiting distributions. If the assumed modular structure of the matrix is completely covered by noise, the empirical eigenvalue spacing distribution is close to the Wigner-Dyson distribution, a curve that approaches zero for small eigenvalue spacings. In the opposite case, when the modular structure of the matrix is prevailing, the empirical eigenvalue spacing distribution comes closer to an Exponential distribution (which represents the distribution of the intervals between two consecutive events in a Poisson process). If the matrix possesses a modular structure (hidden by noise), we expect that the NNSD changes gradually from the Wigner-Dyson case to the Exponential case when the threshold is increased stepwise. This change is monitored in a plot window if the \code{plot.comp} variable is left at it's default (\code{plot.comp = TRUE}). Two additional parameters are critical for proper functioning of the algorithm. For some types of input matrices, it might be necessary to remove the outliers of the eigenvalue distribution, in order to correctly investigate the bulk of the eigenvalue spectrum. This is achieved by the setting \code{discard.outliers = TRUE}, which is the default setting. In some other cases, it might be useful to retain the outliers during analysis. Another critical parameter is \code{discard.zeros}. If set to TRUE, rows and columns exclusively containing zeros outside the main diagonal are removed from the matrix at each threshold. This causes the matrix to shrink during the loop over thresholds. Setting \code{discard.zeros = TRUE} can be especially useful when the NNSD piles up at the left tail of the histogram shown during program execution. For very fast computers, the argument \code{wait.seconds} can be set to a non-zero value, in order to enable the user to follow that change of the NNSD visually. The distance between NNSD and the limiting distributions is not calculated over the whole range of eigenvalue spacings but over the interval \code{(0, max.ev.spacing)}. At a spacing of zero, the difference between the Wigner-Dyson distribution and the Exponential distribution is mostly pronounced. The maximum spacing considered in the distance calculation is determined by the parameter \code{max.ev.spacing}. This parameter should not be lower than \eqn{\sqrt(2 / \pi)}, where the peak of the Wigner-Dyson distribution lies. On the other hand, is does not make sense to choose too high values for \code{max.ev.spacing}, because both the Wigner-Dyson and the Exponential distribution assume rather low values in the right tail (which might cause numerical errors). If the algorithm works well, a relatively sharp transition from the Wigner-Dyson case to the Exponential case should become apparent. The latter is (hopefully) confirmed in a plot which is shown after completion of the loop over the sequence of suppositional thresholds. This plot shows the calculated distance between the NNSD and both limiting distributions versus applied threshold. The user can interactively choose candidate thresholds by clickung with the left mouse button on the points of the red-coloured curve. The selection is terminated by a right mouse-click somewhere on the plot. Likewise, candidate thresholds can be chosen in a plot showing the p-values for a Kolmogorov-Smirnov test (testing for exponentiality), and in a plot showing the Sum of Squared Errors between the empirical NNSD and the Exponential distribution versus threshold. The hereby chosen candidate thresholds are returned by the function. The analysis can now be refined by narrowing down the search interval for the thresholds by specifying the \code{interval} parameter.    
}


\value{
A list containing the following entries:
  \item{unfold.method}{ The method that was used for eigenvalue unfolding. Either 'gaussian' or 'spline'.}
  \item{dist.method}{ The method that was chosen to estimate the distance to the limiting distributions. Either 'LL' (Log Likelihood) or 'KLD' (Kullback-Leibler Divergence).}  
  \item{tested.thresholds}{ A vector containing the probed thresholds.}
  \item{dist.Wigner}{ A vector containing the numerical values of the estimated distance between empirical eigenvalue spacing distribution and Wigner-Dyson distribution.}  
  \item{dist.Expon}{ A vector containing the numerical values of the estimated distance between empirical eigenvalue spacing distribution and Exponential distribution.} 
  \item{nr.zeros}{ A vector containing the number of zero-valued matrix elements for each threshold.}    
  \item{nr.uniq.ev}{ An integer vector indicating the number of unique eigenvalues at each threshold (when degenerated eigenvalues are counted as one).}
  \item{max.ev.mult}{ An integer vector indicating the maximum eigenvalue multiplicity (degeneracy) at each threshold. } 
  \item{nr.spacings}{ An integer vector containing the number of eigenvalue spacings at each threshold. This number is smaller than the matrix dimension if eigenvalues are degenerated.}
  \item{nr.small.spacings}{ An integer vector containing the number of small spacings (\code{< max.ev.spacing/1000}) for each probed threshold.}
  \item{perc.small.spacings}{ A vector containing the percentage of small spacings for each probed threshold.}  
  \item{eff.dimension}{ An integer number specifying the 'effective dimension' of the matrix after thresholding, i.e. the number of non-zero rows and columns.}
  \item{comparison.plots}{ A character vector containing the names of the plots comparing the empirical eigenvalue spacing distribution with both limiting distributions (only if \code{plot.comp = T}).}   
  \item{rm.dimension}{ An integer vector indicating the dimension of the matrix after each thresholding (only if \code{discard.zeros = TRUE}).}
  \item{nr.outliers.removed}{ An integer vector containing the number of outliers in the eigenvalue spectrum removed at each threshold. Only if \code{discard.outliers = TRUE}.} 
  \item{p.ks}{ A vector containing the p-values for the Kolmogorov-Smirnov test at each probed threshold.}
  \item{sse.exp}{ A vector containing the Sum of Squared Errors (SSE) between observed NNSD and Exponential distribution for each probed threshold.}  
  \item{number.zeros.plot}{ A character string with the name of the plot depicting the number of zero-valued matrix elements versus threshold.}   
  \item{number.uniq.ev.plot}{ The name of the plot showing the number of unique eigenvalues vs. threshold.}  
  \item{max.ev.mult.plot}{ The name of the plot showing the maximum eigenvalue multiplicity vs. threshold.}
  \item{mat.dimension.plot}{ The name of the plot showing the matrix dimension after each thresholding step (only if \code{discard.zeros = TRUE}).}
  \item{num.ev.spacings.plot}{ The name of the plot showing the number of eigenvalue spacings vs. threshold.}   
  \item{distance.plot}{ The name of the (main) plot showing the distance of the empirical eigenvalue spacing distribution to the limiting distributions at each threshold.}   
  \item{cumfit.plots}{ A character vector containing the names of the plots showing the spline-function fitting of the cumulative eigenvalue spacing distribution (only if \code{save.fit = TRUE}).}  
  \item{space.plots}{ A character vector containing the names of the scatter plots with the eigenvalue spacing (only if \code{plot.spacing = TRUE}).}   
  \item{chosen.thresholds}{ A vector containing the potential thresholds chosen by the user from the distance plot.} 
  \item{p.ks.plot}{ The name of the plot showing the p-values for the Kolmogorov-Smirnow test versus probed thresholds.}   
  \item{p.ks.chosen}{ A vector containing the candidate thresholds chosen by the user based on the Kolmogorov-Smirnow test.}     
  \item{sse.plot}{ The name of the plot showing the SSE between observed NNSD and Exponential distribution versus probed thresholds.}  
  \item{sse.chosen}{ A vector containing the candidate thresholds chosen by the user based on the SSE test.}    
}

 


\references{
  \url{https://en.wikipedia.org/wiki/Random_matrix} 
  
  Wigner, E. P. , \emph{Characteristic vectors of bordered matrices with infinite dimensions}, Ann. Math. 62, 548-564, 1955.
  
  Mehta, M., \emph{Random Matrices}, 3nd edition. Academic Press, 2004.
  
  Furht, B. and Escalante, A. (eds.), \emph{Handbook of Data Intensive Computing}, Springer Science and Business Media, 2011.
  
  Luo, F. et al., \emph{Constructing gene co-expression networks and predicting functions of unknown genes by random matrix theory.} BMC Bioinformatics, 2007.
}

\author{ Uwe Menzel \email{uwemenzel@gmail.com} }

\note{It is recommended to check the input matrix using \code{rm.matrix.validation} before running this function. If the histogram of empirical eigenvalue spacings piles up at zero, one should set \code{discard.zeros = TRUE}}


\seealso{ Creating a random matrix: \code{\link{create.rand.mat}}}


\examples{

## Run with pre-defined random matrix:
set.seed(777)
random.matrix <- create.rand.mat(size = 1000, distrib = "norm")$rand.matr
dim(random.matrix)		# 1000 1000

\dontrun{
res <- rm.get.threshold(random.matrix)	# threshold might be 3.21
str(res)								# List of 26
rm.show.plots(res$comparison.plots)  	# watch sequence of plots once more
}

## Try other parameters:
\dontrun{
res <- rm.get.threshold(random.matrix, unfold.method = "spline")	
res <- rm.get.threshold(random.matrix, dist.method = "KLD")			
res <- rm.get.threshold(random.matrix, discard.outliers = FALSE) # might cause problems
res <- rm.get.threshold(random.matrix, wait.seconds = 2)	 	 # slow down
res <- rm.get.threshold(random.matrix, discard.zeros = TRUE)    	
res <- rm.get.threshold(random.matrix, discard.zeros = TRUE, dist.method = "KLD")	
}

## Refine analysis by choosing narrower threshold range 
\dontrun{
res <- rm.get.threshold(random.matrix, interval = c(2.5, 3.5))		
}

## Apply the identified threshold to the matrix
cleaned.matrix <- rm.denoise.mat(random.matrix, threshold = 3.21)	
cleaned.matrix <- rm.discard.zeros(cleaned.matrix)	 					  
dim(cleaned.matrix)	# smaller

## Find the clusters in the thresholded matrix:
\dontrun{
  library(igraph)
  g  <- graph.adjacency(cleaned.matrix, mode = "undirected")
  clusters(g)		
}


\dontrun{

  ## Create modular matrix and validate:
  matlist = list()
  for (i in 1:4) matlist[[i]] = get.adjacency(erdos.renyi.game(250, 0.1))	
  mat <- bdiag(matlist)					# create block-diagonal matrix 		 
  rm.matrix.validation(as.matrix(mat))	# Exponential case, modular matrix
  
  ## Add noise:
  mat1 = add.Gaussian.noise(as.matrix(mat), mean = 0, stddev = 0.1)
  
  ## Find threshold, denoise, reconstruct the modules:
  res <- rm.get.threshold(mat1)	# threshold possibly about 0.46
  # a smaller interval had been ok as well:
  res <- rm.get.threshold(mat1, interval = c(0, 0.8)) 
  cleaned <- rm.denoise.mat(mat1, 0.5)
  matr <- cleaned != 0
  g  <- graph.adjacency(matr, mode = "undirected")
  clusters(g)	# 4 clusters reconstructed
  
}					  

}



% Add one or more standard keywords, see file 'KEYWORDS' in the R documentation directory.
\keyword{ models }
\keyword{ Random Matrix Theory }
\keyword{ Correlation matrix }
\keyword{ Signal noise separation }



