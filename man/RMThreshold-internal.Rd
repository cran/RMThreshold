\name{RMThreshold-internal}
\alias{kb.distance}
\alias{rm.exp.distrib}
\alias{kld}
\alias{rm.get.file.extension}
\alias{rm.ev.unfold}
\alias{rm.spacing.scatter}
\alias{rm.trapez.int}
\alias{rm.reorder.ev}
\alias{rm.get.sparseness}
\alias{wigner.surmise}
\alias{wigner.semi.circle}
\alias{rm.get.distance}
\alias{rm.likelihood.plot}
\alias{rm.distance.plot}
\alias{rm.unfold.gauss}
\alias{rm.unfold.spline}
\alias{rm.sse.plot}
\alias{rm.sse}
\alias{rm.show.test}


\title{ Internal functions for the RMThreshold package }

\description{ Internal functions for the RMThreshold package}


\usage{
kb.distance(histo) 
rm.exp.distrib(x)
kld(observed, expected, plot) 
rm.get.file.extension(plotnames)
rm.ev.unfold(rand.mat, unfold.method, bandwidth, nr.fit.points, 
			discard.outliers, fn, pop.up, silent)  
rm.spacing.scatter(ev.spacing, title, pop.up, fn) 
rm.trapez.int(x, y)
rm.reorder.ev(eigenvalues, eigenvec) 
rm.get.sparseness(mat) 
wigner.surmise(x)
wigner.semi.circle(x)
rm.get.distance(ev.spacing, dist.method, nr.breaks)
rm.likelihood.plot(thresholds, log.le, log.lw, smooth.par, fn, interactive)
rm.distance.plot(thresholds, dist.Expon, dist.Wigner, smooth.par, fn, interactive)
rm.unfold.gauss(eigenvalues, bandwidth, fn, pop.up, silent)
rm.unfold.spline(eigenvalues, nr.fit.points, fn, pop.up)
rm.sse(ev.spacing, bandwidth, nr.points, N)
rm.show.test(thresholds, p.values, main, fn, interactive)
rm.sse.plot(thresholds, sse.values, main, fn, interactive)
}


\arguments{	
  \item{nr.fit.points}{ Number of supporting points used for the cubic spline to the empirical cumulative distribution function.}
  \item{discard.outliers}{ A logical variable that determines if outliers are discarded from the spectrum of eigenvalues.}  
  \item{silent}{ A logical variable that decides if a function outputs runtime messages or not.}  
  \item{histo}{ An R object of class 'histogram'. Output of the \code{hist} function. }
  \item{x}{ A real-valued number.}  
  \item{eigenvalues}{ A numerical vector containing the eigenvalues of a matrix.}
  \item{wigner}{ A logical variable that determines if the Wigner semi-circle or the Wigner surmise is added to a plot.}  
  \item{title}{ String variable containing the title of a plot.}
  \item{pop.up}{ A logical variable that determines if a plot window is supposed tp pop up during function execution.}  
  \item{fn}{ A filename.}  
  \item{observed}{ A numerical vector with the observed frequency of eigenvalues in each bin of the distribution function. }
  \item{expected}{ A numerical vector with the expected frequency of eigenvalues in each bin. The expected values refer to a limiting distribution, either Wigner-Dyson or Exponential distribution.}  
   \item{plot}{ A logical variable that determines if a plot should be created. }  
  \item{plotnames}{ A string or character vector of filenames for plots to be viewed. }  
  \item{rand.mat}{ A symmetric, real-valued random matrix. }  
  \item{fit}{ Name of the plot showing the cubic spline fit to the cumulative distribution. No plot will be made if \code{fit = NULL}} 
  \item{ev.spacing}{ A vector with the normalized eigenvalue spacings of a random matrix. }  
  \item{y}{ A numerical vector defining the y-values of a function. }  
  \item{eigenvec}{ A matrix containing the eigenvectors of a matrix (in columns).}  
  \item{mat}{ A real-valued matrix.} 
  \item{nr.breaks}{ Number of bins used in the histogram to subdivide the empirical eigenvalue spacing distribution.}   
  \item{log.le}{ Log likelihood of the observed eigenvalue distances when an exponential distribution is assumed.}   
  \item{log.lw}{ Log likelihood of the observed eigenvalue distances when a Wigner-Dyson distribution is assumed.}   
  \item{smooth.par}{ Parameter controlling the degree of smoothing in the loess regression curve presented in the final plot (distance vs. threshold).}     
  \item{interactive}{ A logical variable that determines if thresholds can be chosen by mouse clicks in the final plot (distance vs. threshold).} 
  \item{bandwidth}{ Bandwidth used to calculate the Gaussian kernel density. See description of the \code{density} function.} 
  \item{unfold.method}{ A string that decides which method is used for eigenvalue unfolding. One of 'gaussian' or 'spline'.}
  \item{dist.method}{ A string that determines which method is used to estimate the distance to the limiting distributions. One of 'LL' (Log Likelihood) or 'KLD' (Kullback-Leibler Distance).} 
  \item{thresholds}{ A numerical vector containing the values of the thresholds probed in the main function (\code{rm.get.threshold}).}
  \item{dist.Expon}{ A numerical vector containing the estimated distances to the Exponential distribution calculated in the main function (\code{rm.get.threshold}).}  
  \item{dist.Wigner}{ A numerical vector containing the estimated distances to the Exponential distribution calculated in the main function (\code{rm.get.threshold}).} 
  \item{p.values}{ A numerical vector containing the p-values for a significance test.}
  \item{sse.values}{ A numerical vector containing the Sum of Squared Errors between observed NNSD and Exponential function.}  
  \item{nr.points}{ Number of supporting points used to approximate the density function.}
  \item{N}{ Number of sections used to calculate the Sum of Squared Errors (SSE).}     
  \item{main}{ String variable containing the title of a plot.}
}


\details{ These functions are not intended to be called by the user.}




