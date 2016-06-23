\name{rm.spacing.distribution}
\alias{rm.spacing.distribution}
\title{Plot the empirical distribution of the eigenvalue spacings}

\description{
A histogram of the empirical distribution of the eigenvalue spacings is plotted. Optionally, curves illustrating the Wigner surmise and/or the Exponential distribution are added. 
}
			
\usage{ rm.spacing.distribution(ev.spacing, nr.breaks = 51, 
	wigner = TRUE, expo = TRUE, 
	title = "Eigenvalue spacing distribution (NNSD)", 
	threshold = NA, dist.Wigner = NA, 
	dist.Expon = NA, pop.up = TRUE, fn = NULL) 	
}	
	
		
\arguments{
  \item{ev.spacing}{ A numeric vector containing the spacings of the eigenvalues.}
  \item{nr.breaks}{ Number of bins used in the histogram.}  
  \item{wigner}{ A logical variable that determines if the Wigner-Dyson distribution (Wigner surmise) is to be added to a plot.}    
  \item{expo}{ A logical variable that determines if the Exponential distribution is to be added to the plot.}   
  \item{title}{ String containing the title of the plot.}  
  \item{threshold}{ If not NA, this value will be displayed in the plot, labeled 'threshold'.}    
  \item{dist.Wigner}{ If not NA, this value will be added to the plot, with a text indicating that it is the numerical value of the Kullback-Leibler distance between the empirical eigenvalue spacing distribution function and the Wigner-Dyson distribution function. }  
  \item{dist.Expon}{ If not NA, this value will be added to the plot, with a text indicating that it is the numerical value of the Kullback-Leibler distance between the empirical eigenvalue spacing distribution function and the Exponential distribution.}  
  \item{pop.up}{ A logical variable that determines if the plot is shown in a plot window.}    
  \item{fn}{ A string determining the filename for storage. Must have extension 'png' or 'pdf'.}  
}
  

\value{ The name of the plot filename chosen, or NULL.}

\author{ Uwe Menzel \email{uwemenzel@gmail.com} }

\seealso{ Plotting the eigenvalue distribution: \code{\link{rm.ev.density}}}

\examples{
## Plot histogram of the spacings of the unfolded eigenvalues of a random matrix:
set.seed(777)
random.matrix <- create.rand.mat(size = 1000, distrib = "norm")$rand.matr
res <- rm.ev.unfold(random.matrix)					  
rm.spacing.distribution(res$ev.spacing) 
}


