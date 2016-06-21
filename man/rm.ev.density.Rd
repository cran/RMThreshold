\name{rm.ev.density}
\alias{rm.ev.density}
\title{Create a density plot and a histogram of the eigenvalue distribution}

\description{
The function creates a density plot of the empirical distribution of the eigenvalues, combined with a histogram. Optionally, a curve illustrating the Wigner semi-circle can be added. The plot can be saved or shown in a plot window. Marks on the x-axis can be added optionally.
}
		
\usage{ rm.ev.density(eigenvalues, nr.breaks = 51, min.bw = 0.01, wigner = TRUE, 
	mark.on.x = NULL, title = "Eigenvalue density distribution", 
	pop.up = TRUE, fn = NULL) 
}

		
\arguments{
  \item{eigenvalues}{ A numeric vector containing the eigenvalues.}
  \item{nr.breaks}{ Number of bins used in the histogram.} 
  \item{min.bw}{ Minimum bandwidth for the calculation of the density curve. If the automatically calculated bandwidth gets too low, it is replaced by this value. That prevents the density curve from being too cliffy.}
  \item{wigner}{ A logical variable that determines if the Wigner semi-circle is to be added to the plot.} 
  \item{mark.on.x}{ A numeric vector or NULL. If not NULL, marks will be added on the x-axis at the positions given by the vector.}   
  \item{title}{ String containing the title of the plot.}  
  \item{pop.up}{ A logical variable that determines if the plot is to be shown in a plot window.}    
  \item{fn}{ A string determining the filename for storage. Must have extension 'png' or 'pdf'.}  
}
  
\value{ The name of the plot filename chosen, or NULL.}

\author{ Uwe Menzel \email{uwemenzel@gmail.com} }

\note{This function plots the density of the eigenvalues. For illustration of their spacings, use \code{rm.spacing.distribution}.}

\seealso{ Plotting the eigenvalue spacing distribution: \code{\link{rm.spacing.distribution}}}

\examples{
## Plot histogram of the spacings of the unfolded eigenvalues of a random matrix:
load(system.file("data", "random.matrix.RData", package="RMThreshold"))
res <- rm.ev.unfold(random.matrix)					  
rm.ev.density(res$eigenvalues, wigner = TRUE) 
}


