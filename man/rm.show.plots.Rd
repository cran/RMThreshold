\name{rm.show.plots}
\alias{rm.show.plots}

\title{ Display a sequence of plots on screen}

\description{ The function displays a sequence of plots saved beforehand (e.g. by \code{\link{rm.get.threshold}}).}

\usage{ rm.show.plots(plotnames)}


\arguments{
  \item{plotnames}{ A character string or -vector containing the names of the plots to be viewed.}
}


\details{ Can be useful when the user wants to inspect the sequence of empirical eigenvalue spacing distributions (repeatedly). }


\value{ No return values are being created.}


\author{ Uwe Menzel <uwemenzel@gmail.com>}

\note{ Might cause problems in some environments, only tested on Linux.}


\examples{

\dontrun{
set.seed(777)
random.matrix <- create.rand.mat(size = 1000, distrib = "norm")$rand.matr
dim(random.matrix)
res <- rm.get.threshold(random.matrix)					  
rm.show.plots(res$comparison.plots) # watch sequence of plots once more
} 
    
}



