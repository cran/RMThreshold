.packageName <- "RMThreshold"



## Check input random matrix for sparseness, rank, symmetry, ... and create diagnostic plots
"rm.matrix.validation" <- 
function(rand.mat, unfold.method = "gaussian", bandwidth = "nrd0", nr.fit.points = 51, discard.outliers = TRUE) 
{ 		
	if(!is.matrix(rand.mat))  stop("\n\n rm.matrix.validation: 'rand.mat' must be a matrix.\n\n") 
  	if(!is.numeric(rand.mat)) stop("\n\n rm.matrix.validation: 'rand.mat' must be a numeric matrix.\n\n") 
	if(!unfold.method %in% c("gaussian", "spline")) stop("\n\n  rm.matrix.validation: parameter 'unfold.method' must be one of 'gaussian' or 'spline'.\n\n")	 
 
  	N = nrow(rand.mat)
  	M = ncol(rand.mat)
  	cat(paste("\n ", N, "times", M, "matrix read.\n\n"))
	cat(paste("  Unfolding:", unfold.method, "\n"))
	if (unfold.method == "gaussian") cat(paste("  Bandwidth:", bandwidth, "\n"))
	if (unfold.method == "spline")   cat(paste("  Number of fit points:", nr.fit.points, "\n"))	
  
  	if(nrow(rand.mat) < 100) cat("\n   WARNING: 'rand.mat' seems to be quite small for this approach.\n\n")  
  	if(N != M) stop("\n rm.matrix.validation: 'rand.mat' must be quadratic.\n\n")
  	if(!isSymmetric(rand.mat)) stop("\n rm.matrix.validation: 'rand.mat' must be symmetric.\n\n") 
  	if(is.complex(rand.mat)) stop("\n rm.matrix.validation: 'rand.mat' must be real-valued.\n\n") 
 
  	nr.nonzero = sum(rand.mat != 0) 
  	perc.nonzero = nr.nonzero/N/N*100
  	sparseness = rm.get.sparseness(rand.mat)
  	rank.mat = Matrix::rankMatrix(rand.mat)[1]  
   
  	cat(paste("  Number of non-zero matrix elements:", nr.nonzero, "(", perc.nonzero, "% )\n")) 
  	cat(paste("  Sparseness:", sparseness, "\n")) 
  	cat(paste("  Rank:", rank.mat, "\n")) 
  
  	fn.uf = "fit.unfold.png"
	res <- rm.ev.unfold(rand.mat, unfold.method = unfold.method, bandwidth = bandwidth, nr.fit.points = nr.fit.points, 
						discard.outliers = discard.outliers, fn = fn.uf, pop.up = FALSE, silent = TRUE) 
	
	cat("\n")
  	eigenvalues = res$eigenvalues
  	ev.spacing = res$ev.spacing
   
  	## 0n-the-fly plot window with validation plots:
  	dev.new(units = "in", width = 10, height = 10, pointsize = 12, res = 96)
  	par(mfrow = c(2,2))
        
    ## Histogram of matrix elements:
    h1 = hist(rand.mat, breaks = 40, plot = F)
    plot(h1, main = "Magnitudes of matrix elements", font.main = 1, xlab = "magnitude", col = "red")

    ## Eigenvalue density (density plot)
    rm.ev.density(eigenvalues, wigner = F)
    
    ## Eigenvalue density (scatterplot with linear fit)
    rm.spacing.scatter(ev.spacing)
         
    ## Eigenvalue spacing density (NNSD):
    rm.spacing.distribution(ev.spacing)
    
  	par(mfrow = c(1,1)) 
	
	fn = "matrix.validation.plot.png"
	dev.copy(png, file = fn, units = "in", width = 10, height = 10, pointsize = 12, res = 96); dev.off()
		
	cat(paste("  Validation plot saved to:", fn, "\n"))
	cat(paste("  Unfolding plot saved to:", fn.uf, "\n\n"))	
	      
  	## Results:
  	results = list()
	results[["sparseness"]] = sparseness
	results[["rank"]] = rank.mat	
	results[["validation.plot"]] = fn
	results[["unfold.plot"]] = fn.uf
	results[["nr.outliers.removed"]] = res$nr.outliers.removed		
  	return(results)
} 








## Main function: loop through thresholds and record distance of NNSD to both limiting distributions (Wigner-Dyson & Exponential)
"rm.get.threshold" <- 
function(rand.mat, nr.thresholds = 51, unfold.method = "gaussian", bandwidth = "nrd0", nr.fit.points = 51, 
	dist.method = "LL", nr.breaks = 51, discard.outliers = TRUE, discard.zeros = FALSE, 
	min.mat.dim = 40, max.ev.spacing = 3, interval = NULL, interactive = TRUE, smooth.par = 0.5,
	plot.comp = TRUE, save.fit = FALSE, plot.spacing = FALSE, wait.seconds = 0) 
{
  	if(!is.matrix(rand.mat))   stop("\n\n rm.get.threshold: 'rand.mat' must be a matrix.\n\n") 
	if(nrow(rand.mat) != ncol(rand.mat)) stop("\n\n rm.get.threshold: 'rand.mat' must be a quadratic matrix.\n\n")
  	if(!isSymmetric(rand.mat)) stop("\n\n rm.get.threshold: 'rand.mat' must be a symmetric matrix.\n\n")
	if(max.ev.spacing <= sqrt(2/pi)) stop("\n\n rm.get.threshold: 'max.ev.spacing' should not be lower than the abscissa of the peak of the Wigner-Dyson distribution (which is approximately 0.8).\n\n")
    if(max.ev.spacing > 5)   cat("\n   WARNING: parameter 'max.ev.spacing' is quite big. Might provoke errors.\n\n")
  	if(nrow(rand.mat) < 100) cat("\n   WARNING: matrix 'rand.mat' seems to be quite small for this approach but let's give it a try ...\n\n")
	if(!unfold.method %in% c("gaussian", "spline")) stop("\n\n  rm.get.threshold: parameter 'unfold.method' must be 'gaussian' or 'spline'.\n\n")	 
	if(!dist.method %in% c("LL", "KLD")) stop("\n\n  rm.get.threshold: parameter 'dist' must be 'LL' (log likelihood) or 'KLD' (Kullback-Leibler).\n\n")	
	
	min.cell = min(abs(rand.mat[upper.tri(rand.mat, diag = F)]))
	max.cell = max(abs(rand.mat[upper.tri(rand.mat, diag = F)]))
	
	if(!is.null(interval)) {
		if(!is.numeric(interval) | (length(interval) != 2)) stop("\n\n rm.get.threshold: 'interval' must be a two-component numeric vector.\n\n")
		request  = ((min(interval) >= min.cell) || (min(interval) == 0)) && (max(interval) <= max.cell)
        if(!request) {
			cat(paste("\n  Range of the absolute values of the matrix elements:", signif(min.cell,5), " to ", signif(max.cell,5), "\n\n"))
			stop("\n  rm.get.threshold: parameter 'interval' must be inside the range of the absolute values of the matrix elements.\n\n")
		}
		thresholds = seq(min(interval), max(interval), len = nr.thresholds)
	} else {
    	thresholds = seq(min.cell, max.cell, len = nr.thresholds)	 
  	} 
 
  	N = nrow(rand.mat)
  	cat(paste("\n ", N, "times", N, "symmetric matrix read.\n\n")) 
  	        
  	results = list()
	results[["unfold.method"]] = unfold.method
	results[["dist.method"]] = dist.method
  	results[["tested.thresholds"]] = numeric(0) 
  	results[["dist.Wigner"]] = numeric(0)	
  	results[["dist.Expon"]] = numeric(0) 	
  	results[["nr.zeros"]] = integer(0)
	results[["nr.uniq.ev"]]	= integer(0)
	results[["max.ev.mult"]] = integer(0)
	results[["nr.spacings"]] = integer(0)
	results[["nr.small.spacings"]] = integer(0)	
	results[["perc.small.spacings"]] = integer(0)		
	results[["eff.dimension"]] = integer(0)		  
  	if(plot.comp) results[["comparison.plots"]] = character(0)  
  	if(save.fit) results[["cumfit.plots"]]      = character(0)  
  	if(plot.spacing) results[["space.plots"]]   = character(0) 
	if(discard.zeros) results[["rm.dimension"]] = integer(0)
	if(discard.outliers) results[["nr.outliers.removed"]] = integer(0)
	results[["p.ks"]] = integer(0)
	results[["sse.exp"]] = integer(0)	

			  
  	for (i in 1:nr.thresholds) {		    
    	loop.start.time <- Sys.time()  
    	thres = thresholds[i] 					
    	cat(paste(" ---> Loop =", i, "  threshold =", signif(thres,3), "\n"))

    	diagon = diag(rand.mat)					
    	rand.mat[which(abs(rand.mat) < abs(thres), arr.ind = T)] = 0	    
    	diag(rand.mat) = diagon					    

        eff.mat = rm.discard.zeros(rand.mat, silent = T)		
		cat(paste("      Effective matrix size =", nrow(eff.mat), "\n"))
								
		if(nrow(eff.mat) < min.mat.dim) {
			cat(paste("\n  Remaining number of non-zero rows & columns is below", min.mat.dim, "...\n  Breaking loop.\n\n")) 
			break
		}
        if(discard.zeros) rand.mat = eff.mat				 
				
    	if(save.fit) fn.fit = paste("RMT.Fit", i, "png", sep=".") else fn.fit <- NULL
		res <- rm.ev.unfold(rand.mat, unfold.method = unfold.method, bandwidth = bandwidth, nr.fit.points = nr.fit.points, discard.outliers = discard.outliers, fn = fn.fit, pop.up = FALSE, silent = TRUE) 		
    	ev.spacing = res$ev.spacing
		ev = res$eigenvalues
		

		l1 = length(ev.spacing)
		if(!is.null(max.ev.spacing)) ev.spacing = ev.spacing[ev.spacing <= max.ev.spacing]	
		l2 = length(ev.spacing)
		cat(paste("      Number of large spacings not considered ( larger than", max.ev.spacing, ") :", l1-l2, "\n"))				
	
		
		## Trace pileup at zero-NNSD:
		epsilon = max.ev.spacing/1000	## experimental 
		nr.small.spacings = sum(ev.spacing < epsilon)
		perc.small.spacings = nr.small.spacings/length(ev.spacing) * 100
		cat(paste("      Percentage of small spacings ( <", epsilon, ") =", round(perc.small.spacings,2), "\n"))		
		
		
		## Kolmogorov-Smirnov test:
		# p.val.ks.test = suppressWarnings(ks.test(ev.spacing, 'pexp', 1)$p.value)	# warnings might be generated because of ties
		p.val.ks.test = ks.test(unique(ev.spacing), 'pexp', 1)$p.value				# unique(ev.spacing) removes ties		

		
		## SSE between NNSD and exponential density:
		sse.exp = rm.sse(ev.spacing)	

												
		results[["eff.dimension"]][i] = nrow(eff.mat)
		results[["nr.small.spacings"]][i] = nr.small.spacings 
		results[["perc.small.spacings"]][i] = perc.small.spacings
		results[["nr.spacings"]][i] = length(ev.spacing) 		
		results[["cumfit.plots"]][i] = fn.fit					
	    results[["tested.thresholds"]][i] = thres						
    	results[["nr.zeros"]][i] = sum(rand.mat == 0)					
		results[["nr.uniq.ev"]][i] = length(unique(ev))					
		results[["max.ev.mult"]][i] = max(table(ev))					
		if(discard.zeros) results[["rm.dimension"]][i] = nrow(rand.mat)	
		if(discard.outliers) results[["nr.outliers.removed"]][i] = res$nr.outliers.removed 
		results[["p.ks"]][i] = p.val.ks.test
		results[["sse.exp"]][i] = sse.exp						
		
    	if(plot.spacing) {	
      		fn = paste("RMT.Spaceplot", i, "png", sep = ".")
      		rm.spacing.scatter(ev.spacing, pop.up = F, fn = fn) 
      		results[["space.plots"]][i] = fn
    	}  

 
    	dres <- rm.get.distance(ev.spacing, dist.method = dist.method, nr.breaks = nr.breaks)
		
		results[["dist.Wigner"]][i] = dres$dist.Wigner
		results[["dist.Expon"]][i]  = dres$dist.Expon
		

   		if(plot.comp) { 	
      		fn = paste("RMT.Spacing", i, "png", sep = ".")
      		rm.spacing.distribution(ev.spacing, threshold = thres, dist.Wigner = dres$dist.Wigner, dist.Expon = dres$dist.Expon, fn = fn)
      		results[["comparison.plots"]][i] = fn
    	}

    	loop.end.time <- Sys.time()
    	loop.time = as.numeric(loop.end.time - loop.start.time)	
    	if(loop.time < wait.seconds) Sys.sleep(wait.seconds - loop.time) 
  	}  	# End of loop
  
  
  	thresholds   = results[["tested.thresholds"]] 		# abbreviations
	nr.unique.ev = results[["nr.uniq.ev"]] 
	max.ev.mult  = results[["max.ev.mult"]]
  	nr.zeros     = results[["nr.zeros"]] 
	nr.spacings  = results[["nr.spacings"]]
  	dist.Wigner  = results[["dist.Wigner"]]
  	dist.Expon   = results[["dist.Expon"]]
	sum.sq.err   = results[["sse.exp"]]
	if(discard.zeros) mat.dim = results[["rm.dimension"]] 
  
  
  	## Plot number of zero matrix-elements vs. threshold 
  	fn = paste("RMT.numzeros", "png", sep = ".")
  	png(fn)	
	mtxt = "Number of zero matrix-elements vs. threshold"				
  	plot(thresholds, nr.zeros, col = "blue", main = mtxt, font.main = 1, xlab = "threshold", ylab = "nr zeros")
  	dev.off()
  	results[["number.zeros.plot"]] = fn 
	
	## Plot number of unique eigenvalues vs. threshold
  	fn = paste("RMT.num.uniq.ev", "png", sep = ".")	
  	png(fn)		
	mtxt = "Number of unique eigenvalues vs. threshold"	
  	plot(thresholds, nr.unique.ev, col = "blue", main = mtxt, font.main = 1, xlab = "threshold", ylab = "nr unique ev")	
  	dev.off()	
  	results[["number.uniq.ev.plot"]] = fn	
	
	## Plot maximum eigenvalue multiplicity vs. threshold
  	fn = paste("RMT.max.ev.mult", "png", sep = ".")	
  	png(fn)	
	mtxt = "Maximum eigenvalue multiplicity vs. threshold"	
  	plot(thresholds, max.ev.mult, col = "blue", main = mtxt, font.main = 1, xlab = "threshold", ylab = "max. ev multiplicity")		
  	dev.off()	
  	results[["max.ev.mult.plot"]] = fn	
	
	## Plot dimension of matrix vs. threshold
	if(discard.zeros) {
  		fn = paste("RMT.mat.dimension", "png", sep = ".")	
		png(fn)	
		mtxt = "Dimension of non-zero matrix vs. threshold"
		plot(thresholds, mat.dim, col = "blue", main = mtxt, font.main = 1, xlab = "threshold", ylab = "matrix dimension")
		dev.off()
		results[["mat.dimension.plot"]] = fn	
	}
	
	## Plot number of eigenvalue spacings vs. threshold    
	fn = paste("RMT.num.ev.spacings", "png", sep = ".")
	png(fn)	
	mtxt = "Number of ev spacings vs. threshold"
	plot(thresholds, nr.spacings , col = "blue", main = mtxt, font.main = 1, xlab = "threshold", ylab = "nr. ev spacings")	
	dev.off()
	results[["num.ev.spacings.plot"]] = fn
	
			
  	## Distance of emprical eigenvalues spacing distribution to both limiting distributions (Wigner surmise and Exp., resp.) versus chosen threshold:   	
	fn = paste("RMT.Dist.vs.Thres", "png", sep = ".") 	
	if(dist.method == "LL")  main.res <- rm.likelihood.plot(thresholds,  log.le = dist.Expon,      log.lw = dist.Wigner, smooth.par = smooth.par, fn = fn, interactive = interactive)
  	if(dist.method == "KLD") main.res <- rm.distance.plot(thresholds, dist.Expon = dist.Expon, dist.Wigner= dist.Wigner, smooth.par = smooth.par, fn = fn, interactive = interactive)
	cat(paste("\n  Distance plot saved to '", fn, "' \n\n"))
    results[["distance.plot"]] = fn		
	results[["chosen.thresholds"]]  = main.res$chosen.thresholds  
		
	
	## p-value for KS-Test versus thresholds:	
	fn = "RMT.pks.vs.Thres.png"
	maintxt = "p-values for KS-test"
	tres <- rm.show.test(thresholds = thresholds, p.values = results[["p.ks"]], main = maintxt, fn = fn)
	cat(paste("\n  Plot with KS-test saved to '", fn, "' \n\n"))
	results[["p.ks.plot"]]  = fn
	results[["p.ks.chosen"]] = tres$chosen


	## SSE between observed NNSD and Exponential distribution:
	fn = "RMT.SSE.vs.Thres.png"
	maintxt = "SSE for NNSD <--> Exponential"
	sse.res <- rm.sse.plot(thresholds = thresholds, sse.values = sum.sq.err, main = maintxt, fn = fn)	
	cat(paste("\n  SSE plot saved to '", fn, "' \n\n"))
	results[["sse.plot"]]  = fn
	results[["sse.chosen"]] = sse.res $chosen
		
	return(results)
} 





## Unfold eigenvalue spacing distribution using spline function fit 
"rm.ev.unfold" <- 
function(rand.mat, unfold.method = "gaussian", bandwidth = "nrd0", nr.fit.points = 51, discard.outliers = TRUE, fn = NULL, pop.up = FALSE, silent = TRUE) 
{
  	if(!is.matrix(rand.mat)) stop("\n\n  rm.ev.unfold: 'rand.mat' must be a matrix.\n\n")
  	if(!isSymmetric(rand.mat)) stop("\n\n  rm.ev.unfold: 'rand.mat' must be symmetric.\n\n") 
	if(!is.null(fn)) if(rm.get.file.extension(fn) != "png") stop("\n\n  rm.ev.unfold: 'fn' must be 'NULL' or some filename with extension '.png'.\n\n")	
	if(!unfold.method %in% c("gaussian", "spline")) stop("\n\n  rm.ev.unfold: 'unfold.method' must be one of 'gaussian' or 'spline'.\n\n")

	remove.outliers <- function(x, factor = 1.5) {
		q25 = quantile(x, probs = 0.25); q75 = quantile(x, probs = 0.75)		
		iqr = unname(q75 - q25)				
		lower.threshold = q25 - (iqr * factor); upper.threshold = q75 + (iqr * factor)	
		return(x[(x >= lower.threshold) & (x <= upper.threshold)])
	}

  	eigenvalues = eigen(rand.mat, only.values = T)$values			
  	eigenvalues = eigenvalues[order(eigenvalues)]/max(abs(eigenvalues))
   	
  	## Remove eigenvalue outliers:
	if(discard.outliers) {
		orig.number.ev = length(eigenvalues)
		eigenvalues = remove.outliers(unique(eigenvalues))
		new.number.ev = length(eigenvalues)
		nr.outliers.removed = orig.number.ev - new.number.ev
	} else {
		nr.outliers.removed = 0
	}
	if(!silent) cat(paste("  Number of discarded outlier eigenvalues:", nr.outliers.removed, "\n"))	
	
			
	if(unfold.method == "gaussian") {	## Use density function obtained by Gaussian broadening
		uf <- rm.unfold.gauss(eigenvalues, bandwidth = bandwidth, fn = fn, pop.up = pop.up, silent = TRUE)		
	}
			
	if(unfold.method == "spline") {		## Use cumulative distribution function together with cubic spline interpolation
	    nr.fit.points = min(nr.fit.points, floor(0.75*length(eigenvalues)))	## experimental: the fit must be smooth, must NOT be too close to the step function (cumulative density) 
		uf <- rm.unfold.spline(eigenvalues, nr.fit.points = nr.fit.points, fn = fn, pop.up = pop.up)	
	}
	
	## Results:
  	results = list() 
	results[["unfold.method"]] = unfold.method   
  	results[["eigenvalues"]] = uf$eigenvalues
	if(unfold.method == "gaussian") results[["bandwidth"]]  = uf$bandwidth
	if(unfold.method == "spline") results[["nr.fit.points"]] = uf$nr.fit.points
  	results[["unfolded.ev"]] = uf$unfolded.ev  	   
  	results[["ev.spacing"]]  = uf$ev.spacing
	results[["plot"]] = uf$plot
	results[["nr.outliers.removed"]] = nr.outliers.removed
  	return(results)   
} 





"rm.unfold.gauss" <-
function(eigenvalues, bandwidth = "nrd0", fn = NULL, pop.up = FALSE, silent = TRUE)
{
	if(!is.numeric(eigenvalues)) stop("\n\n  rm.unfold.density: 'eigenvalues' must be a numeric vector.\n\n")
	if(!is.null(fn)) if(rm.get.file.extension(fn) != "png") stop("\n\n  rm.unfold.density: filename must have extension '.png'. \n\n")

	dens = density(eigenvalues, bw = "nrd0", adjust = 1, n = 512, kernel = "gaussian", give.Rkern = FALSE)
	bandw = dens$bw
	if(!silent) cat(paste("  Bandwith for Gaussian broadening =", signif(bandw,5), "\n"))

	midpoints <- function(x) return(x[-length(x)] + 0.5*diff(x))
	scale.function = approx(dens$x, dens$y, xout = midpoints(eigenvalues))
		
	ev.spacing = diff(eigenvalues)
	ev.spacing = ev.spacing * scale.function$y
	ev.spacing = ev.spacing / mean(ev.spacing)
	
	uf.ev = cumsum(ev.spacing)
	uf.ev = uf.ev/max(uf.ev)				

	## Plots:	
	make.plot <- function() {
	
		## Calculated density and the scale points derived from it (should agree)
		plot(scale.function$x, scale.function$y, xlab = "eigenvalues", ylab = "density", main = "Eigenvalue density", font.main = 1, cex = 0.8, col = "blue", pch = 4)
		lines(dens, col = "red", lwd = 1.5) 
						
		## Empirical density distribution of eigenvalues (histogram) with density:
		br = seq(min(eigenvalues), max(eigenvalues), len = 51)
		yl = "PDF of eigenvalues"
		hist(eigenvalues, breaks = br, freq = F, col = "darkolivegreen2", main = "Eigenvalue density distribution", font.main = 1, xlab = "eigenvalue", ylab = yl)	
		lines(dens, col = "red", lwd = 1.5) 	
				
		## original and unfolded eigenvalues, ordered, plotted vs. index
		ylim = c(min(eigenvalues, uf.ev), max(eigenvalues, uf.ev)) 
		plot(eigenvalues, cex=0.4, main = "Ordered eigenvalues", font.main = 1, ylim = ylim)	# nonlinear
		points(uf.ev, cex=0.4, col = "red")														# linear
		legend("topleft", c("original", "unfolded"), col = c("black", "red"), lty = 1, lwd = 2)
		
		## Eigenvalue spacing distribution, Wigner-Dyson distribution, Exponential distribution
		rm.spacing.distribution(ev.spacing) 
	}
	
	if(pop.up) {
 		dev.new(units = "in", width = 10, height = 10, pointsize = 12, res = 96)
  		par(mfrow = c(2,2))
		make.plot()	
		par(mfrow = c(1,1))
	}
	
	if(!is.null(fn)) {	
		png(filename = fn, units = "in", width = 10, height = 10, pointsize = 12, res = 96)
		par(mfrow = c(2,2))
			make.plot()	
		par(mfrow = c(1,1))					
		dev.off()	
	}
	
	results = list()
	results[["eigenvalues"]] = eigenvalues
	results[["bandwidth"]] = bandw
	results[["unfolded.ev"]] = uf.ev
	results[["ev.spacing"]]  = ev.spacing 
	if(!is.null(fn)) results[["plot"]] = fn
	return(results)
}






"rm.unfold.spline" <-
function(eigenvalues, nr.fit.points = 51, fn = NULL, pop.up = FALSE)
{
	if(!is.numeric(eigenvalues)) stop("\n\n  rm.unfold.spline: 'eigenvalues' must be a numeric vector.\n\n")
	if(!is.null(fn)) if(rm.get.file.extension(fn) != "png") stop("\n\n  rm.unfold.spline: filename must have extension '.png'. \n\n")
	
	cumulative.distrib = ecdf(eigenvalues)		
	support = seq(min(eigenvalues), max(eigenvalues), len = nr.fit.points)									               
	D = splinefun(support, cumulative.distrib(support) , method = "hyman") 	 		

	uf.ev = D(eigenvalues)						# unfolded eigenvalues						
	ev.spacing = diff(uf.ev)					# spacing betweem them							
	ev.spacing = ev.spacing / mean(ev.spacing)	# normalized spacing		

	## Plots:	
	make.plot <- function() {
	
    	## Eigenvalue density (scatterplot with linear fit)
    	rm.spacing.scatter(ev.spacing)
 	
		## Cumulative distribution function with spline fit
		maintxt = "Cumulative distribution function with spline fit"	
		mtxt = 	paste("nr. eigenvalues =", length(eigenvalues), "  nr. fit points =", nr.fit.points)			
   		plot(eigenvalues, cumulative.distrib(eigenvalues), pch = 20, cex = 0.6, col = "red", main = maintxt, font.main = 1, xlab = "eigenvalue", ylab = "F(eigenvalue)")		
    	curve(D, from = min(eigenvalues), to = max(eigenvalues), col = "blue", add = T, cex = 1.2)
    	mtext(side = 3, mtxt, col = "blue")    
					
		## Original and unfolded eigenvalues    ordered, plotted over index
		ylim = c(min(eigenvalues, uf.ev), max(eigenvalues, uf.ev)) 
		plot(eigenvalues, cex=0.4, main = "Ordered eigenvalues", font.main = 1, ylim = ylim)	# nonlinear
		points(uf.ev, cex=0.4, col = "red")														# linear
		legend("topleft", c("original", "unfolded"), col = c("black", "red"), lty = 1, lwd = 2)
		
		## Eigenvalue spacing distribution, Wigner-Dyson distribution, Exponential distribution
		rm.spacing.distribution(ev.spacing) 
	}
	
	if(pop.up) {
 		dev.new(units = "in", width = 10, height = 10, pointsize = 12, res = 96)
  		par(mfrow = c(2,2))
		make.plot()	
		par(mfrow = c(1,1))
	}
	
	if(!is.null(fn)) {	
		png(filename = fn, units = "in", width = 10, height = 10, pointsize = 12, res = 96)
		par(mfrow = c(2,2))
			make.plot()	
		par(mfrow = c(1,1))					
		dev.off()	
	}
	
	results = list()
	results[["eigenvalues"]] = eigenvalues
	results[["nr.fit.points"]] = nr.fit.points
	results[["unfolded.ev"]] = uf.ev
	results[["ev.spacing"]]  = ev.spacing 
	if(!is.null(fn)) results[["plot"]] = fn
	return(results)
}






"rm.get.distance" <-
function(ev.spacing, dist.method = "LL", nr.breaks = 51)
{
	if(!dist.method %in% c("LL", "KLD")) stop("\n\n  rm.get.distance: parameter 'dist.method' must be 'LL' (log likelihood) or 'KLD' (Kullback-Leibler).\n\n") 

	results = list()
	results[["dist.method"]] = dist.method
	if(dist.method == "KLD") results[["nr.breaks"]] = nr.breaks
			
	if(dist.method == "LL") {						## Log Likelihood (per eigenvalue)   
		evs = ev.spacing[ev.spacing != 0]
		N = length(evs)
		log.LE = -sum(evs)/N
		log.LW = log(pi/2) + sum(log(evs))/N - 0.25*pi*sum(evs^2)/N 
		results[["dist.Wigner"]] = log.LW 			## log likelihood when Wigner distribution is assumed
		results[["dist.Expon"]]  = log.LE			## log likelihood when Exponential distribution is assumed
	}
	
	if(dist.method == "KLD") {						## Kullback-Leibler Divergence
		hsp <- hist(ev.spacing, breaks = seq(min(ev.spacing), max(ev.spacing), len = nr.breaks), plot = F)	 
		res.dist = kb.distance(hsp)
		results[["dist.Wigner"]] = res.dist$klw		## Kullback-Leibler divergence  observed -- Wigner surmise 
		results[["dist.Expon"]]  = res.dist$klp		## Kullback-Leibler divergence  observed -- Exponential distribution    
	}
	
	return(results)
}





"rm.sse" <-
function(ev.spacing, bandwidth = "nrd0", nr.points = 1000, N = 20)  # subdivide into N section with equal area
{
	dens = density(ev.spacing, bw = bandwidth, kernel = "gaussian", n = 512)		
	x = seq(min(ev.spacing), max(ev.spacing), len = nr.points)	
	observed = approx(dens$x, dens$y, xout = x)$y					# calculate density at the support points
	expected = exp(-x)	
	A = exp(-min(ev.spacing)) - exp(-max(ev.spacing))		
	
	xs <- numeric(N+1)
	xs[1] = min(ev.spacing)
	for (i in 1:N) xs[i+1] = -log(exp(-xs[i]) - A/N)

	area = numeric(N) 	## area under density curve (for observed) for each section
	for (i in 1:N) {	
		xsec = x[(x > xs[i]) & (x < xs[i+1])]
		xsec = c(xs[i], xsec, xs[i+1])								
		ysec = approx(dens$x, dens$y, xout = xsec)$y
		area[i] = rm.trapez.int(xsec, ysec)
	}
		
	SSE = sum((area[i] - A/N)^2)		
	return(SSE)
}







"rm.likelihood.plot" <-
function(thresholds, log.le, log.lw, smooth.par = 0.5, fn = NULL, interactive = TRUE)
{
	t = thresholds
  	le = -log.le	## display positive values
  	lw = -log.lw

	matplot(t, cbind(lw, le), col = c("red", "blue"), type = "p", pch = c(1,2), ylab = "ln(L)", xlab = "threshold", cex = 0.8)
  	title(main = "Log Likelihoods", font.main = 1)
  	legend("left", c("to Wigner", "to Exponential"), col = c("red", "blue"), pch = c(1,2), bty = "n") 

	lo.le <- loess(le ~ t, span = smooth.par)		
	lo.lw <- loess(lw ~ t, span = smooth.par)	

	xl <- seq(min(t), max(t), (max(t) - min(t))/1000)
	lines(xl, predict(lo.le, xl), col="blue", lwd=1, lty=3)
	lines(xl, predict(lo.lw, xl), col="red",  lwd=1, lty=3)
	mtext(side = 3, "(Negative log-likelihood per eigenvalue)", col = "darkgreen", cex = 0.9)
		
	if(interactive) {
    	cat("\n  Mark proper threshold by left-clicking on red point(s), stop with right-click\n\n")
    	coord <- suppressWarnings(sapply(list(t, lw),"[",identify(t, lw, labels = as.character(signif(t,4))))) 	
 		if(is.matrix(coord)) chosen = coord[,1] else chosen = coord[1]
		if(is.list(chosen)) chosen = "not chosen"
	}

    if(!is.null(fn)) {
		dev.copy(png, file = fn); dev.off()	
	}
    
	results = list()
	results[["loess.exp"]] = lo.le
	results[["loess.wig"]] = lo.lw
	if(interactive) results[["chosen.thresholds"]] = chosen else results[["chosen.thresholds"]] = "none"
	return(results) 
}





"rm.distance.plot" <-
function(thresholds, dist.Expon, dist.Wigner, smooth.par = 0.5, fn = NULL, interactive = TRUE)
{	
	t = thresholds
	le = dist.Expon
	lw = dist.Wigner
	
	matplot(t, cbind(lw, le), col = c("red", "blue"), type = "p", pch = c(1,2), ylab = "KLD", xlab = "threshold", cex = 0.8)
  	title(main = "Kullback-Leibler distance", font.main = 1)
  	legend("topleft", c("to Wigner", "to Exponential"), col = c("red", "blue"), pch = c(1,2), bty = "n") 

	lo.le <- loess(le ~ t, span = smooth.par)	
	lo.lw <- loess(lw ~ t, span = smooth.par)	

	xl <- seq(min(t), max(t), (max(t) - min(t))/1000)
	lines(xl, predict(lo.le, xl), col="blue", lwd=1, lty=3)
	lines(xl, predict(lo.lw, xl), col="red",  lwd=1, lty=3)
	mtext(side = 3, "(with loess regression curves)", col = "darkgreen", cex = 0.9)
		
	if(interactive) {
    	cat("\n  Mark proper threshold by left-clicking on red point(s), stop with right-click\n\n")
    	coord <- suppressWarnings(sapply(list(t, lw),"[",identify(t, lw, labels = as.character(signif(t,4))))) 	
 		if(is.matrix(coord)) chosen = coord[,1] else chosen = coord[1]
		if(is.list(chosen)) chosen = "not chosen"
	}

    if(!is.null(fn)) {
		dev.copy(png, file = fn); dev.off()	
	}
    
	results = list()
	results[["loess.exp"]] = lo.le
	results[["loess.wig"]] = lo.lw	
	if(interactive) results[["chosen.thresholds"]] = chosen else results[["chosen.thresholds"]] = "none"
	return(results) 
}






## Empirical eigenvalue spacing distribution and both limitig distributions
"rm.spacing.distribution" <- 
function(ev.spacing, nr.breaks = 51, wigner = TRUE, expo = TRUE, title = "Eigenvalue spacing distribution (NNSD)", threshold = NA, dist.Wigner = NA, dist.Expon = NA, pop.up = TRUE, fn = NULL) 
{   	
	if(!is.vector(ev.spacing)) stop("\n\n  rm.spacing.distribution: argument 'ev.spacing' must be a numeric vector.\n\n")
	if(class(ev.spacing) != "numeric") stop("\n\n  rm.spacing.distribution: argument 'ev.spacing' must be a numeric vector.\n\n")	
  	if((pop.up == F) && is.null(fn)) stop("\n\n  rm.spacing.distribution: Nothing to do. Enter filename or set 'pop.up = T'\n\n") 
  	if(!is.null(fn)) {			
    	ftype = rm.get.file.extension(fn) 	
    	if(!ftype %in% c("png", "pdf")) stop("\n\n  rm.spacing.distribution: Output file must have extension 'png' or 'pdf'.\n\n")
  	}  
	
	histo <- hist(ev.spacing, breaks = seq(min(ev.spacing), max(ev.spacing), len = nr.breaks), plot = F)  
  	make.plot <- function(histo = histo, title = title, threshold = threshold, dist.Wigner = dist.Wigner, dist.Expon = dist.Expon) {
    	plot(histo, freq = F, col = "darkolivegreen2", main = title, font.main = 1, xlab = "eigenvalue spacing", ylab = "PDF of eigenvalue spacing")			
    	
		if(expo) {
			actual.ymax = par("yaxp")[2]      	
    		x0 = -log(actual.ymax*0.98)	    	
    		curve(rm.exp.distrib,	from = max(x0,min(histo$breaks)), to = max(histo$breaks), n = 1001, add = T, type = "l", lty = 1, col = "red", lwd = 2)
    	}		
		if(wigner) curve(wigner.surmise, from = min(histo$breaks), to = max(histo$breaks), n = 1001, add = T, type = "l", lty = 1, col = "blue3", lwd = 2)
    	if((!is.na(dist.Wigner))&&(!is.na(dist.Expon))) mtext(side=3, paste("Distance Obs <=> Wigner =", signif(dist.Wigner,3), "    Distance Obs <=> Exponential =", signif(dist.Expon,3)), col = "brown", cex = 0.8)  
    	if(!is.na(threshold)) mtext(side = 4, paste("threshold =", signif(threshold,4)))
    	legend("topright", inset = 0.05, c("Exponential", "Wigner surm."), col = c("red", "blue"), lty = 1, lwd = 2, cex = 0.8)  
  	}   
  
  	if(pop.up) {				
    	make.plot(histo = histo, title = title, threshold = threshold, dist.Wigner = dist.Wigner, dist.Expon = dist.Expon)
  	} 
  
  	if(!is.null(fn)) {   			
    	if(ftype == "png") png(fn) else pdf(fn, width = 8, height = 8)
    	make.plot(histo = histo, title = title, threshold = threshold, dist.Wigner = dist.Wigner, dist.Expon = dist.Expon)  
    	dev.off()
  	}          
  	invisible(fn)
} 





## Scatterplot of empirical eigenvalue spacing and linear fit  
"rm.spacing.scatter" <- 
function(ev.spacing, title = "Scatterplot of eigenvalue spacing", pop.up = TRUE, fn = NULL) 
{
  	if(!is.numeric(ev.spacing)) stop("\n\n  rm.spacing.scatter: 'ev.spacing' must be a numerical vector.\n\n")  
  	if((pop.up == F) && is.null(fn)) stop("\n\n  rm.spacing.scatter: Nothing to do. Enter filename or set 'pop.up = TRUE'\n\n") 
  
  	if(!is.null(fn)) {
    	ftype = rm.get.file.extension(fn) 	
    	if(!ftype %in% c("png", "pdf")) stop("\n\n  rm.spacing.scatter: Output file must have extension 'png' or 'pdf'.\n\n")
  	}   
  
  	## Linear fit
  	index = 1:length(ev.spacing) 
  	linear.fit = lm(ev.spacing ~ index) 
  	intercept = linear.fit$coefficients[1] 
  	if(!(abs(intercept - 1) < 0.01)) cat(paste(" WARNING: intercept is", round(intercept,4), " - should be 1.\n"))
  	slope = linear.fit$coefficients[2]
  	if(!(abs(slope) < 0.001)) cat(paste(" WARNING: slope is", round(slope,4), " - should be 0.\n\n"))
  
  	make.plot <- function(index = index, ev.spacing = ev.spacing, title = title, linear.fit = linear.fit) {
    	plot(index, ev.spacing, xlab = "index", ylab = "spacing", main = title, font.main = 1, col = "blue", cex = 0.8)  
    	abline(linear.fit, col = "red", lwd = 2.0, lty = 2)
    	mtext(side = 3, "Linear fit", col = "red") 
  	}
   
  	if(pop.up) {
    	make.plot(index = index, ev.spacing = ev.spacing, title = title, linear.fit = linear.fit)
  	}
  
  	if(!is.null(fn)) { 
    	if(ftype == "png") png(fn) else pdf(fn, width = 8, height = 8)
    	make.plot(index = index, ev.spacing = ev.spacing, title = title, linear.fit = linear.fit)    
    	dev.off()
  	}
  
  	invisible(fn)  
}  





## Eigenvalue density plot (Wigner semi-circle)
"rm.ev.density" <- 
function(eigenvalues, nr.breaks = 51, min.bw = 0.01, wigner = TRUE, mark.on.x = NULL, title = "Eigenvalue density distribution", pop.up = TRUE, fn = NULL) 
{   
  	if(!is.numeric(eigenvalues)) stop("\n  rm.ev.density: argument 'eigenvalues' must be a numeric vector.\n\n")
  	if(!is.null(mark.on.x) & !is.numeric(mark.on.x)) stop("\n  rm.ev.density: argument 'mark.on.x' must be numeric.\n\n")
  	if((pop.up == F) && is.null(fn)) stop("\n\n  rm.ev.density: Nothing to do. Enter filename or set 'pop.up = T'\n\n") 
  
  	if(!is.null(fn)) {
    	ftype = rm.get.file.extension(fn) 	
    	if(!ftype %in% c("png", "pdf")) stop("\n\n  rm.ev.density: Output file must have extension 'png' or 'pdf'.\n\n")
  	}   
   
  	make.plot <- function(eigenvalues = eigenvalues, nr.breaks = nr.breaks, wigner = wigner, mark.on.x = mark.on.x, title = title) {   	
		dens = density(eigenvalues, bw = "SJ")
		if(dens$bw < min.bw) dens = density(eigenvalues, bw = min.bw)	## repeat if bandwidth too small
  		dx = max(dens$x) - min(dens$x)		
  		dy = max(dens$y) - min(dens$y)		
  		fill = rm.trapez.int(dens$x, dens$y)/(dx*dy)
		br = seq(min(eigenvalues), max(eigenvalues), len = nr.breaks)
		yl = "PDF of eigenvalues"
		hist(eigenvalues, breaks = br, freq = F, col = "darkolivegreen2", main = title, font.main = 1, xlab = "eigenvalue", ylab = yl)	
		lines(dens, col = "red", lwd = 1.5) 	
    	if(!is.null(mark.on.x)) points(x = mark.on.x, y = rep(0, length(mark.on.x)), col = "darkolivegreen", pch = 2) 
    	mtext(side = 3, paste(" Fill factor =", signif(fill,3)), col = "blue") 
    	if(wigner) {
      		curve(wigner.semi.circle, from = -1, to = +1, n = 1001, add = T, type = "l", lty = 2, col = "blue", lwd = 1.2)  
      		legend("topright", c("Density", "Wigner"), col = c("red", "blue"), lty = c(1,2), lwd = 2, cex = 0.8)
    	}	
  	}
  
  	if(pop.up) {
    	make.plot(eigenvalues = eigenvalues, nr.breaks = nr.breaks, wigner = wigner, mark.on.x = mark.on.x , title = title)    
  	}
  
  	if(!is.null(fn)) { 
    	if(ftype == "png") png(fn) else pdf(fn, width = 8, height = 8)
    	make.plot(eigenvalues = eigenvalues, nr.breaks = nr.breaks, wigner = wigner, mark.on.x = mark.on.x , title = title)     
    	dev.off()
  	}
  
  	invisible(fn)    
} 




"rm.get.sparseness" <- 
function(mat) 
{
  	if(!is.matrix(mat)) stop("\n\n  rm.get.sparseness: argument 'mat' must be a matrix.\n\n")
  	nr.zeros = sum(mat == 0)	
  	nr.cells = nrow(mat)*ncol(mat)
  	sparseness = signif(nr.zeros/nr.cells,4)
  	return(sparseness)
}	



"rm.trapez.int" <- 
function(x, y)  
{	
  	if(length(x) != length(y)) stop("\n\n  rm.trapez.int: vectors 'x' and 'y' must have the same length.\n\n")
  	ind = 2:length(x)
  	integral = as.double((x[ind] - x[ind-1]) %*% (y[ind] + y[ind-1]))/2
  	return(integral)
}


"wigner.surmise" <- function(x) pi/2*x*exp(-pi/4*x^2)		

"rm.exp.distrib" <- function(x) exp(-x)			

"wigner.semi.circle" <- function(x) 2/pi*sqrt(1-x^2)	



"rm.get.file.extension" <- 
function(plotnames) 
{
  	dot.pos = regexpr("\\.([[:alnum:]]+)$", plotnames)	
  	extension = unique(ifelse(dot.pos > -1L, substring(plotnames, dot.pos + 1L), ""))
  	return(extension)
}





"kb.distance" <- 
function(histo) 
{
  	if(class(histo) != "histogram") stop("\n\n  kb.distance: 'histo' must be output of 'hist' function.\n\n")
  	observed = histo$density
  	expected.Wigner  = wigner.surmise(histo$mids)
  	expected.Poisson = rm.exp.distrib(histo$mids)    

  	klw = kld(observed, expected.Wigner,  plot = NULL)$dist   # Kullback-Leibler divergence  observed -- Wigner surmise   
  	klp = kld(observed, expected.Poisson, plot = NULL)$dist   # Kullback-Leibler divergence  observed -- Exponential distrib.
  	return(list(klw = klw, klp = klp))
}  





## Kullback-Leibler divergence 
"kld" <- 
function(observed, expected, plot = NULL)     	
{
  	if(!is.numeric(expected)) stop("\n\n kld: Vector 'expected' must be numeric.\n\n") 
  	if(!is.numeric(observed)) stop("\n\n kld: Vector 'observed' must be numeric.\n\n")   
  	if(all(expected == 0)) stop("\n\n kld: All expected values are zero.\n\n") 
  	if(all(observed == 0)) stop("\n\n kld: All observed values are zero.\n\n")  
  	if(any(expected < 0)) stop("\n\n kld: expected frequency below zero.\n\n")
  	if(any(observed < 0)) stop("\n\n kld: observed frequency below zero.\n\n")  
   
  	result = list()
  	result[["observed"]] = observed
  	result[["expected"]] = expected 
 
  	if(!is.null(plot)) { 
		mtxt = "Distribution of observed and expected values"
    	matplot(cbind(observed, expected), col = c("red", "blue"), main = mtxt, font.main = 1, type = "b", pch = c(1,17), ylab = "PDF", xlab = "index")
    	legend("topright", c("observed", "expected"), col = c("red", "blue"), lwd = 2, pch = c(1,17))
    	dev.copy(png, file = plot); dev.off()
    	result[["plot"]] = plot
  	}
  
  	ind = which(observed <= 0)	# The Kullback-Leibler divergence is defined only if observed == 0 implies expected == 0	   
  	if(length(ind) > 0) {		
    	observed = observed[-ind] 
    	expected = expected[-ind] 
  	}
	
  	expected = expected/sum(expected)  	# normalize
  	observed = observed/sum(observed)
	
	distance = sum(observed * log(observed/expected))	# Inf if any(expected == 0) 

   	result[["dist"]] = distance 
  	return(result)        
} 





## Apply threshold to a matrix, the diagonal is not touched if keep.diag = T
"rm.denoise.mat" <- 
function(mat, threshold, keep.diag = TRUE)	    
{
  	if(!is.matrix(mat)) stop("\n  rm.denoise.mat: argument 'mat' must be a matrix\n\n") 
  	nr.nonzeros.before = sum(mat != 0)
  	if(keep.diag) diagon = diag(mat)						    
  	mat[which(abs(mat) < abs(threshold), arr.ind=T)] = 0    	
  	if(keep.diag)  diag(mat) = diagon					
  	nr.nonzeros.after = sum(mat != 0) 
  	cat(paste("\n  Number of non-zero matrix elements reduced from", nr.nonzeros.before, "to", nr.nonzeros.after, "\n\n")) 		 
  	return(mat)
}



## Remove rows and columns consisting of zeros only
"rm.discard.zeros" <- 
function(mat, tol = 0, silent = FALSE) 
{    
  	if(!is.matrix(mat)) stop("\n\n  discard.zero.rows: argument 'mat' must be a natrix.\n\n")
  	is.null.vector <- function(x) ifelse(all(abs(x) <= tol), TRUE, FALSE)
  	diagon = diag(mat)
  	diag(mat) = 0					 
  	zero.rows = apply(mat, 1, is.null.vector)	
  	zero.cols = apply(mat, 2, is.null.vector)	 
  	mat = mat[!zero.rows, !zero.cols]  		
  	diagon = diagon[!zero.rows]
  	diag(mat) = diagon				
  	if(!silent) if((sum(zero.rows) > 0)|(sum(zero.cols) > 0)) cat(paste("  ", sum(zero.rows), "zero-rows and", sum(zero.cols), "zero-columns removed.\n")) 
  	return(mat)
}



"create.rand.mat" <- 
function(size = 1000, distrib = c("norm", "unif"), mean = 0, stddev = 1)	
{ 
  distrib <- match.arg(distrib, c("norm", "unif"))
  cat(paste("\n  Required distribution of matrix elements:", distrib, "\n"))
  
  if((distrib != "norm") & ((mean!=0) | (stddev!=1))) cat("\n  Parameters 'mean' and 'stddev' ignored, apply only for the Normal distribution.\n\n")    
  if(distrib == "norm") data = rnorm(size*size, mean = mean, sd = stddev)	
  if(distrib == "unif") data = runif(size*size, min = -1, max = 1)			
    
  rand.mat = matrix(data, nrow=size) 
  if(distrib == "norm") rand.mat = (rand.mat + t(rand.mat))/sqrt(2)	# make symmetric, normalize 
  if(distrib == "unif") rand.mat = (rand.mat + t(rand.mat))/2		# make symmetric, normalize
  
  m.diag  = mean(diag(rand.mat))	
  sd.diag = sd(diag(rand.mat))		 
  ut = rand.mat[upper.tri(rand.mat)]
  m.ut = mean(ut)					
  sd.ut = sd(ut)					

  cat(paste("  The mean of the main diagonal is", signif(m.diag, 4), "\n"))
  cat(paste("  The std. deviation of the main diagonal is", signif(sd.diag, 4), "\n"))
  cat(paste("  The mean of the upper triangle is", signif(m.ut, 4), "\n"))
  cat(paste("  The std. deviation of the upper triangle is", signif(sd.ut, 4), "\n"))
  cat("  The matrix is real and symmetric.\n\n")
  
  results = list()  
  results[["mean.diag"]] = signif(m.diag, 4)
  results[["stddev.diag"]] = signif(sd.diag, 4)
  results[["mean.triangle"]] = signif(m.ut, 4)
  results[["stddev.triangle"]] = signif(sd.ut, 4)  
  results[["rand.matrix"]] = rand.mat  
  return(results)  
}  





"rm.show.test" <- 
function(thresholds, p.values, main = "P-value vs. threshold", fn = "pval.vs.threshold.png", interactive = TRUE)
{
	t = thresholds; p = p.values
	if(length(t) != length(p))  stop("\n\n  rm.show.test: 'thresholds' and 'p.values' must have the same length.\n\n")	
	if(rm.get.file.extension(fn) != "png") stop("\n\n  rm.show.test: 'fn' must be 'NULL' or some filename with extension '.png'.\n\n")	

	plot(t, p, col = "red", main = main, font.main = 1, ylab = "p.value", xlab = "threshold")
	mtext(side = 3, "Low p-values indicate that exponentiality is unlikely", col = "blue", cex = 0.9)								
	
	if(interactive) {
		cat("\n  Mark proper threshold by left-clicking on red point(s), stop with right-click\n\n")
		tp <- suppressWarnings(sapply(list(t, p),"[",identify(t, p, labels = as.character(signif(t,4))))) 	
		if(is.matrix(tp)) chosen = tp[,1] else chosen = tp[1]
		if(is.list(chosen)) chosen = "not chosen"
	}
		
	if(!is.null(fn)) {dev.copy(png, file = fn); dev.off()}
	
    results <- list()
	results[["chosen"]] = chosen		
	results[["plot"]] = fn 
	return(results)
}







"rm.sse.plot" <- 
function(thresholds, sse.values, main = "SSE for NNSD <--> Exponential", fn = "RMT.SSE.vs.Thres.png", interactive = TRUE)
{
	t = thresholds; sse = sse.values
	if(length(t) != length(sse))  stop("\n\n  rm.sse.plot: 'thresholds' and 'sse.values' must have the same length.\n\n")	
	if(rm.get.file.extension(fn) != "png") stop("\n\n  rm.sse.plot: 'fn' must be 'NULL' or some filename with extension '.png'.\n\n")	

	plot(t, sse, col = "red", main = main, font.main = 1, ylab = "SSE", xlab = "threshold")
	mtext(side = 3, "NNSD is close to Exponential when SSE is low", col = "blue", cex = 0.8)								
	
	if(interactive) {
		cat("\n  Mark proper threshold by left-clicking on red point(s), stop with right-click\n\n")
		tp <- suppressWarnings(sapply(list(t, sse),"[",identify(t, sse, labels = as.character(signif(t,4))))) 	
		if(is.matrix(tp)) chosen = tp[,1] else chosen = tp[1]
		if(is.list(chosen)) chosen = "not chosen"
	}
		
	if(!is.null(fn)) {dev.copy(png, file = fn); dev.off()}
	
    results <- list()
	results[["chosen"]] = chosen		
	results[["plot"]] = fn 
	return(results)
}




"rm.show.plots" <- 
function(plotnames) 
{
  	if(!is.character(plotnames)) stop("\n rm.show.plots: argument 'plotnames' must be a character string or -vector containing the names of the plots.\n\n")  
  	cat(paste("\n  Current folder is:", getwd(), "\n\n"))
    
  	ext = rm.get.file.extension(plotnames)  
  	if((length(ext) > 1) || (ext != "png")) stop("\n\n rm.show.plots: Please use png format only.\n\n")
  
  	for (i in 1:length(plotnames)) {
    	file = plotnames[i]					
    	if(!file.exists(file)) {
      		cat(paste("  File '", file, "' not found.\n"))
      		next
    	} else {
      		img = png::readPNG(file)				
      		plot.new()					
      		usr <- par("usr")					
      		rasterImage(img, usr[1], usr[3], usr[2], usr[4])	
      		cat ("Press [enter] to continue")
      		line <- readline() 
    	}
  	} 
  	cat(paste("\n", i, "plots shown\n\n"))
} 



"add.Gaussian.noise" <- 
function(mat, mean = 0, stddev = 1, symm = TRUE)  
{
  	if(!is.matrix(mat)) stop("\n\n  add.Gaussian.noise: argument 'mat' must be a matrix.\n\n")
  	cat(paste("\n ", nrow(mat), "x", ncol(mat), "matrix read.\n"))
  	if((symm == TRUE) && (nrow(mat) != ncol(mat))) stop("\n\n  add.Gaussian.noise: matrix not quadratic - cannot be symmetrized.\n\n")
  
  	nr.nonzero = sum(mat != 0)
  	cat(paste("  This matrix contains", nr.nonzero, "non-zero cells.\n\n"))
  
  	noise = matrix(rnorm(nrow(mat)*ncol(mat), mean = mean, sd = stddev), nrow=nrow(mat))
  	mat = mat + noise 
  	if(symm) mat = (mat + t(mat))/sqrt(2)	   
  
  	cat(paste("  Gaussian noise with mean", mean, "and standard deviation", stddev, "added.\n"))
  	cat(paste("  The noise ranges from", signif(min(noise),4), "to", signif(max(noise),4), "\n"))
  	if(symm) cat("  The output matrix has been symmetrized.\n\n") else cat("\n")
  	return(mat)
}






## Reorder eigenvalues and eigenvectors accordingly
"rm.reorder.ev" <- 
function(eigenvalues, eigenvec = NULL) 	
{
  if(!is.null(eigenvec)) {
    if(!is.matrix(eigenvec)) stop("\n  rm.reorder.ev: 'eigenvec' must be a matrix (eigenvectors in columns)\n\n")
    if(ncol(eigenvec) != length(eigenvalues))  stop("\n  rm.reorder.ev: number of cols of 'eigenvec' must be equal to length of 'eigenvalues'\n\n")
  }    
  new.order = order(eigenvalues)		# smallest first
  if(!is.null(eigenvec)) eigenvec = eigenvec[, new.order]
  eigenvalues = eigenvalues[new.order] 
  results = list()
  results[["eigenvalues"]] = eigenvalues
  if(!is.null(eigenvec)) results[["eigenvec"]] = eigenvec
  return(results)
} 




"rm.connections" <- 
function(mat, nr.list = 30, abs.val = TRUE, fn = NULL)  
{
  	if(!is.matrix(mat)) stop("\n\n  rm.connections: 'mat' must be a matrix.\n\n") 
  	if(isSymmetric(mat)) cat("\n  Matrix is symmetric.\n\n") else cat("\n  Matrix is not symmetric.\n\n")   
  	if(isSymmetric(mat)) 
    	if(nr.list > sum(upper.tri(mat))) stop("\n\n  rm.connections: desired number of list elements ('nr.list') exceeds number of unique matrix elements.\n\n") 
  	if(!isSymmetric(mat)) 
    	if(nr.list > nrow(mat)*ncol(mat)) stop("\n\n  rm.connections: desired number of list elements ('nr.list') exceeds number of matrix elements.\n\n")
 	if(!is.null(fn)) 
    	if(rm.get.file.extension(fn) != "txt") stop("\n\n  rm.connections: output file must have extension 'txt' (can only be saved as text file).\n\n")
  
  	if(abs.val) temp.mat = abs(mat) else temp.mat = mat   
  	if(isSymmetric(temp.mat)) {
    	v = temp.mat[upper.tri(temp.mat)]
    	value = mat[upper.tri(mat)]			
    	rows = integer(0); for (i in 1:(nrow(temp.mat)-1)) rows = c(rows, 1:i)
    	cols = integer(0); for (i in 2:ncol(temp.mat)) cols = c(cols, rep(i,i-1))
  	} else {
    	v = as.vector(temp.mat)  
    	value = as.vector(mat)			
    	rows = rep(1:nrow(temp.mat), ncol(temp.mat))
    	cols = rep(1:ncol(temp.mat), each = nrow(temp.mat))
  	}
  	df = data.frame(row = rows, col = cols, v = v, value = value)  
  	df = df[order(df$v, decreasing = T),]
  	df$v <- NULL
  	df = df[1:nr.list,]
  	rownames(df) = 1:nrow(df)
  	if(!is.null(rownames(mat))) df = cbind(df, row.name = rownames(mat)[df$row])
  	if(!is.null(colnames(mat))) df = cbind(df, col.name = colnames(mat)[df$col])       
  
  	print(df) 
  	if(!is.null(fn)) write.table(df, file = fn, quote = F, row.names = F, sep = "\t") 
  	return(df)   
}


## uwemenzel@gmail.com




















