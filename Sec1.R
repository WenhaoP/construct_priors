################################################################################
# Sec1.R                                                                       #
#    Code for Section 1 of the paper                                           #
################################################################################

# Set to TRUE to save .eps figures instead of drawing on screen
do.print <- TRUE

# Set seed
set.seed(10)

# Choose locations
n = 100
loc = runif(n)

# Initialize storage
x = matrix(0, nrow = n, ncol = 3)  

# Generate iid realization
z = rnorm(n)

# Generate realizations with structure
for(i in 1:3){
    # Set parameters
	rho = 10^((i-1)*3)
	sig2 = 10^((i-1)*3)

    # Calculate covariance matrix
    D = as.matrix(dist(loc))
	Sigma = sig2*exp(-D/rho)

    # Generate realization
	L = chol(Sigma)
	x[,i] = t(L)%*%z
}

# Find appropriate y limits
yLims = apply(x, 2, range)
yRs = diff(yLims)
yR = max(yRs)

# Make plots
for(i in 1:3){
    if (do.print) {
        setEPS()
        postscript(paste("Figures/sampExp_", i, ".eps", sep = ""))
    }
    yMean = mean(yLims[,i])
    yLim = c(yMean-0.7*yR, yMean+0.7*yR)
    plot(loc, x[,i], bty = "l", lwd = 3L, cex.lab = 3, cex.axis = 3, cex =1.5, xlab = "", ylab = "", ylim = yLim)
    
    if(do.print) {
        dev.off()
    }
}
