################################################################################
# SecS1.R                                                                      #
#    Bounded domain examples                                                   #
################################################################################

# Load libraries
library(INLA)

# Construct Matern covariance matrix with tau fixed
makeCov = function(kappa, loc, nu){
    d = dim(loc)[2]
    Sig = inla.matern.cov(nu = nu, kappa = kappa, x = as.matrix(dist(loc)), d = d)
    return(Sig)
}

# Calculate KLD  
gaussKLD = function(Sig, Sig0){
    tmpMat = solve(Sig0, Sig)
    KLD = 0.5*(sum(diag(tmpMat)) - dim(Sig0)[1] - determinant(tmpMat)$modulus)
    
    return(KLD)
}

# Calculate prior wrt base model
calcKappaPi = function(kappa, loc, nu, rho0, lam = 1, eps = 1e-3){
    # Calculate covariance matrix of base model
    kappa0 = sqrt(8*nu)/rho0
    Sig0 = makeCov(kappa0, loc, nu)
    
    # Need derivative of distance
    kappaS = kappa + eps
    
    # Iterate through different kappa's
    KLD = kappa
    KLDS = kappaS
    for(kIdx in 1:length(kappa)){
        print(kIdx)
        
        # Calculate covariance matrix of alternative model
        Sig = makeCov(kappa[kIdx], loc, nu)
        KLD[kIdx] = gaussKLD(Sig, Sig0)
        
        # Calculate covariance matrix of shifted kappas
        Sig = makeCov(kappaS[kIdx], loc, nu)
        KLDS[kIdx] = gaussKLD(Sig, Sig0)
    }
    
    # Calculate distance
    dd = sqrt(2*KLD)
    ddS = sqrt(2*KLDS)
    
    # Calculate derivative
    dFac = (ddS-dd)/(kappaS-kappa)
    
    # Calculate numerical prior
    kVec = kappa
    piKap = lam*exp(-lam*dd)*dFac
    
    return(list(piKappa = piKap, dd = dd, dFac = dFac, kappa = kVec))
}

# Set to TRUE to save .eps figures instead of drawing on screen
do.print <- TRUE

## Exact prior calculation for d = 1, \nu = 0.5
    d = 1
    nu = 0.5
    
    # Base model
    Rs = c(0.5, 1, 2, 4)
    YLIM = c(2.2, 2.2, 1.6, 1.6)
    for(i in 1:4){
        
        rho0 = Rs[i]
        if(do.print){
            postscript(paste("Figures/1DExpExact_", i, ".eps", sep =""))
        }
        plot(NULL, xlim = c(0, rho0), ylim = c(0, YLIM[i]), xlab = "rho", ylab = "Probability density", bty = "l", cex.axis = 1.5, cex.lab = 1.5, cex = 1.5)
        
        # Analytic expression for distance as a function of range
        dFun = function(rho){
            ans = sqrt(rho/rho0 - 1 - log(rho) + log(rho0) + sqrt(8*nu)*(rho/(2*rho0^2) - 1/rho0 + 1/(2*rho)))
        }
        
        # Analytic expression for the derivative of distance as a function of range
        derFun = function(rho){
            ans = 1/(2*dFun(rho))*(1/rho - 1/rho0 + sqrt(8*nu)*(1/(2*rho^2) - 1/(2*rho0^2)))
        }
        
        # Grid for evaluating prior
        #rho = exp(seq(-10, log(rho0)-1e-3, length.out = 100))
        rho = seq(exp(-10), rho0-1e-3, length.out = 1000)
        
        # Calculate distances
        dd = dFun(rho)
        
        # Calculate scaling factors
        dFac = abs(derFun(rho))
        
        # Prior specification
        alpha = 0.05
        R0 = 0.05
        
        # Calculate lamda to achieve desired calibration
        lam = -log(alpha)/dFun(R0)
        piR = lam*exp(-lam*dd)*dFac
        lines(rho, piR, lwd = 2)
        
        # Calculate asymptotic expression limited to same domain
        #     lam2 = -log(alpha)/(R0^-0.5-rho0^-0.5)
        #        piAR = d/2*lam2*rho^-1.5*exp(-lam2*rho^-0.5)/exp(-lam2*rho0^-0.5)
        #       lines(rho, piAR, col = "red", lwd = 1.5, lty = 2)
        
        # Calculate unlimited asymptotic
        lam3 = -log(alpha)/R0^-0.5
        piFR = d/2*lam3*rho^-1.5*exp(-lam3*rho^-0.5)
        lines(rho, piFR, lwd = 2, lty = 2)
        if(do.print){
            dev.off()
        }
    }

### Approximate calculations for d = 2 and \nu = 1.5
# Base model
nu = 1.5
d = 2

# grid for calculation of prior
rho0 = 4
lrVec = seq(log(1e-5), log(rho0)-1e-4, length.out = 100)
rVec = exp(lrVec)
kappa = sqrt(8*nu)/rVec
kappa0 = kappa[1]

# Set calibration of prior
R0 = 0.05
alpha = 0.05

# Iterate through different number of observation locations
numX = c(10, 20, 40)
piR = matrix(0, nrow = length(rVec), ncol = length(numX))
ddR = piR
for(i in 1:length(numX)){
    # Observation locations
    x = seq(0, 1, length.out = numX[i])
    loc = matrix(c(rep(x, numX[i]), rep(x, each = numX[i])), ncol = 2)
    
    # Calculate lambda
    tmpPrior = calcKappaPi(sqrt(8*nu)/R0, loc, nu, rho0)
    tmpDist = tmpPrior$dd
    lam = -log(alpha)/tmpDist
    
    # Calculate prior
    piKappaStr = calcKappaPi(kappa, loc, nu, rho0, lam)
    
    # Extract prior for kappa
    piKap = piKappaStr$piKappa
    ddR[, i] = piKappaStr$dd
    
    # Calculate prior for rho
    piR[, i] = piKap*sqrt(8*nu)/rVec^2
}

# Calculate unlimited asymptotic
d = 2
lam3 = -log(alpha)/R0^-(d/2)
piFR = d/2*lam3*rVec^(-d/2-1)*exp(-lam3*rVec^-(d/2))

# Make plot of in original scale
if(do.print){
    postscript("Figures/boundedOrig.eps")
}
plot(NULL, xlim = c(1e-6, rho0), ylim = c(1e-5, 4.0), xlab = "rho", ylab = "Probability density", bty = "l", cex = 1.5, cex.lab = 1.5, cex.axis = 1.5)
for(i in 1:3)
    lines(rVec, piR[,i], lwd = 2, lty = i)
lines(rVec, piFR, lwd = 2, lty = 4)
if(do.print){
    dev.off()
}

# Make plot in log-log scale
if(do.print){
    postscript("Figures/boundedLogLog.eps")
}
plot(NULL, xlim = c(1e-6, rho0), ylim = c(1e-5, 25.0), xlab = "rho", ylab = "Probability density", bty = "l", cex = 1.5, cex.lab = 1.5, cex.axis = 1.5, log = "xy")
for(i in 1:3)
    lines(rVec, piR[,i], lwd = 2, lty = i)
lines(rVec, piFR, lwd = 2, lty = 4)
if(do.print){
    dev.off()
}
