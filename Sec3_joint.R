################################################################################
# Sec3_joint.R                                                                 #
#    Simulate from joint posterior for PriorPC and PriorJe                     #
################################################################################

# Set to TRUE to save .eps figures instead of drawing on screen
do.print <- TRUE

# Graphics settings
bty = "l"
lwd = 2L
cex.lab = 1.4
cex.axis = 1.4

# Load stored realization
load('locations.RData')
load('observation.RData')

# Plot spatial locations
if (do.print) {
    setEPS()
    postscript("Figures/spatialDesign.eps")
}
plot(loc, xlim = c(0,1), ylim=c(0,1), xlab = "", ylab = "", cex.axis = cex.axis, bty = bty, lwd = lwd)
if (do.print) {
    dev.off()
}

# Load MCMC sampler
source('MCMCsampler.R')

# Set random seed
set.seed(100)

# Jeffreys' rule prior
numSamples = 1000000
modelObj = intPrior.Exp.makeModel(obs = obs, loc = loc)
modelObj = intPrior.Exp.addPrior(modelObj = modelObj, name = "Jeffrey")
modelObj = intPrior.Exp.addTarget(modelObj)
inits = c(0, 0)
scale = c(0.6, 0.08)
acc.rate = 0.4
adapt = TRUE
chainJe = intPrior.addHyperSamples(modelObj = modelObj, numSamples = numSamples, inits = inits, scale = scale, acc.rate = acc.rate, adapt = adapt)

# Remove burn-in and thinning
chainJe$hyper$samples = chainJe$hyper$samples[seq(25001, numSamples, 1),]
chainJe$hyper$samples = chainJe$hyper$samples[floor(seq(1, numSamples-25000, length.out = 10000)),]
        
# PC prior
numSamples = 250000
a1  = 0.05
rho0 = 0.1
a2 = 0.05
sig0 = 10
modelObj = intPrior.Exp.makeModel(obs = obs, loc = loc)
pars = list(a1 = a1, rho0 = rho0, a2 = a2, sig0 = sig0)
modelObj = intPrior.Exp.addPrior(modelObj = modelObj, name = "PC-prior", pars = pars)
modelObj = intPrior.Exp.addTarget(modelObj)
inits = c(0, 0)
scale = c(0.7, 0.6)
acc.rate = 0.4
adapt = TRUE
chainPC = intPrior.addHyperSamples(modelObj = modelObj, numSamples = numSamples, inits = inits, scale = scale, acc.rate = acc.rate, adapt = adapt)

# Remove burn-in and thin
chainPC$hyper$samples = chainPC$hyper$samples[seq(5001, numSamples, 1),]
chainPC$hyper$samples = chainPC$hyper$samples[floor(seq(1, numSamples - 5000, length.out = 10000)),]

 # Bivariate plot in log-scale
if (do.print) {
    postscript('Figures/bivariateLog.eps')
}
par(cex.lab=cex.lab, cex.axis=cex.axis)
idxJe = floor(seq(1, dim(chainJe$hyper$samples)[1], length.out = 1000))
plot(exp(as.matrix(chainJe$hyper$samples[idxJe, ])), xlab = "Range", ylab = "Standard deviation", log = "xy", lwd = lwd, bty = bty)
idxPC = floor(seq(1, dim(chainPC$hyper$samples)[1], length.out = 1000))
points(exp(as.matrix(chainPC$hyper$samples[idxPC, ])), col = 'red', lwd = lwd)
if (do.print) {
    dev.off()
}

# Marginal posterior for log(range)
densJeff = density(chainJe$hyper$samples[,1])
limJeff = quantile(chainJe$hyper$samples[,1], probs = c(0.025, 0.975))
densPC   = density(chainPC$hyper$samples[,1])
limPC = quantile(chainPC$hyper$samples[,1], probs = c(0.025, 0.975))
if (do.print) {
    setEPS()
    postscript('Figures/logNomRange.eps')
}
plot(densJeff$x, densJeff$y, xlab = 'Log(Range)', ylab = 'Density', type = 'l', xlim = c(-3, 12), ylim = c(0, 0.8), bty = bty, lwd = lwd)
abline(v = limJeff[1], lwd = lwd)
abline(v = limJeff[2], lwd = lwd)
lines(densPC$x, densPC$y, lty = 2, lwd = lwd)
abline(v = limPC[1], lty = 2, lwd = lwd)
abline(v = limPC[2], lty = 2, lwd = lwd)
if (do.print) {
    dev.off()   
}

# Marginal posterior for log(std.dev.)
densJeff = density(chainJe$hyper$samples[,2])
limJeff = quantile(chainJe$hyper$samples[,2], probs = c(0.025, 0.975))
densPC   = density(chainPC$hyper$samples[,2])
limPC = quantile(chainPC$hyper$samples[,2], probs = c(0.025, 0.975))
if (do.print){
    setEPS()
    postscript('Figures/logStdDev.eps')
}
plot(densJeff$x, densJeff$y, xlab = 'Log(Standard deviation)', ylab = 'Density', type = 'l', xlim = c(-2, 6), ylim = c(0, 1.5), bty = bty, lwd = lwd)
abline(v = limJeff[1], lwd = lwd)
abline(v = limJeff[2], lwd = lwd)
lines(densPC$x, densPC$y, lty = 2, lwd = lwd)
abline(v = limPC[1], lty = 2, lwd = lwd)
abline(v = limPC[2], lty = 2, lwd = lwd)
if (do.print){
    dev.off()

}
