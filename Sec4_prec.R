################################################################################
# Sec4_prec.R                                                                  #
#    Analyze precipitation                                                     #
################################################################################

### Initialize
    ## Library
        # Load libraries
        library(MASS)
#        library(mixOmics)
        library(PBSmapping)
        library(INLA)
        source('deformationModel.R')
    
    ## Directory to store results
        bDir = "Results/Precipitation/"
    
    ## Data
        # Load data
        load("data.RData")
        
        # Re-scale observation	
        prec = data$precipitation/1000
        
        # Observation data
        obs = prec
        loc = cbind(data$x, data$y)
        
        # Get covariates from digital elevation model (GLOBE)
        linCov = cbind(1, data$altitude/1e3)
        source('readHeight.R')
        tmpCov = getCovariates(mesh)
        covHeight = tmpCov$covHeight
        covHeight = covHeight/1e3
        covGrad = tmpCov$covGrad
        covGrad = covGrad/1e2
        covCov = cbind(as.vector(covHeight), as.vector(covGrad))
#       covCov = matrix(as.vector(covHeight), ncol = 1)
        covCov = inla.spde3.deform.centerCovariates(Cov = covCov, mesh = mesh)

### Stationary model
    ## Set-up the model
        # Set observations and covariates
        statModel = inla.spde3.deform.makeModel(obs = obs,
                                                loc = loc,
                                                mesh = mesh,
                                                linCov = linCov,
                                                covRange = NULL,
                                                covStd = NULL)
        
        # Add prior
        statModel = inla.spde3.deform.setPrior(modelObj = statModel,
                                               uN = 3, aN = 0.05,
                                               lR = 10, aR = 0.05,
                                               uS = 3, aS = 0.05)
        
        # Add constraints
        constr = cbind(1, covCov[,1])
        statModel = inla.spde3.deform.setTarget(modelObj = statModel,
                                                constr = constr)
    ## INLA
        # Fit
        guess = c(4.09, 5.389, -0.322)
        result.stat = inla.spde3.deform.analyze(modelObj = statModel,
                                                control.mode = list(theta = guess, restart = TRUE))
    ## MCMC
        # Sample hyperparameters
        numSamples = 300000
        load('statScale.RData')
        scale = scaleStat
        acc.rate = 0.4
        adapt = TRUE
        inits = as.vector(result.stat$mode$theta)
        mcmc.stat = mcmc.addHyperSamples(statModel, numSamples, inits, scale, acc.rate, adapt)
        
        # Remove burn-in
        mcmc.stat$hyper$samples = mcmc.stat$hyper$samples[-(1:5e4),]
        mcmc.stat$hyper$log.p   = mcmc.stat$hyper$log.p[-(1:5e4)]
        
        # Thin
        idxThin = floor(seq(1, length(mcmc.stat$hyper$log.p), length.out = 15000))
        mcmc.stat$hyper$samples = mcmc.stat$hyper$samples[idxThin, ]
        mcmc.stat$hyper$log.p   = mcmc.stat$hyper$log.p[idxThin]
        
    ## Plotting
        # Plot settings
        bty = "l"
        lwd = 2L
        cex.lab = 1.4
        cex.axis = 1.4
        new.windows = TRUE
        write.new.figures = FALSE
        
        # Noise
        if(new.windows)
            inla.dev.new()
        inla.spde3.deform.plotPosterior(result.stat, 1, function(x){1/sqrt(x)}, xlab = "Nugget standard deviation", ylab = "Density")
        if (write.new.figures) {
            dev.print(postscript,  file=paste(bDir, "statNugget.eps", sep = ""))
        }
        
        # Spatial Range
        if(new.windows)
            inla.dev.new()
        inla.spde3.deform.plotPosterior(result.stat, 2, function(x){exp(x)},  xlab = "Range", ylab = "Density")
        if (write.new.figures) {
            dev.print(postscript,  file=paste(bDir, "statRange.eps", sep = ""))
        }
        
        # Noise
        if(new.windows)
            inla.dev.new()
        inla.spde3.deform.plotPosterior(result.stat, 3, function(x){exp(x)}, xlab = "Nugget standard deviation", ylab = "Density")
        if (write.new.figures) {
            dev.print(postscript,  file=paste(bDir, "statStdDev.eps", sep = ""))
        }
        
        # Estimated parameters
        cat("Stationary model:\n")
        cat(paste("\tNugget std.dev. : ", exp(-0.5*result.stat$mode$theta[1]), "\n", sep = ""))
        cat(paste("\tRange:          : ", exp(result.stat$mode$theta[2]), "\n", sep = ""))
        cat(paste("\tProcess std.dev : ", exp(result.stat$mode$theta[3]), "\n\n", sep = ""))
        
    ## Model evaluation
        # Re-calculate with INLA returning CPO and WAIC
        init = result.stat$mode$theta
        result.stat.pred = inla.spde3.deform.analyze(statModel,
                                                     compute = TRUE,
                                                     control.mode = list(theta = init, restart = TRUE), 
                                                     control.compute = list(cpo = TRUE,
                                                                            dic = TRUE,
                                                                            mlik = TRUE,
                                                                            waic = TRUE,
                                                                            po = TRUE))
        result.stat.cpo = inla.cpo(result.stat.pred)
        
        # Sample from predictive distributions
        mcmc.stat = mcmc.addPredSamples(modelObj = mcmc.stat, loc = loc, covValues = linCov, thinning = 1)
        
        # Calculate predictive scores (WAIC and hold-one-out log-predictive scores)
        mcmc.stat = mcmc.LOOP(mcmc.stat, WAIC = TRUE, CRPS = TRUE)
  
        # Compare MCMC and INLA
        plot(mcmc.stat$pred$LOOLpost, log(result.stat.cpo$cpo$cpo))
        
### Non-stationary model
    ## Set-up the model
        # Set observations and covariates
        nonStatModel = inla.spde3.deform.makeModel(obs = obs,
                                                   loc = loc,
                                                   mesh = mesh,
                                                   linCov = linCov,
                                                   covRange = covCov,
                                                   covStd = covCov)
        
        # Add prior
        lamR = 20
        lamS = 20
        nonStatModel = inla.spde3.deform.setPrior(modelObj = nonStatModel, 
                                                  uN = 3, aN = 0.05,
                                                  lR = 10, aR = 0.05,
                                                  uS = 3, aS = 0.05,
                                                  lamRange = lamR, lamStdDev =  lamS)	
        print(mcmc.priorProb(lamR, 100000, nonStatModel)) # Check lambda value
        
        # Add constraints to spatial field
        nonStatModel = inla.spde3.deform.setTarget(modelObj = nonStatModel, constr = constr)

    ## INLA
        # Fit
        guess = c(4.26, 5.79, -0.5, 3.49, 6.06, -2.85, 3.00, 0, 0)
        # NOTE: THE NON-STATIONARY MODEL REQUIRES A SPECIAL VERSION OF INLA
        result.ns = inla.spde3.deform.analyze(nonStatModel,
                                              control.mode = list(theta = guess, restart = TRUE))
    
        # Use sampling to get posteriors for non-linear combinations
        theta = inla.hyperpar.sample(n = 100000, result.ns, intern = TRUE)
        theta[,5:6] = theta[,5:6]*exp(-0.5*theta[,4])
        theta[,8:9] = theta[,8:9]*exp(-0.5*theta[,7])
        
        # Store results
        t = Sys.time()
        tmpStmp = strftime(t, "%Y-%m-%d-%H-%M-%S")
        fName = paste(bDir, 'inlaFit_time_', tmpStmp, '_lamRange_', lamR, '_lamStdDev_', lamS, '.RData', sep = "")
        save(result.ns, theta, file = fName)
        
    ## MCMC
        # Sample hyperparameters
        numSamples = 600000
        load('scaleNS.RData')
        scale = scaleNS
        acc.rate = 0.4
        adapt = TRUE
        inits = result.ns$mode$theta
        mcmc.ns = mcmc.addHyperSamples(nonStatModel, numSamples, inits, scale, acc.rate, adapt)
        
        # Remove burn-in
        mcmc.ns$hyper$samples = mcmc.ns$hyper$samples[-(1:5e4),]
        mcmc.ns$hyper$log.p   = mcmc.ns$hyper$log.p[-(1:5e4)]
        
        # Thin
        idxThin = floor(seq(1, length(mcmc.ns$hyper$log.p), length.out = 7000))
        mcmc.ns$hyper$samples = mcmc.ns$hyper$samples[idxThin, ]
        mcmc.ns$hyper$log.p   = mcmc.ns$hyper$samples[idxThin]
        
        # Extract samples
        samp.mcmc = mcmc.ns$hyper$samples
    
    ## Plotting
        # Plot prior and posterior
        pr.spde = function(x, theta, S){
            int.fun = function(x, theta, S){
                    # hyperprior
                    lTau = x
                    lam = 20
                    log.prior = inla.pc.dprec(prec = exp(lTau), lam = lam, log = TRUE) + lTau
                    
                    # Prior
                    S = S*exp(0.5*lTau)
                    log.prior = log.prior + inla.spde3.deform.logGauss(theta, S)
                    
                    exp(log.prior)
            }
            
            ths = seq(-2,2, length.out = 1000)
            val = ths
            for(i in 1:length(ths)){
                val[i] = integrate(Vectorize(int.fun), lower = 0, upper = 30, theta = ths[i], S = S)$value
            }
            
        }
        
        # Posteriors for hyperparameters in range INLA
        if(new.windows)
            inla.dev.new()
        par(cex.lab = cex.lab, cex.axis = cex.axis)
        dPar1 = hist(theta[,5], 100, plot = FALSE)
        dPar2 = hist(theta[,6], 100, plot = FALSE)
        xLim1 = range(c(dPar1$breaks, dPar2$breaks))
        yLim1 = range(c(dPar1$density, dPar2$density))
        yLim1[1] = 0
        plot(dPar1, freq = F, col = rgb(0.8, 0, 0, 0.25), xlim = xLim1, ylim = yLim1, main = "", xlab = "Parameter")
        plot(dPar2, freq = F, col = rgb(0, 0, 0.8, 0.25), add = TRUE)
        if (write.new.figures) {
            dev.copy2eps(file=paste(bDir, "nonStatRangeINLA.eps", sep = ""))
        }
        
        # Posteriors for hyperparameters in std.dev. INLA
        if(new.windows)
            inla.dev.new()
        par(cex.lab = cex.lab, cex.axis = cex.axis)
        dPar1 = hist(theta[,8], 100, plot = FALSE)
        dPar2 = hist(theta[,9], 100, plot = FALSE)
        xLim2 = range(c(dPar1$breaks, dPar2$breaks))
        yLim2 = range(c(dPar1$density, dPar2$density))
        yLim2[1] = 0
        plot(dPar1, freq = F, col = rgb(0.8, 0, 0, 0.25), xlim = xLim2, ylim = yLim2, main = "", xlab = "Parameter")
        plot(dPar2, freq = F, col = rgb(0, 0, 0.8, 0.25), add = TRUE)
        if (write.new.figures) {
            dev.copy2eps(file=paste(bDir, "nonStatStdDevINLA.eps", sep = ""))
        }
        
        # Posteriors for hyperparameters in range MCMC
        if(new.windows)
            inla.dev.new()
        par(cex.lab = cex.lab, cex.axis = cex.axis)
        dPar1 = hist(samp.mcmc[,5], 100, plot = FALSE)
        dPar2 = hist(samp.mcmc[,6], 100, plot = FALSE)
#         xLim = range(c(dPar1$breaks, dPar2$breaks))
#         yLim = range(c(dPar1$density, dPar2$density))
#        yLim[1] = 0
        plot(dPar1, freq = F, col = rgb(0.8, 0, 0, 0.25), xlim = xLim1, ylim = yLim1, main = "", xlab = "Parameter")
        plot(dPar2, freq = F, col = rgb(0, 0, 0.8, 0.25), add = TRUE)
        if (write.new.figures) {
            dev.copy2eps(file=paste(bDir, "nonStatRangeMCMC.eps", sep = ""))
        }
        
        # Posteriors for hyperparameters in std.dev. MCMC
        if(new.windows)
            inla.dev.new()
        par(cex.lab = cex.lab, cex.axis = cex.axis)
        dPar1 = hist(samp.mcmc[,8], 100, plot = FALSE)
        dPar2 = hist(samp.mcmc[,9], 100, plot = FALSE)
#         xLim = range(c(dPar1$breaks, dPar2$breaks))
#         yLim = range(c(dPar1$density, dPar2$density))
#        yLim[1] = 0
        plot(dPar1, freq = F, col = rgb(0.8, 0, 0, 0.25), xlim = xLim2, ylim = yLim2, main = "", xlab = "Parameter")
        plot(dPar2, freq = F, col = rgb(0, 0, 0.8, 0.25), add = TRUE)
        if (write.new.figures) {
            dev.copy2eps(file=paste(bDir, "nonStatStdDevMCMC.eps", sep = ""))
        }

    ## Model evaluation
        # Re-run INLA with calculation of CPO and WAIC
        init = result.ns$mode$theta
        # NOTE: THE NON-STATIONARY MODEL REQUIRES A SPECIAL VERSION OF INLA
        result.ns.pred = inla.spde3.deform.analyze(modelObj = nonStatModel,
                                                   compute = TRUE,
                                                   control.mode = list(theta = init, restart = TRUE), 
                                                   control.compute = list(cpo = TRUE,
                                                                          dic = TRUE,
                                                                          mlik = TRUE,
                                                                          waic = TRUE,
                                                                          po = TRUE))
        result.ns.cpo = try(inla.cpo(result.ns.pred))
        
        # Find predictive distributions
        mcmc.ns = mcmc.addPredSamples(modelObj = mcmc.ns, loc = loc, covValues = linCov, thinning = 1)
        
        # Calculate predictive scores (WAIC and hold-one-out log-predictive scores)
        mcmc.ns = mcmc.LOOP(mcmc.ns, WAIC = TRUE, CRPS = TRUE)

### Compare predictions
    ## Leave-one-out prediction intervals (Gaussian approx)
        stat.cent = mcmc.stat$pred$PP
        stat.up   = stat.cent + 2*mcmc.stat$pred$Psig
        stat.lo   = stat.cent - 2*mcmc.stat$pred$Psig
        
        ns.cent = mcmc.ns$pred$PP
        ns.up   = ns.cent + 2*mcmc.ns$pred$Psig
        ns.lo   = ns.cent - 2*mcmc.ns$pred$Psig
        
        # Plot limits
        ylim = range(c(stat.up, stat.lo, ns.up, ns.lo))
        plot(stat.cent, type = 'p', ylim = ylim)
        for(i in 1:length(ns.up)){
            lines(c(i, i), c(stat.lo[i], stat.up[i]), col = 'green')
        }
        
        plot(ns.cent, type = 'p')
        for(i in 1:length(ns.up)){
            lines(c(i, i), c(ns.lo[i], ns.up[i]), col = 'green')
        }
    
    ## Point-wise leave-one-out
        # Log-predictive
        plot(mcmc.stat$pred$LOOLpost, mcmc.ns$pred$LOOLpost, xlim = c(-8, 1.2), ylim = c(-8, 1.2))
        lines(c(-10, 2), c(-10, 2), col = 'red')        
        print(mean(mcmc.stat$pred$LOOLpost))
        print(mean(mcmc.ns$pred$LOOLpost))
        dev.copy2eps(file=paste(bDir, "nonStatPredLOG.eps", sep = ""))
    
        # CRPS
        plot(mcmc.stat$pred$CRPS, mcmc.ns$pred$CRPS, xlim = c(0, 0.65), ylim = c(0, 0.65))
        lines(c(0, 0.7), c(0, 0.7), col = 'red')
        print(mean(mcmc.stat$pred$CRPS))
        print(mean(mcmc.ns$pred$CRPS))
        dev.copy2eps(file=paste(bDir, "nonStatPredCRPS.eps", sep = ""))
        
### Checking the implied structure
    ## Select test locations
    idx = c(365, 
            936,
            120,
            316, 
            385)
    
    ## Initialize storage
    covMat = matrix(0, nrow = mesh$n, ncol = length(idx))
    sdMat  = matrix(0, nrow = mesh$n, ncol = 1)
    corMat = covMat
    
    ## Average our thetas
    numThetas = 100
    cLength = dim(mcmc.ns$hyper$samples)[1]
    tIdx = floor(seq(1, cLength, length.out = numThetas))
    for(iter in tIdx){
        # Current iteration
        print(iter)
        
        # Get theta
        theta = mcmc.ns$hyper$samples[iter, -1]
        
        # Get unrestricted precision matrix
        Qxx = INLA:::inla.spde3.precision(spde = mcmc.ns$data$spde, theta = theta)
        
        # Factorize
        L = Cholesky(Qxx)
        
        # Marginal standard devations
        sdUn = sqrt(diag(inla.qinv(Qxx)))
        
        ## Precomputation
            # Covariance matrix of constraints
            A = mcmc.ns$constr[, -c(1236,1237)]
            
            # TMP
            constr.Sigma = A%*%solve(L, t(A))
            
            # Transfer
            tr = solve(L, t(A))
        
        # Correct sd according to constraints
        tmpSD = sqrt(sdUn^2 - diag(tr%*%solve(constr.Sigma, t(tr))))
        sdMat = sdMat + tmpSD
        
        # Run through each location
        for(iLoc in 1:length(idx)){
            # Unconstrained covariance
            rhs = rep(0, mesh$n)
            rhs[idx[iLoc]] = 1
            covUn = solve(L, rhs)
            
            # Unconstrained correlation
            corUn = covUn/(sdUn*sdUn[idx[iLoc]])
            
            # Correct covariance for constraint
            tmpRHS = A%*%covUn
            tmpCov = covUn - tr%*%solve(constr.Sigma, tmpRHS)
            covMat[,iLoc] = covMat[,iLoc] + as.vector(tmpCov)
            
            # Correct correlation
            corMat[,iLoc] = corMat[,iLoc] + as.vector(tmpCov/(tmpSD*tmpSD[idx[iLoc]]))
        }
    }
    covMat = covMat/numThetas
    sdMat  = sdMat/numThetas
    corMat = corMat/numThetas
    
    # Print the stuff in Matlab
    write.table(covMat, 'delCovMat.txt', row.names = F, col.names = F)
    write.table(corMat, 'delCorMat.txt', row.names = F, col.names = F)
    write.table(sdMat, 'delSDMat.txt', row.names = F, col.names = F)
    
### Calculate MAP structure
    ## Find map
    idx.MAP = which.max(mcmc.ns$hyper$log.p)
    theta.MAP = mcmc.ns$hyper$samples[idx.MAP,]
    
    ## Select test locations
    idx = c(365, 
            936,
            120,
            316, 
            385)
    
    ## Initialize storage
    covMat = matrix(0, nrow = mesh$n, ncol = length(idx))
    sdMat  = matrix(0, nrow = mesh$n, ncol = 1)
    corMat = covMat
    
    ## Calculate
        # Get theta
        theta = theta.MAP[-1]
        
        # Get unrestricted precision matrix
        Qxx = INLA:::inla.spde3.precision(spde = mcmc.ns$data$spde, theta = theta)
        
        # Factorize
        L = Cholesky(Qxx)
        
        # Marginal standard devations
        sdUn = sqrt(diag(inla.qinv(Qxx)))
        
        ## Precomputation
        # Covariance matrix of constraints
        A = mcmc.ns$constr[, -c(1236,1237)]
        constr.Sigma = A%*%solve(L, t(A))
        
        # Transfer
        tr = solve(L, t(A))
        
        # Correct sd according to constraints
        tmpSD = sqrt(sdUn^2 - diag(tr%*%solve(constr.Sigma, t(tr))))
        sdMat = sdMat + tmpSD
        
        # Run through each location
        for(iLoc in 1:length(idx)){
            # Unconstrained covariance
            rhs = rep(0, mesh$n)
            rhs[idx[iLoc]] = 1
            covUn = solve(L, rhs)
            
            # Unconstrained correlation
            corUn = covUn/(sdUn*sdUn[idx[iLoc]])
            
            # Correct covariance for constraint
            tmpRHS = A%*%covUn
            tmpCov = covUn - tr%*%solve(constr.Sigma, tmpRHS)
            covMat[,iLoc] = covMat[,iLoc] + as.vector(tmpCov)
            
            # Correct correlation
            corMat[,iLoc] = corMat[,iLoc] + as.vector(tmpCov/(tmpSD*tmpSD[idx[iLoc]]))
        }
    
    # Print the stuff in Matlab
    write.table(covMat, 'delMAPCovMat.txt', row.names = F, col.names = F)
    write.table(corMat, 'delMAPCorMat.txt', row.names = F, col.names = F)
    write.table(sdMat, 'delMAPSDMat.txt', row.names = F, col.names = F)

## Get average range and std.dev.
    avRange = 1:mesh$n
    avRange = 0*avRange
    avStdDev = avRange
    for(i in 1:length(mcmc.ns$hyper$log.p)){
        avRange  = avRange  + exp(mcmc.ns$hyper$samples[i, 2] + covCov%*%mcmc.ns$hyper$samples[i, 5:6])
        avStdDev = avStdDev + exp(mcmc.ns$hyper$samples[i, 3] + covCov%*%mcmc.ns$hyper$samples[i, 8:9])
    }
    avRange = avRange/length(mcmc.ns$hyper$log.p)
    avStdDev = avStdDev/length(mcmc.ns$hyper$log.p)
    
    # Print in Matlab
    write.table(avRange, 'delAvRange.txt', row.names = F, col.names = F)
    write.table(avStdDev, 'delAvStdDev.txt', row.names = F, col.names = F)
    
## Posteriors for non-stationarity
    # Generate realizations of non-stat
    num = 7000
    lRange = matrix(0, nrow = num, ncol = length(idx))
    lStdDev = lRange
    idxT = 1:7000
    locIdx = idx
    for(i in 1:num){
        cIdx = idxT[i]
        lRange[i,]  = exp(covCov[locIdx,]%*%mcmc.ns$hyper$samples[cIdx, 5:6])
        lStdDev[i,] = exp(covCov[locIdx,]%*%mcmc.ns$hyper$samples[cIdx, 8:9])
    }
    
    # Generate plots
    numToPlot = 1
    hist(lRange[,numToPlot], freq = FALSE, 20, xlab = "Non-stationary effect", main = "")
    dev.copy2eps(file = paste(bDir, "nonStatEffectRange.eps", sep = ""))
    
    hist(lStdDev[, numToPlot], freq = FALSE, 20, xlab = "Non-stationary effect", main = "")
    dev.copy2eps(file = paste(bDir, "nonStatEffectStdDev.eps", sep = ""))

    # Get covariate values on mesh
    predCovHeight = getHeight(mcmc.ns$data$mesh$loc[, c(1,2)])
    predCovV = cbind(1, predCovHeight$hCov/1e3)
    
    # Set sea equal to zero
    predCovV[predCovV[,2]<0, 2] = 0;
    predCovV[is.na(predCovV[,2]), 2] = 0;
    
    # Make predictions on triangulation
    meshPred = colMeans(mcmc.ns$pred$samplesLatent)
    predOnMesh = meshPred[1:1235] + predCovV%*%meshPred[1236:1237]
    write.table(predOnMesh, 'delPredOnMesh.txt', row.names = F, col.names = F)
