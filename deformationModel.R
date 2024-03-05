################################################################################
# deformationModel.R                                                           #
#    Code for the non-stationary SPDE model based on a spatially               #
#    varying metric tensor.                                                    #
#                                                                              #
#    Author: Geir-Arne Fuglstad <geirarne.fuglstad@gmail.com>                  #
#    Affiliation: Norwegian University of Science and Technology               #
################################################################################

# Check that INLA library is loaded
require(INLA)

# Gaussian prior
inla.spde3.deform.logGauss = function(x, S){
    # Calculate log density
    #require(mvtnorm)
    S = as.matrix(S)
    nPar = dim(S)[2]
    logDens = -0.5*nPar*log(2*pi) + 0.5*determinant(S, logarithm = TRUE)$modulus - 0.5*t(x)%*%S%*%x
    #logDens = dmvnorm(x, mean = rep(0, nPar), sigma = solve(S), log = TRUE)
    
    return(as.numeric(logDens))
}

# PC prior for range and standard deviation
inla.spde3.deform.logStat = function(x, par){
    # Range
    lambdaRange = par[1]
    logDens = log(lambdaRange) - 2*x[1] - lambdaRange*exp(-x[1]) + x[1]
    
    # Standard deviation
    lambdaStd = par[2]
    logDens = logDens + log(lambdaStd) - lambdaStd*exp(x[2]) + x[2]
    
    return(logDens)
}

## TMP: This wrapper function will be called by the embedded R engine
log.PC.inla = function(x){
    return(inla.spde3.deform.logPC(x, MCMC = FALSE))
}

## TMP: New wrapper not loading stuff from disk
log.PC.inla.noRead = function(x){
    return(inla.spde3.deform.logPC(x, MCMC = FALSE, READ = FALSE))
}

# Full prior for model
inla.spde3.deform.logPC = function(x, modelObj = NULL, MCMC = TRUE, READ = TRUE, hyper = TRUE){
    # Load libraries
    require(INLA)
    
    # Load Gramians
    if(!MCMC){
        if(READ){
            load(file = "tmpNonameSmatrices.RData")
            load(file = "tmpNonameHyperparameters.RData")
        }
        
        # Calculate stuff
        d = 2
        lambdaRange =-lR^(d/2)*log(aR)
        
        # Parameter for log(std.dev.)
        lambdaStdDev = -log(aS)/uS
        
        # Penalisation rates for stationary model
        parStat = c(lambdaRange, lambdaStdDev)
        
    } else{
        ## Prior on noise
        uN = modelObj$prior$data$noise$U
        aN = modelObj$prior$data$noise$alpha
        
        ## Prior on stationary SPDE (log(range) and log(Std.Dev.)
        # Parameter for log(range)
        lR = modelObj$prior$data$range$L
        aR = modelObj$prior$data$range$alpha
        d = 2
        lambdaRange =-lR^(d/2)*log(aR)
        
        # Parameter for log(std.dev.)
        uS = modelObj$prior$data$stdDev$U
        aS = modelObj$prior$data$stdDev$alpha
        lambdaStdDev = -log(aS)/uS
        
        # Penalisation rates for stationary model
        parStat = c(lambdaRange, lambdaStdDev)
        
        ## Non-stat prior
        Sr = modelObj$prior$data$nonStat$Sr
        lamRange = modelObj$prior$data$nonStat$lamRange
        Ss = modelObj$prior$data$nonStat$Ss
        lamStdDev = modelObj$prior$data$nonStat$lamStdDev
    }
    
    # Prior on Gaussian noise
    log.prior = inla.pc.dprec(prec = exp(x[1]), u = uN, alpha = aN, log = TRUE)+x[1]
    x = x[-1]
    
    # Stationary prior
    xStat = x[1:2]
    x = x[-(1:2)]
    log.prior = log.prior+inla.spde3.deform.logStat(xStat, parStat)
    
    # Range
    if(!is.null(Sr)){
        # Calculate prior on log-scaling
        if(hyper){
            lTau = x[1]
            x = x[-1]
            lam = lamRange
            log.prior = log.prior + inla.pc.dprec(prec = exp(lTau), lam = lam, log = TRUE) + lTau
        } else{
            Sr = Sr*lamRange
        }
        
        # Calculate prior on range parameters
        nPar = dim(Sr)[1]
        xRange = x[1:nPar]
        x = x[-(1:nPar)]
        log.prior = log.prior + inla.spde3.deform.logGauss(xRange, Sr)
    }
    
    # Standard deviation
    if(!is.null(Ss)){
        # Calculate prior on log-scaling
        if(hyper){
            lTau = x[1]
            x = x[-1]
            lam = lamStdDev
            log.prior = log.prior + inla.pc.dprec(prec = exp(lTau), lam = lam, log = TRUE) + lTau
        } else {
            Ss = Ss*lamStdDev
        }
        
        # Calculate prior on standard deviation
        nPar = dim(Ss)[1]
        xStd = x[1:nPar]
        x = x[-(1:nPar)]
        log.prior = log.prior + inla.spde3.deform.logGauss(xStd, Ss)
    }
    
    return(log.prior)
}

# Create spde object of new model using the experimental SPDE3 interface
inla.spde3.deform.makeSPDE = function(mesh, covRange, covStd){
    # Fix input (NOTICE MINUS SIGN)
    if(!is.null(covRange)){
        covRange = -as.matrix(covRange)
    }
    if(!is.null(covStd)){
        covStd   = -as.matrix(covStd)
    }
    
    # Make sure INLA library is loaded
    require(INLA)
    
    # Generate FEM matrices
    femMatrix = inla.mesh.fem(mesh)
    C = femMatrix$c0
    G = femMatrix$g1
    
    # Make M matrices
    M0 = C
    M1 = 0*C
    M2 = G+t(G)
    M3 = diag(1/sqrt(diag(C)))%*%G
    
    # Precalculate parameters and fix smoothness
    alpha = 2
    d = 2
    nu = alpha - d/2
    kappa0 = log(sqrt(8*nu))
    tau0   = -0.5*log(4*pi)
    
    # Need to add extra column that has no effect to use
    # as hyperparameter in prior
    if(is.null(covRange)){
        rangeTau = NULL
    } else{
        rangeTau = 0
    }
    
    if(is.null(covStd)){
        stdTau = NULL
    } else {
        stdTau = 0
    }
    
    # Assemble B matrices
    B0 = cbind(   tau0,  0, -1, rangeTau, 0*covRange, stdTau,   covStd)
    B1 = cbind( kappa0, -1,  0, rangeTau,   covRange, stdTau, 0*covStd)
    B2 = cbind(      0,  0,  0, rangeTau, 0*covRange, stdTau, 0*covStd)
    B3 = -2*B1
    
    ## Set temporary prior
        #    This prior will NOT be used
        numPar = 0
        if(!is.null(covRange)){
            numPar = numPar + dim(covRange)[2]+1
        }
        if(!is.null(covStd)){
            numPar = numPar + dim(covStd)[2]+1
        }
        
        theta.mu = c(3, 0, rep(0,numPar))
        Q = diag(1, nrow = length(theta.mu), ncol = length(theta.mu))
    
    # Make spde object
    res = INLA:::inla.spde3.generic(M0 = M0, M1 = M1, M2 = M2, M3 = M3,
                                    B0 = B0, B1 = B1, B2 = B2, B3 = B3,
                                    theta.mu = theta.mu, theta.Q = Q,
                                    transform = "log")
    return(res)
    
}

# Centre covariates
inla.spde3.deform.centerCovariates = function(Cov,
                                              mesh){
    # Extract integration matrix
    C0 = inla.mesh.fem(mesh)$c0
    
    # Make it integrate to zero
    Cov = Cov - matrix(rep(colSums(C0%*%Cov)/sum(diag(C0)), each = dim(Cov)[1]), ncol = dim(Cov)[2])
}

# Collect information in an object
inla.spde3.deform.makeModel = function(obs,
                                       loc,
                                       mesh,
                                       linCov,
                                       covRange,
                                       covStd){
    # Construct A matrix to map observations to mesh
    A = inla.spde.make.A(mesh = mesh, loc = loc)
    
    # Construct SPDE object
    spde = inla.spde3.deform.makeSPDE(mesh, covRange, covStd)
    
    ## Avoid annoying situations
        # Fix input
        if(!is.null(covRange))
            covRange = as.matrix(covRange)
        if(!is.null(covStd))
            covStd   = as.matrix(covStd)
    
    # Store data
    data = list(obs = obs,
                loc = loc,
                linCov = linCov,
                covRange = covRange,
                covStd = covStd,
                mesh = mesh,
                A = A,
                spde = spde)
    
    # Store everything
    modelObj = list(data = data)
    
    return(modelObj)
}

# Configure priors
inla.spde3.deform.setPrior = function(modelObj,
                                      tauB = 1e-4,
                                      uN = NULL, aN = NULL, 
                                      lR = NULL, aR = NULL, 
                                      uS = NULL, aS = NULL, 
                                      lamRange = NULL, lamStdDev = NULL){
    # TODO: How to set these values in practice?
    
    # Fixed effects
    fixedPar = list(tauB = tauB)
    
    # Gaussian measurement noise
    if(is.null(uN)){
        uN = 3
        aN = 0.05
    }
    noisePar = list(U = uN, alpha = aN)
    
    # Spatial range
    if(is.null(lR)){
        lR = 10
        aR = 0.05
    }
    rangePar = list(L = lR, alpha = aR)
    
    # Process variance
    if(is.null(uS)){
        uS = 3
        aS = 0.05
    }
    stdDevPar = list(U = uS, alpha = aS)
    
    ## Non-stationarity
        # Calculate the Gramian for Range
        C0 = inla.mesh.fem(mesh)$c0
        if(!is.null(modelObj$data$covRange)){
            covRange = modelObj$data$covRange
            covInt = C0%*%covRange
            Sr = t(covInt)%*%covRange
            Sr = Sr/sum(diag(C0))
        } else {
            Sr = NULL
        }
        
        # Calculate the Gramian for Std.Dev.
        if(!is.null(modelObj$data$covStd)){
            covStd = modelObj$data$covStd
            covInt = C0%*%covStd
            Ss = t(covInt)%*%covStd
            Ss = Ss/sum(diag(C0))
        } else {
            Ss = NULL
        }
    
        # Set range parameter
        if(is.null(lamRange)){
            lamRange = 1
        }
        
        # Set standard deviation parameter
        if(is.null(lamStdDev)){
            lamStdDev = 1
        }
        
        # Collect
        nonStatPar = list(lamRange = lamRange, lamStdDev = lamStdDev, Sr = Sr, Ss = Ss)
    
    ## This is only used by MCMC sampler
        # Set prior function
        fun = inla.spde3.deform.logPC
    
    # Collect everything
    data = list(fixed = fixedPar,
                noise = noisePar,
                range = rangePar,
                stdDev = stdDevPar,
                nonStat = nonStatPar)
    
    # Add to model object
    modelObj$prior = list(data = data,
                          fun = fun)
    
    return(modelObj)
    
}

# Function to collect priors in a prior object
inla.spde3.deform.setTarget = function(modelObj,
                                       constr = NULL){
    
    if(is.null(constr)){
        modelObj$constr = NULL
        modelObj$target = nonStat.logTarget
    } else {
        # Extract integration matrix
        C0 = inla.mesh.fem(modelObj$data$mesh)$c0
        
        # Set constraints
        modelObj$constr = C0%*%constr/sum(diag(C0))
        
        # Constraint is on spatial field
        pCol = dim(constr)[2]
        pRow = dim(modelObj$data$linCov)[2]
        modelObj$constr = rBind(modelObj$constr, Matrix(0, ncol = pCol, nrow = pRow))
        
        # Correct rotation
        modelObj$constr = t(modelObj$constr)
        
        # Target
        modelObj$target = nonStat.logTarget
    }
    
    return(modelObj)
}

inla.spde3.deform.analyze = function(modelObj, compute = FALSE, ...){
    # Extract spde stuff
    spde.ns = modelObj$data$spde
    A.sites = modelObj$data$A
    mesh = modelObj$data$mesh
    
    # Formula to describe model
    linCov = modelObj$data$linCov
    if(!is.null(modelObj$constr)){
        Aconst = as.matrix(modelObj$constr)
        Aconst = Aconst[,1:(dim(Aconst)[2]-dim(linCov)[2]), drop = FALSE]
        econst = matrix(0, ncol = 1, nrow = dim(Aconst)[1])
        extraconstr = list(A = Aconst, e = econst)
    } else{
        extraconstr = NULL
    }
    
    if(!is.null(linCov)){
        formula.ns = y ~ -1 + linCov + f(field, model = spde.ns, extraconstr = extraconstr)
    } else{
        formula.ns = y ~ -1 + f(field, model = spde.ns, extraconstr = extraconstr)
    }
    
    ## Collect data in a stack
    # Estimation stack
    obs = as.vector(modelObj$data$obs)
    if(!is.null(linCov)){
        data.stack.est = inla.stack(tag = "est",
                                    data = list(y = obs),
                                    A = list(A.sites,1),
                                    effects = list(list(field = 1:mesh$n),
                                                   list(linCov = linCov)))
    } else {
        data.stack.est = inla.stack(tag = "est",
                                    data = list(y = obs),
                                    A = list(A.sites),
                                    effects = list(list(field = 1:mesh$n)))
    } 
    
    # Combine stacks:
    data.stack = inla.stack(data.stack.est)
    
    ## Store stuff on disk for INLA prior calculations
        # S matrices
        Sr = modelObj$prior$data$nonStat$Sr
        Ss = modelObj$prior$data$nonStat$Ss

        # Hyperparameters stationary part of the model
        uN = modelObj$prior$data$noise$U
        aN = modelObj$prior$data$noise$alpha
        lR = modelObj$prior$data$range$L
        aR = modelObj$prior$data$range$alpha
        uS = modelObj$prior$data$stdDev$U
        aS = modelObj$prior$data$stdDev$alpha
        
        # Hyperparameters for non-stationary part
        lamRange = modelObj$prior$data$nonStat$lamRange
        lamStdDev = modelObj$prior$data$nonStat$lamStdDev
        
        # Store variables needed for prior
        fname = tempfile()
        save(Sr, Ss, uN, aN, lR, aR, uS, aS, lamRange, lamStdDev, file = fname)
        
    # Need specially compiled version of INLA that
    # can run priors from separate R files
    # INLA:::inla.my.update(binaries=TRUE)
    
    # Prior for fixed effects
    tauB = modelObj$prior$data$fixed$tauB
    
    # Estimate model
    result.ns = inla(formula.ns,
                     family = "gaussian",
                     data = inla.stack.data(data.stack),
                     control.fixed = list(prec = tauB, prec.intercept = tauB),
                     control.predictor = list(A = inla.stack.A(data.stack), compute = compute),
                     control.expert = list(jp.Rfile = "deformationModel.R",
                                           jp.func  = "log.PC.inla.noRead",
                                           jp.RData = fname),
                     ...)
    
    # Delete file
    unlink(fname)
    
    
    return(result.ns)
}

inla.spde3.deform.simulateStatData = function(n = 1, inlaRes, modelObj){
    # Extract estimated parameters
    theta = inlaRes$mode$theta
    
    obs = c()
    for(i in 1:n){
        # Simulate latent field
        numSpatial = dim(modelObj$data$A)[2]
        xSample = mcmc.sampleLatent(n = 1, x = theta, modelObj = modelObj, cond = FALSE)
        
        # Combine into observations
        xObs = modelObj$data$A%*%xSample[1:numSpatial]
        
        # Add measurement noise
        xObs = xObs + rnorm(length(xObs))*exp(-0.5*theta[1])
        
        # Combine
        obs = cbind(obs, as.vector(xObs))
    }
    
    return(obs)
}

inla.spde3.deform.plotPosterior = function(inlaRes, num, tFun = NULL, add = FALSE, ...){
    # No transform specified
    if(is.null(tFun)){
        tFun = function(x){x}
    }
    
    # Plot settings
    bty = "l"
    lwd = 2L
    cex.lab = 1.4
    cex.axis = 1.4
    save.figure = FALSE
    
    # Plot
    par(cex.lab = cex.lab, cex.axis = cex.axis)
    marg = inla.smarginal(inla.tmarginal(tFun, inlaRes$marginals.hyperpar[[num]]))
    if(add){
        lines(marg)
    } else{
        plot(marg, type = 'l', ...)
    } 
}

## Functions used by MCMC sampler
    # Sample hyperparameters of the model
    mcmc.addHyperSamples = function(modelObj,
                                    numSamples,
                                    inits,
                                    scale,
                                    acc.rate,
                                    adapt = TRUE,
                                    hyper = TRUE){
        # Load sampler library
        require(adaptMCMC)
        
        # Run
        hmm = MCMC(p = modelObj$target,
                   n = numSamples, 
                   init = inits, 
                   scale = scale, 
                   adapt = adapt,
                   acc.rate = acc.rate,
                   modelObj = modelObj,
                   hyper = hyper)
        
        # Transform to desired parametrization
        if(hyper){
            # Parameters in log(range)
            idxEnd = 3
            numRange = dim(modelObj$data$covRange)[2]
            if(is.null(numRange)){
                numRange = 0
            } else{
                idxStart = idxEnd+2
                idxEnd   = idxStart+numRange-1
                hmm$samples[,idxStart:idxEnd] = hmm$samples[,idxStart:idxEnd]*exp(-0.5*hmm$samples[,idxStart-1])
            }
            
            numStdDev = dim(modelObj$data$covStd)[2]
            if(is.null(numStdDev)){
                numStdDev = 0
            } else{
                idxStart = idxEnd+2
                idxEnd   = idxStart+numStdDev-1
                hmm$samples[,idxStart:idxEnd] = hmm$samples[,idxStart:idxEnd]*exp(-0.5*hmm$samples[,idxStart-1])
            }
        } else{
            # Add fixed taus
            numRange = dim(modelObj$data$covRange)[2]
            idxEnd = 3
            tLen = dim(hmm$samples)[2]
            if(is.null(numRange)){
                numRange = 0
            } else{
                hmm$samples = cbind(hmm$samples[,1:idxEnd, drop = FALSE], 0, hmm$samples[,(idxEnd+1):tLen, drop = FALSE])
                idxEnd = idxEnd+1+numRange
            }
            
            # Insert extra parameter
            numStdDev = dim(modelObj$data$covStd)[2]
            tLen = dim(hmm$samples)[2]
            if(is.null(numStdDev)){
                numStdDev = 0
            } else{
                hmm$samples = cbind(hmm$samples[,1:idxEnd, drop = FALSE], 0, hmm$samples[,(idxEnd+1):tLen, drop = FALSE])
            }
            
        }
        
        # Add result to object
        modelObj$hyper = hmm
        
        return(modelObj)
    }

    # Sample latent field
    mcmc.sampleLatent = function(n = 1, x, modelObj, cond = TRUE){
        # Prior on fixed effects
        tauB = modelObj$prior$data$fixed$tauB
        
        # Calculate model distribution
        theta = x[-1]
        Qxx = INLA:::inla.spde3.precision(spde = modelObj$data$spde, theta = theta)
        p = dim(modelObj$data$linCov)[2]
        if(is.null(p))
            p = 0
        m = dim(Qxx)[1]
        Qbb = diag(rep(tauB, p), ncol = p, nrow = p)
        Qxb = Matrix(0, ncol = p, nrow = m)
        Q = rBind(cBind(Qxx, Qxb), cBind(t(Qxb), Qbb))
        
        # Calculate A matrix
        A = cBind(modelObj$data$A, modelObj$data$linCov)
        
        
        obs = modelObj$data$obs
        numLatent = dim(Q)[1]
        
        if(cond){
            # Calculate conditional distribution
            tau0 = exp(x[1])
            Qc = Q + tau0*t(A)%*%A
            
            # Precompute cholesky factors
            Lc = Cholesky(Qc, LDL = FALSE)
            
            # Find conditional mean
            muC = solve(Lc, tau0*t(A)%*%obs)
            
            # Permute
            muCperm = muC[Lc@perm+1,,drop=FALSE]
            
            # Sample
            meshValues = as.vector(solve(Lc, rnorm(numLatent), system = 'Lt') + muCperm)
        } else {
            Lc = Cholesky(Q, LDL = FALSE)
            meshValues = as.vector(solve(Lc, rnorm(numLatent), system = 'Lt'))
        }
            
        # Permute back
        meshValues = meshValues[invPerm(as.integer(Lc@perm+1))]
        
        if(!is.null(modelObj$constr)){
            # Sample unconstrained
            xSample = meshValues
            
            # Correct
            Acon = modelObj$constr
            Sigma1 = Acon%*%solve(Lc, t(Acon))
            tmpX = t(Acon)%*%solve(Sigma1, Acon%*%xSample)
            xStar = xSample - solve(Lc, tmpX)
            
            # Store
            meshValues = as.vector(xStar)
        }
        
        return(meshValues)
    }
    
    
    # Get prediction samples
    mcmc.addPredSamples = function(modelObj, loc, covValues, thinning){
        # Convert locations to mesh combinations
        Apred = inla.spde.make.A(mesh = modelObj$data$mesh, loc = loc)
        Apred = cbind(Apred, covValues)
        
        # Extract thetas
        thetaAll = modelObj$hyper$samples
        
        # Thin
        thetaAll = thetaAll[seq(1, dim(thetaAll)[1], thinning), ]
        numSamples = dim(thetaAll)[1]
        
        # Store
        numPred = dim(loc)[1]
        modelObj$pred = list(samples = matrix(0, nrow = numSamples, ncol = numPred))
        modelObj$pred$samplesNN = modelObj$pred$samples
        p = dim(modelObj$data$linCov)[2]
        if(is.null(p))
            p = 0
        numLat = modelObj$data$mesh$n + p
        modelObj$pred$samplesLatent = matrix(0, nrow = numSamples, ncol = numLat)
        modelObj$pred$hSamples = thetaAll
        
        # Sample
        idxBreaks = floor(seq(1, dim(thetaAll)[1], length.out = 21))
        idxBreaks[1] = 0
        for(per in 1:20){
            print(paste((per-1)*5, "% done", sep = ""))
            if(idxBreaks[per] == idxBreaks[per+1]){
                next
            }
            for(idx in (idxBreaks[per]+1):idxBreaks[per+1]){
                # Get current sample of hyperparameters
                x = thetaAll[idx, ]
                
                # Get sample of latent field conditional on hyperpars
                meshValues = mcmc.sampleLatent(n = 1, x, modelObj)
                modelObj$pred$samplesLatent[idx,] = meshValues
                
                # Convert to prediction locations
                meshValues = as.vector(Apred%*%meshValues)
                
                # Calculate values at prediction locations without noise
                modelObj$pred$samplesNN[idx,] = meshValues
                
                # Predict with noise
                tau0 = exp(x[1])
                modelObj$pred$samples[idx,] = modelObj$pred$samplesNN[idx,] + rnorm(numPred)/sqrt(tau0)
                
            }
        }
        
        return(modelObj)
    }
    
    # Post-calculate leave-one-out predictions
    mcmc.LOOP = function(modelObj, WAIC = FALSE, CRPS = FALSE){
        # Extract stuff
        obs = modelObj$data$obs
        hSamples = modelObj$pred$hSamples
        pSamples = modelObj$pred$samplesNN
        
        # Get dimensions
        modelObj$pred$LOOLpost = modelObj$pred$samples*0
        numSamples = dim(modelObj$pred$LOOLpost)[1]
        numPred = dim(modelObj$pred$LOOLpost)[2]
        
        # Temporary store CRPS values
        tmpCRPS = modelObj$pred$LOOLpost
        
        # Temporarp store point predictions
        tmpPP = modelObj$pred$LOOLpost
        
        # Temporary store E[y_i^2| y_{-i}]
        tmpE2 = tmpPP
        
        # Calculate all densities \pi(y_i| x^s, \theta^s), where
        # (x^s, \theta^s) are drawn from \pi(x, \theta | y)
        for(oIdx in 1:numPred){
            for(idx in 1:numSamples){
                modelObj$pred$LOOLpost[idx, oIdx] = dnorm(obs[oIdx], mean = pSamples[idx, oIdx], sd = exp(-0.5*hSamples[idx,1]))
                
                # Get mean and variance of predictive dist given x and \theta
                mu = pSamples[idx, oIdx]
                sig = exp(-0.5*hSamples[idx,1])
                      
                # Calculate point predictions given x and \theta
                tmpPP[idx, oIdx] = mu
                
                # Calculate expected value of square
                tmpE2[idx, oIdx] = sig^2+mu^2
            }
        }
        
        ## WAIC
            # Estimate efficient number of parameters
            varWAIC = apply(log(modelObj$pred$LOOLpost), 2, var)
            pWAIC = sum(varWAIC)
            
            # Log-predictive density
            lpdWAIC = sum(log(colMeans(modelObj$pred$LOOLpost)))
            
            # In deviance scale
            WAIC = -2*(lpdWAIC-pWAIC)
            
            # Store it
            modelObj$pred$pWAIC = pWAIC
            modelObj$pred$lpdWAIC = lpdWAIC
            modelObj$pred$WAIC = WAIC
        
        ## CRPS (Not a very nice way to do it...)
            if(CRPS){
                CRPS = matrix(0, ncol = 1, nrow = numPred)
                
                for(cObs in 1:numPred){
                    print(cObs)
                    # Calculate with importance sampling E|Y-y|
                    allMu = pSamples[, cObs] - obs[cObs]
                    allSD = exp(-0.5*hSamples[, 1])
                    W = 1/modelObj$pred$LOOLpost[, cObs]
                    E1 = sum(W*abs(rnorm(numSamples, mean = allMu, sd = allSD)))/sum(W)
                    
                    # Calculate with importance sampling E|Y-Y'|
                    allMu = rep(pSamples[, cObs], each = numSamples) - rep(pSamples[, cObs], numSamples)
                    allSD = rep(exp(-0.5*hSamples[, 1]), each = numSamples) + rep(exp(-0.5*hSamples[,1]), numSamples)
                    W = rep(1/modelObj$pred$LOOLpost[, cObs], each = numSamples)*rep(1/modelObj$pred$LOOLpost[, cObs], numSamples)
                    E2 = sum(W*abs(rnorm(numSamples^2, mean = allMu, sd = allSD)))/sum(W)
                    
                    CRPS[cObs] = E1-0.5*E2
                    print(CRPS[cObs])
                }
                    
                modelObj$pred$CRPS = CRPS
            }
        
        ## Point predictions
            # Importance sampling
            modelObj$pred$PP = colMeans(tmpPP/modelObj$pred$LOOLpost)/colMeans(1/modelObj$pred$LOOLpost)
            
        ## Prediction std.devs.
            # Importance sampling
            E2 = colMeans(tmpE2/modelObj$pred$LOOLpost)/colMeans(1/modelObj$pred$LOOLpost)
            modelObj$pred$Psig = sqrt(E2-modelObj$pred$PP^2)
            
        ## LOO score		
            modelObj$pred$LOOLpost = log(1/colMeans(1/modelObj$pred$LOOLpost))
        
        res = modelObj
    }
    
    # Evaluate posterior using INLA
    logTargetINLA = function(x, modelObj){
        # Prior on fixed effects
        tauB = modelObj$prior$data$fixed$tauB
        
        # Calculate model distribution
        theta = x[-1]
        Qxx = INLA:::inla.spde3.precision(spde = modelObj$spde, theta = theta)
        p = dim(modelObj$linCov)[2]
        if(is.null(p))
            p = 0
        m = dim(Qxx)[1]
        Qbb = diag(rep(tauB, p), ncol = p, nrow = p)
        Qxb = Matrix(0, ncol = p, nrow = m)
        Q = rBind(cBind(Qxx, Qxb), cBind(t(Qxb), Qbb))
        
        # Calculate A matrix
        A = cBind(modelObj$A, modelObj$A%*%modelObj$linCov)
        
        # Calculate conditional distribution
        tau0 = exp(x[1])
        Qc = Q + tau0*t(A)%*%A
        
        # Calculate linear part
        b = tau0*t(A)%*%obs
        
        
        # Check for constraints
        if(!is.null(modelObj$constr)){
            constr = list(A = modelObj$constr, e = rep(0, dim(modelObj$constr)[1]))
        } else {
            constr = NULL
        }
        
        # Evaluate in x = 0
        xEval = as.matrix(rep(0, dim(Q)[1]))
        
        # pi(x|theta)
        lDens1 = inla.qsample(sample = xEval, Q = Q,  b = xEval, constr = constr)
        
        # pi(y| x, theta)
        lDens2 = 0.5*length(modelObj$obs)*(log(tau0)-log(2*pi)) - 0.5*t(modelObj$obs)%*%modelObj$obs*tau0
        
        # pi(x|theta, y)
        lDens3 = inla.qsample(sample = xEval, Q = Qc, b = b, constr = constr)
        
        # Prior
        lPrior = modelObj$prior$priorFun(x, pars = modelObj$prior)
        
        # Total density
        lDens = lDens1$logdens+lDens2-lDens3$logdens+lPrior
        
        return(as.numeric(lDens))
    }
    
    # Evaluate posterior
    nonStat.logTarget = function(x, modelObj, hyper = TRUE){
        # Prior on fixed effects
        tauB = modelObj$prior$data$fixed$tauB
        
        # Calculate model distribution
        theta = x[-1]
        
        # Reparametrize
        if(hyper){
            numRange = dim(modelObj$data$covRange)[2]
            idxEnd = 2
            if(is.null(numRange)){
                numRange = 0
            } else{
                idxStart = idxEnd+2
                idxEnd   = idxStart+numRange-1
                theta[idxStart:idxEnd] = theta[idxStart:idxEnd]*exp(-0.5*theta[idxStart-1])
            }
    
            numStdDev = dim(modelObj$data$covStd)[2]
            if(is.null(numStdDev)){
                numStdDev = 0
            } else{
                idxStart = idxEnd+2
                idxEnd   = idxStart+numStdDev-1
                theta[idxStart:idxEnd] = theta[idxStart:idxEnd]*exp(-0.5*theta[idxStart-1])
            }
        } else {
            # Insert extra parameter
            numRange = dim(modelObj$data$covRange)[2]
            idxEnd = 2
            if(is.null(numRange)){
                numRange = 0
            } else{
                theta = c(theta[1:idxEnd], 0, theta[(idxEnd+1):length(theta)])
                idxEnd = idxEnd+1+numRange
            }
            
            # Insert extra parameter
            numStdDev = dim(modelObj$data$covStd)[2]
            if(is.null(numStdDev)){
                numStdDev = 0
            } else{
                theta = c(theta[1:idxEnd], 0, theta[(idxEnd+1):length(theta)])
            }
        }

        # This function requires the extra tau parameters, but they are not actually used
        Qxx = INLA:::inla.spde3.precision(spde = modelObj$data$spde, theta = theta)
        p = dim(modelObj$data$linCov)[2]
        if(is.null(p))
            p = 0
        m = dim(Qxx)[1]
        Qbb = diag(rep(tauB, p), ncol = p, nrow = p)
        Qxb = Matrix(0, ncol = p, nrow = m)
        Q = rBind(cBind(Qxx, Qxb), cBind(t(Qxb), Qbb))
        
        # Calculate A matrix
        A = modelObj$data$A
        A = cBind(A, modelObj$data$linCov)
        
        # Calculate conditional distribution
        tau0 = exp(x[1])
        Qc = Q + tau0*t(A)%*%A
        
        # Precompute cholesky factors
        L = Cholesky(Q)
        Lc = Cholesky(Qc)
        
        # Find conditional mean
        obs = modelObj$data$obs
        muC = solve(Lc, tau0*t(A)%*%obs)
        
        # Likelihood: Part 1
        lDetQ = 2*determinant(L, logarithm = TRUE)$modulus
        lDetQc = 2*determinant(Lc, logarithm = TRUE)$modulus
        lDetQn = length(obs)*log(tau0)
        lDet = as.numeric(0.5*lDetQ+0.5*lDetQn-0.5*lDetQc)
        
        # Likelihood: Part 2
        lQuad = 0.5*t(muC)%*%Qc%*%muC -0.5*t(obs)%*%obs*tau0
        lLike = -length(obs)/2*log(2*pi) + lDet + lQuad
        
        # Check for constraints
        if(!is.null(modelObj$constr)){
            Acon = modelObj$constr
            
            # Correction for \pi(x|\theta, Ax = 0)
            muCor = Acon%*%rep(0, dim(Q)[1])
            SigmaCor = Acon%*%solve(L, t(Acon))
            lLike = lLike - (-0.5*determinant(SigmaCor, logarithm = TRUE)$modulus - 0.5*t(muCor)%*%solve(SigmaCor, muCor))
            
            # Correction for \pi(x|\theta, y, Ax = 0)
            muCor = Acon%*%muC
            SigmaCor = Acon%*%solve(Lc, t(Acon))
            lLike = lLike - 0.5*determinant(SigmaCor, logarithm = TRUE)$modulus - 0.5*t(muCor)%*%solve(SigmaCor, muCor)
        } 
        
        # Prior
        lPrior = modelObj$prior$fun(x, modelObj = modelObj, hyper = hyper)
        
        # Total density
        lDens = lLike+lPrior
        
        return(as.numeric(lDens))
    }
    
    mcmc.priorProb = function(lambda, N, modelObj){
        # Simulate hyperhyperparameter
        tau = inla.pc.rprec(n = N, lambda = lambda)
        
        # Simulate hyperparameters
        Sigma = solve(modelObj$prior$data$nonStat$Sr)
        xStar = mvrnorm(n = N, mu = rep(0, dim(Sigma)[1]), Sigma = Sigma)
        
        # Check each realization of covariates
        count = 0
        for(iter in 1:N){
            # Scale according to tau
            x = as.vector(xStar[iter,])/sqrt(tau[iter])
            
            # Calculate values at mesh nodes
            mVals = modelObj$data$covRang%*%x
            
            # Find maximum
            U = max(abs(mVals))
            
            # Check if max relative difference exceeded
            if(U > log(2)){
                count = count + 1
            }
        }
        
        return(count/N)
    }

    mcmc.priorProbGauss = function(lambda, N, modelObj){
        # Simulate hyperhyperparameter
        tau = lambda
        
        # Simulate hyperparameters
        Sigma = solve(modelObj$prior$data$nonStat$Sr)
        xStar = mvrnorm(n = N, mu = rep(0, dim(Sigma)[1]), Sigma = Sigma)
        
        # Check each realization of covariates
        count = 0
        for(iter in 1:N){
            # Scale according to tau
            x = as.vector(xStar[iter,])/sqrt(tau)
            
            # Calculate values at mesh nodes
            mVals = modelObj$data$covRang%*%x
            
            # Find maximum
            U = max(abs(mVals))
            
            # Check if max relative difference exceeded
            if(U > log(2)){
                count = count + 1
            }
        }
        
        return(count/N)
    }
    
    
    
    
## DEBUG FUNCTIONS: REMOVE LATER
    myCheck = function(modelObj, ...){
        # Extract spde stuff
        spde.ns = modelObj$data$spde
        A.sites = modelObj$data$A
        
        # Formula to describe model
        linCov = modelObj$data$linCov
        Aconst = as.matrix(modelObj$constr)
        Aconst = Aconst[,1:(dim(Aconst)[2]-dim(linCov)[2]), drop = FALSE]
        econst = matrix(0, ncol = 1, nrow = dim(Aconst)[1])
        formula.ns = y ~ -1 + linCov + f(field, model = spde.ns, extraconstr = list(A = Aconst, e = econst))
        
        ## Collect data in a stack
        # Estimation stack
        obs = as.vector(modelObj$data$obs)
        data.stack.est = inla.stack(tag = "est",
                                    data = list(y = obs),
                                    A = list(A.sites,1),
                                    effects = list(list(field = 1:mesh$n),
                                                   list(linCov = linCov)))
        
        # Combine stacks:
        data.stack = inla.stack(data.stack.est)
        
        ## Store stuff on disk for INLA prior calculations
        # S matrices
        Sr = modelObj$prior$data$nonStat$Sr
        Ss = modelObj$prior$data$nonStat$Ss
        save(Sr, Ss, file = "tmpNonameSmatrices.RData")
        
        # Hyperparameters stationary part of the model
        uN = modelObj$prior$data$noise$U
        aN = modelObj$prior$data$noise$alpha
        lR = modelObj$prior$data$range$L
        aR = modelObj$prior$data$range$alpha
        uS = modelObj$prior$data$stdDev$U
        aS = modelObj$prior$data$stdDev$alpha
        
        # Hyperparameters for non-stationary part
        lamRange = modelObj$prior$data$nonStat$lamRange
        lamStdDev = modelObj$prior$data$nonStat$lamStdDev
        
        # Store for use in prior calculations
        save(uN, aN, lR, aR, uS, aS, lamRange, lamStdDev, file = "tmpNonameHyperparameters.RData")
        
        # Need specially compiled version of INLA that
        # can run priors from separate R files
        # INLA:::inla.my.update(binaries=TRUE)
        
        # Prior for fixed effects
        tauB = modelObj$prior$data$fixed$tauB
        
        # Estimate model
        result.ns = inla(formula.ns,
                         family = "gaussian",
                         data = inla.stack.data(data.stack),
                         control.fixed = list(prec = tauB, prec.intercept = tauB),
                         control.predictor = list(A = inla.stack.A(data.stack)),
                         control.compute = list(dic = TRUE, cpo = TRUE, waic = TRUE),
                         verbose = TRUE,
                         control.expert = list(jp.Rfile = "deformationModel.R",
                                               jp.func  = "log.PC.inla"),
                         ...)
        
        
        return(result.ns)
    }
