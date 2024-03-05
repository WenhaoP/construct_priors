################################################################################
# Sec4_sim.R                                                                   #
#    Check coverage for different values of the hyperparameters                #
################################################################################

## Library
    # Load libraries
    library(MASS)
    library(INLA)
    source('deformationModel.R')

## Data
    # Load data
    load("data.RData")
    
    # Set mesh
    mesh = mesh
    
    # Re-scale observation	
    prec = data$precipitation/1000
    
    # Observation data
    obs = prec
    loc = cbind(data$x, data$y)
    
    # Get covariates from digital elevation model (GLOBE)
    linCov = cbind(1, data$altitude/1e3)
    covCov = NULL

### Part 1: Run stationary model on data
    ## Set up problem
    # Make model object
    statObj = inla.spde3.deform.makeModel(obs = obs,
                                          loc = loc,
                                          mesh = mesh,
                                          linCov = linCov,
                                          covRange = covCov,
                                          covStd = covCov)
    # Add prior settings
    statObj = inla.spde3.deform.setPrior(statObj,
                                         uN = 3, aN = 0.05,
                                         lR = 10, aR = 0.05,
                                         uS = 3, aS = 0.05)	
    
    # Add constraints to spatial field
    constr = as.matrix(rep(1, length(h)))
    statObj = inla.spde3.deform.setTarget(modelObj = statObj, constr = constr)

## Fit model
    # INLA
    guess = c(4.2, 5.5, 0)
    result = inla.spde3.deform.analyze(statObj, control.mode = list(theta = guess, restart = TRUE))

### Part 2: Simulate stationary data according to estimated model
###         and fit non-stationary model with *only* spatial field and noise?
    ## Covariates in first-order structure
    linCov = NULL
    
    ## Covariates in second-order structure
    source('readHeight.R')
    tmpCov = getCovariates(mesh)
    covHeight = tmpCov$covHeight
    covHeight = covHeight/1e3
    covGrad = tmpCov$covGrad
    covGrad = covGrad/1e2
    covCov = cbind(as.vector(covHeight), as.vector(covGrad))
    covCov = inla.spde3.deform.centerCovariates(Cov = covCov, mesh = mesh)
    
    ## Set up problem
    # Make model object
    nonStatObj = inla.spde3.deform.makeModel(obs = obs,
                                             loc = loc,
                                             mesh = mesh,
                                             linCov = linCov,
                                             covRange = covCov,
                                             covStd = covCov)
    # Add prior settings
    lamR = 20
    lamS = 20
    nonStatObj = inla.spde3.deform.setPrior(nonStatObj,
                                            uN = 3, aN = 0.05,
                                            lR = 10, aR = 0.05,
                                            uS = 3, aS = 0.05,
                                            lamRange = lamR, lamStdDev = lamS)	
    mcmc.priorProb(lamR, 100000, nonStatObj) # Check lambda value
    
    # Add constraints to spatial field
    constr = NULL
    nonStatObj = inla.spde3.deform.setTarget(modelObj = nonStatObj, constr = constr)
    
    ## Run through simulations
    guess = result$summary.hyperpar$mode
    guess[1] = log(guess[1])
    guess = c(guess, rep(1e-3, 6))
    
    # Set seed for reproducability
    require(MASS)
    #      set.seed(3448534)
    
    # Load required packages
    require(doParallel)
    cl = makeCluster(8)
    registerDoParallel(cl)
    require(foreach)
    numReal = 200
    load('SimStudyStartValue.RData')
    allRes = foreach(iter = 1:numReal, .combine = 'c', .packages = c("INLA")) %dopar% {
        print(iter)
        
        # Generate observation
        simObj = statObj;
        simObj$constr = NULL;
        simObj$data$linCov = NULL;
        obs = inla.spde3.deform.simulateStatData(n = 1, inlaRes = result, modelObj = simObj)
        
        # Change observation
        nonStatObj$data$obs = obs
        
        # Re-fit
        # NOTE: THE NON-STATIONARY MODEL REQUIRES A SPECIAL VERSION OF INLA
        res = try(inla.spde3.deform.analyze(nonStatObj, num.threads = 1, control.mode = list(theta = init, restart = TRUE)))
        
        if(class(res) == "try-error"){
            res = NULL
        }
        
        # Return
        res$obsUsed = obs
        list(res)
    }
    
    trueVal = result$summary.hyperpar$mode
    trueVal[1] = log(trueVal[1])
    trueVal = c(trueVal, rep(0,6))
    
    # Make a useful name stamp
    t = Sys.time()
    tmpStmp = strftime(t, "%Y-%m-%d-%H-%M-%S")
    fName = paste('Results/simulationStudyData_NumReal_', numReal, '_time_', tmpStmp, '_lamRange_', lamR, '_lamStdDev_', lamS, '.RData', sep = "")
    save(allRes, trueVal, file = 'Results/simulationStudyData.Rdata')

### Part 3: Process results
    ## Check coverage
    # Initialize
    freqCov = vector("double", length = 9)
    thetas = c()
    trueVal = result$summary.hyperpar$mode
    trueVal[1] = log(trueVal[1])
    trueVal = c(trueVal, rep(0,6))
    count = 0
    
    # Join everything together
    for(iter in 1:length(allRes)){
        print(iter)
        if(is.null(allRes[[iter]]$summary.hyperpar))
            next
        
        # Need sampling to get non-linear combinations
        theta = inla.hyperpar.sample(n = 100000, allRes[[iter]], intern = TRUE)
        theta[,5:6] = theta[,5:6]*exp(-0.5*theta[,4])
        theta[,8:9] = theta[,8:9]*exp(-0.5*theta[,7])
        
        # Find quantiles
        qMat = apply(theta, 2, quantile, probs = c(0.025, 0.975))
        for(pNum in 1:9){
            if((qMat[1,pNum] < trueVal[pNum]) && (trueVal[pNum] < qMat[2,pNum])){
                freqCov[pNum] = freqCov[pNum] + 1
            }
        }
        
        # Aggregate thetas to "integrate out" data-generating mechanism
        thetas = rbind(thetas, theta)
        
        count = count+1
    }
    
    # Frequentist coverages
    print(freqCov/count)

    
