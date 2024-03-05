################################################################################
# Sec3_PriorJe.R                                                               #
#    Simulation study for PriorJe                                              #
################################################################################

# Load stored locations
load('locations.RData')

# Set seed
set.seed(100)

# Load libraries
library(MASS)
library(foreach)
library(doMC)
library(boa)
registerDoMC(10)

# Read functions
source('MCMCsampler.R')
source('helperFunctions.R')

# Model parameters
stdDev = 1
varMod = stdDev^2
nomRange = c(0.1, 1)

# Number of realizations
numReal = 1500

# Initialize storage
jePrior  = array(list(NULL), dim = c(2, numReal))

# Run through settings
numSamples = 200000
for(idxR in 1:2){
    cat(sprintf("nomRange: %d", idxR))
    print(Sys.time())
    
    # Run through realizations
    jePrior[idxR, ] = foreach(iter = 1:numReal, .combine = c) %dopar%{
        # Generate a sample
        obs = genRel(loc, nomR = nomRange[idxR], varMod = varMod)
        
        # Make model
        modelObj = intPrior.Exp.makeModel(obs = obs, loc = loc)
        
        # Set prior
        modelObj = intPrior.Exp.addPrior(modelObj = modelObj, name = "Jeffrey")
        
        # Set target
        modelObj = intPrior.Exp.addTarget(modelObj)
        
        # Run sampler
        inits = c(log(nomRange[idxR]), 0)
        scale = c(0.6, 0.08)
        acc.rate = 0.4
        adapt = TRUE
        tmp = try(intPrior.addHyperSamples(modelObj = modelObj, numSamples = numSamples, inits = inits, scale = scale, acc.rate = acc.rate, adapt = adapt))
        
        if(class(tmp) == 'try-error'){
            list(res = list(NULL))
        } else {
            # Remove burn-in
            tmp$hyper$samples = tmp$hyper$samples[seq(25000, numSamples, 1),]
            tmp$hyper$log.p = tmp$hyper$log.p[seq(25000, numSamples, 1)]
            
            # Quantile intervals
            rRange = quantile(tmp$hyper$samples[,1], probs = c(0.025, 0.975))
            sRange = quantile(tmp$hyper$samples[,2], probs = c(0.025, 0.975))
            rMedian = median(tmp$hyper$samples[,1])
            sMedian = median(tmp$hyper$samples[,2])
            rTrue  = nomRange[idxR]
            sTrue  = sqrt(varMod)
            
            # Different highest posterior density intervals
            hpds = list(logRange  = boa.hpd(tmp$hyper$samples[, 1], 0.05),
                        logStdDev = boa.hpd(tmp$hyper$samples[, 2], 0.05),
                        range  = boa.hpd(exp(tmp$hyper$samples[, 1]), 0.05),
                        stdDev = boa.hpd(exp(tmp$hyper$samples[, 2]), 0.05),
                        var    = boa.hpd(exp(tmp$hyper$samples[, 2]*2), 0.05))
            
            # Compose results
            list(res = list(rRange = rRange, sRange = sRange, 
                            rTrue = rTrue, sTrue = sTrue,
                            rMedian = rMedian, sMedian = sMedian,
                            hpds = hpds,
                            obs = obs))
            }
    }
    try(save(jePrior, file = 'Results/jePrior.RData'))
}
try(save(jePrior, file = 'Results/jePrior.RData'))

# Tables for quantile-based credible intervals
jeCovTableR = array(0, dim = c(2))
jeLenTableR = array(0, dim = c(2))
jeCovTableV = array(0, dim = c(2))
jeLenTableV = array(0, dim = c(2))

# Tables for hpd credible intervals
jeCovTableR_hpd = array(0, dim = c(2))
jeLenTableR_hpd = array(0, dim = c(2))
jeCovTableV_hpd = array(0, dim = c(2))
jeLenTableV_hpd = array(0, dim = c(2))

for(ra in 1:2){
    count = 0
    for(real in 1:1500){
        if(is.null(jePrior[ra, real][[1]][[1]])){
            next
        }
        count = count + 1
        
        # Endpoints
        R = jePrior[ra, real][[1]]$rRange
        V = 2*jePrior[ra, real][[1]]$sRange
        R_hpd = jePrior[[ra, real]]$hpds$range
        V_hpd = jePrior[[ra, real]]$hpds$var
        
        # Lengths of intervals
        jeLenTableR[ra] = jeLenTableR[ra] + diff(exp(R))
        jeLenTableV[ra] = jeLenTableV[ra] + diff(exp(V))
        jeLenTableR_hpd[ra] = jeLenTableR_hpd[ra] + diff(R_hpd)
        jeLenTableV_hpd[ra] = jeLenTableV_hpd[ra] + diff(V_hpd)
        
        # Check if inside quantile-based intervals
        if((exp(R[1]) < nomRange[ra]) && (nomRange[ra] < exp(R[2]))){
            jeCovTableR[ra] = jeCovTableR[ra] + 1
        }

        if((exp(V[1]/2) < stdDev) && (stdDev < exp(V[2]/2))){
            jeCovTableV[ra] = jeCovTableV[ra] + 1
        }
        
        # Check if inside hpd intervals
        if((R_hpd[1] < nomRange[ra]) && (nomRange[ra] <R_hpd[2])){
            jeCovTableR_hpd[ra] = jeCovTableR_hpd[ra] + 1
        }
        
        if((V_hpd[1] < varMod) && (varMod < V_hpd[2])){
            jeCovTableV_hpd[ra] = jeCovTableV_hpd[ra] + 1
        }
        if(count == 1000)
            break
    }
    
    # Quantile-based intervals
    jeCovTableR[ra] = jeCovTableR[ra]/count
    jeLenTableR[ra] = jeLenTableR[ra]/count
    jeCovTableV[ra] = jeCovTableV[ra]/count
    jeLenTableV[ra] = jeLenTableV[ra]/count
    
    # HPD intervals
    jeCovTableR_hpd[ra] = jeCovTableR_hpd[ra]/count
    jeLenTableR_hpd[ra] = jeLenTableR_hpd[ra]/count
    jeCovTableV_hpd[ra] = jeCovTableV_hpd[ra]/count
    jeLenTableV_hpd[ra] = jeLenTableV_hpd[ra]/count
}
try(save(jeCovTableR, jeLenTableR, jeCovTableV, jeLenTableV, 
         jeCovTableR_hpd, jeLenTableR_hpd, jeCovTableV_hpd, jeLenTableV_hpd,
         jePrior, file = 'Results/jePrior.RData'))
