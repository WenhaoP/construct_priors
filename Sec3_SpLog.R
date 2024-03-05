################################################################################
# Sec3_SpLog.R                                                                 #
#    Example with spatial logistic regression                                  #
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
nomRange = c(0.1)
size = 20

# Number of realizations
numReal = 1500

# Initialize storage
PCpriorLR  = array(list(NULL), dim = c(1, 4, 4, numReal))

# Run through settings
numSamples = 1e6
for(idxR in c(1)){
    for(idxL1 in 1:4){
        for(idxL2 in 1:4){
            # Prior: range
            a1  = 0.05
            rho0 = nomRange[idxR]*c(0.025, 0.1, 0.4, 1.6)[idxL1]
            
            # Prior: std.dev.
            a2 = 0.05
            sig0 = sqrt(varMod)*c(0.625,2.5,10, 40)[idxL2]
            
            # Current setting
            cat(sprintf("nomRange: %d\nLam1: %d\nLam2: %d\n", idxR, idxL1, idxL2))
            print(Sys.time())
            
            # Run through realizations
            PCpriorLR[idxR, idxL1, idxL2, ] = foreach(iter = 1:numReal, .combine = c) %dopar%{
                # Generate a sample
                tmpObs = genRelLogistic(loc, nomR = nomRange[idxR], varMod = varMod, size = size)
                obs = tmpObs$log.obs
                
                # Make model
                modelObj = intPrior.Log.makeModel(obs = obs, n = size, loc = loc)
                
                # Set prior
                pars = list(a1 = a1, rho0 = rho0, a2 = a2, sig0 = sig0)
                modelObj = intPrior.Log.addPrior(modelObj = modelObj, name = "PC-prior", pars = pars)
                
                # Set target
                modelObj = intPrior.Log.addTarget(modelObj)
                
                # Run sampler
                inits = c(log(nomRange[idxR]), 0, tmpObs$spatial.obs)
                scale = c(0.08, 0.002, rep(0.01, 25))
                acc.rate = 0.3
                adapt = FALSE
                tmp = try(intPrior.addHyperSamples(modelObj = modelObj, numSamples = numSamples, inits = inits, scale = scale, acc.rate = acc.rate, adapt = adapt))
                
                if(class(tmp) == 'try-error'){
                    list(res = list(NULL))
                } else {
                    # Remove burn-in and thin
                    tmp$hyper$samples = tmp$hyper$samples[seq(1e5, numSamples, 1),]
                    tmp$hyper$log.p = tmp$hyper$log.p[seq(1e5, numSamples, 1)]
                    
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
                                    obs = obs,
                                    spatial.obs = tmpObs$spatial.obs))
                }
            }	
            try(save(PCpriorLR, file = 'Results/PCpriorLR.RData'))
        }
    }
}
try(save(PCpriorLR, file = 'Results/PCpriorLR.RData'))

# Store quantile-based results
pcCovTableR = array(0, dim = c(1, 4, 4))
pcLenTableR = array(0, dim = c(1, 4, 4))
pcCovTableV = array(0, dim = c(1, 4, 4))
pcLenTableV = array(0, dim = c(1, 4, 4))

# Store HPD results
pcCovTableR_hpd = array(0, dim = c(1, 4, 4))
pcLenTableR_hpd = array(0, dim = c(1, 4, 4))
pcCovTableV_hpd = array(0, dim = c(1, 4, 4))
pcLenTableV_hpd = array(0, dim = c(1, 4, 4))

for(ra in 1:1){
    for(row in 1:4){
        for(col in 1:4){
            count = 0
            for(real in 1:1500){
                if(is.null(PCpriorLR[ra, row, col, real][[1]][[1]])){
                    next
                }
                count = count + 1
                
                # Quantile-based credible intervals
                R = PCpriorLR[ra, row, col, real][[1]]$rRange
                V = 2*PCpriorLR[ra, row, col, real][[1]]$sRange
                
                # HPD credible intervals
                R_hpd = PCpriorLR[ra, row, col, real][[1]]$hpds$range
                V_hpd = PCpriorLR[ra, row, col, real][[1]]$hpds$var
                
                # Lengths
                pcLenTableR[ra, row, col] = pcLenTableR[ra, row, col] + diff(exp(R))
                pcLenTableV[ra, row, col] = pcLenTableV[ra, row, col] + diff(exp(V))
                pcLenTableR_hpd[ra, row, col] = pcLenTableR_hpd[ra, row, col] + diff(R_hpd)
                pcLenTableV_hpd[ra, row, col] = pcLenTableV_hpd[ra, row, col] + diff(V_hpd)
                
                # Coverage of quantile based credible intervals
                if((exp(R[1]) < nomRange[ra]) && (nomRange[ra] < exp(R[2]))){
                    pcCovTableR[ra, row, col] = pcCovTableR[ra, row, col] + 1
                }
                if((exp(V[1]/2) < stdDev) && (stdDev < exp(V[2]/2))){
                    pcCovTableV[ra, row, col] = pcCovTableV[ra, row, col] + 1
                }
                
                # Coverage of HPD credible intervals
                if((R_hpd[1] < nomRange[ra]) && (nomRange[ra] < R_hpd[2])){
                    pcCovTableR_hpd[ra, row, col] = pcCovTableR_hpd[ra, row, col] + 1
                }
                if((V_hpd[1] < varMod) && (varMod < V_hpd[2])){
                    pcCovTableV_hpd[ra, row, col] = pcCovTableV_hpd[ra, row, col] + 1
                }
                
                if(count == 1000){
                    break
                }
            }
            
            # Quantile-based intervals
            pcCovTableR[ra, row, col] = pcCovTableR[ra, row, col]/count
            pcLenTableR[ra, row, col] = pcLenTableR[ra, row, col]/count
            pcCovTableV[ra, row, col] = pcCovTableV[ra, row, col]/count
            pcLenTableV[ra, row, col] = pcLenTableV[ra, row, col]/count
            
            # HPD intervals
            pcCovTableR_hpd[ra, row, col] = pcCovTableR_hpd[ra, row, col]/count
            pcLenTableR_hpd[ra, row, col] = pcLenTableR_hpd[ra, row, col]/count
            pcCovTableV_hpd[ra, row, col] = pcCovTableV_hpd[ra, row, col]/count
            pcLenTableV_hpd[ra, row, col] = pcLenTableV_hpd[ra, row, col]/count
        }
    }
}

try(save(pcCovTableR, pcLenTableR, pcCovTableV, pcLenTableV,
         pcCovTableR_hpd, pcLenTableR_hpd, pcCovTableV_hpd, pcLenTableV_hpd,
         PCpriorLR, file = 'Results/PCpriorLR.RData'))
