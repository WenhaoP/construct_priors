################################################################################
# Sec3_PriorUn1.R                                                              #
#    Simulation study for PriorUn1                                             #
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
un1Prior = array(list(NULL), dim = c(2, 3, 3, numReal))

# Run through settings
numSamples = 50000
for(idxR in 1:2){
    for(idxLo in 1:3){
        for(idxUp in 1:3){
            # Set bounds for range
            lLim = c(0.05, 0.005, 0.0005)[idxLo]
            uLim = c(2, 20, 200)[idxUp]
            par = list(lower = lLim, upper = uLim)
            # Current setting
            cat(sprintf("nomRange: %d\nLam1: %d\nLam2: %d\n", idxR, idxLo, idxUp))
            print(Sys.time())
            
            # Run through realizations
            un1Prior[idxR, idxLo, idxUp, ] = foreach(iter = 1:numReal, .combine = c) %dopar%{
                # Generate a sample
                obs = genRel(loc, nomR = nomRange[idxR], varMod = varMod)
                
                # Make model
                modelObj = intPrior.Exp.makeModel(obs = obs, loc = loc)
                
                # Set prior
                modelObj = intPrior.Exp.addPrior(modelObj = modelObj, name = "Uniform", pars = par)
                
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
                    # Remove burn-in and thin
                    tmp$hyper$samples = tmp$hyper$samples[seq(5000, numSamples, 1),]
                    tmp$hyper$log.p = tmp$hyper$log.p[seq(5000, numSamples, 1)]
                    
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
            try(save(un1Prior, file = 'Results/un1Prior.RData'))		
        }
    }
}
try(save(un1Prior, file = 'Results/un1Prior.RData'))

# Quantile-based tables
un1CovTableR = array(0, dim = c(2,3,3))
un1LenTableR = array(0, dim = c(2,3,3))
un1CovTableV = array(0, dim = c(2,3,3))
un1LenTableV = array(0, dim = c(2,3,3))

# HPD tables
un1CovTableR_hpd = array(0, dim = c(2,3,3))
un1LenTableR_hpd = array(0, dim = c(2,3,3))
un1CovTableV_hpd = array(0, dim = c(2,3,3))
un1LenTableV_hpd = array(0, dim = c(2,3,3))

for(ra in 1:2){
    for(row in 1:3){
        for(col in 1:3){
            count = 0
            for(real in 1:1500){
                if(is.null(un1Prior[ra, row, col, real][[1]][[1]])){
                    next
                }
                count = count + 1
                
                # Intervals
                R = un1Prior[ra, row, col, real][[1]]$rRange
                V = 2*un1Prior[ra, row, col, real][[1]]$sRange
                R_hpd = un1Prior[ra, row, col, real][[1]]$hpds$range
                V_hpd = un1Prior[ra, row, col, real][[1]]$hpds$var
                
                # Lengths
                un1LenTableR[ra, row, col] = un1LenTableR[ra, row, col] + diff(exp(R))
                un1LenTableV[ra, row, col] = un1LenTableV[ra, row, col] + diff(exp(V))
                un1LenTableR_hpd[ra, row, col] = un1LenTableR_hpd[ra, row, col] + diff(R_hpd)
                un1LenTableV_hpd[ra, row, col] = un1LenTableV_hpd[ra, row, col] + diff(V_hpd)
                
                # Quantile-based credible intervals
                if((exp(R[1]) < nomRange[ra]) && (nomRange[ra] < exp(R[2]))){
                    un1CovTableR[ra, row, col] = un1CovTableR[ra, row, col] + 1
                }
                
                if((exp(V[1]/2) < stdDev) && (stdDev < exp(V[2]/2))){
                    un1CovTableV[ra, row, col] = un1CovTableV[ra, row, col] + 1
                }
                
                # HPD credible intervals
                if((R_hpd[1] < nomRange[ra]) && (nomRange[ra] < R_hpd[2])){
                    un1CovTableR_hpd[ra, row, col] = un1CovTableR_hpd[ra, row, col] + 1
                }
                
                if((V_hpd[1] < varMod) && (varMod < V_hpd[2])){
                    un1CovTableV_hpd[ra, row, col] = un1CovTableV_hpd[ra, row, col] + 1
                }
                if(count == 1000){
                    break
                }
            }
            
            # Quantile-based intervals
            un1CovTableR[ra, row, col] = un1CovTableR[ra, row, col]/count
            un1LenTableR[ra, row, col] = un1LenTableR[ra, row, col]/count
            un1CovTableV[ra, row, col] = un1CovTableV[ra, row, col]/count
            un1LenTableV[ra, row, col] = un1LenTableV[ra, row, col]/count
            
            # HPD intervals
            un1CovTableR_hpd[ra, row, col] = un1CovTableR_hpd[ra, row, col]/count
            un1LenTableR_hpd[ra, row, col] = un1LenTableR_hpd[ra, row, col]/count
            un1CovTableV_hpd[ra, row, col] = un1CovTableV_hpd[ra, row, col]/count
            un1LenTableV_hpd[ra, row, col] = un1LenTableV_hpd[ra, row, col]/count
        }
    }
}
try(save(un1CovTableR, un1LenTableR, un1CovTableV, un1LenTableV,
         un1CovTableR_hpd, un1LenTableR_hpd, un1CovTableV_hpd, un1LenTableV_hpd,
         un1Prior, file = 'Results/un1Prior.RData'))
