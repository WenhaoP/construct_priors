#########################################################################
# MCMCsampler.R                                                         #
#    Use adaptMCMC package to do sampling with different likelihoods    #
#    and priors.                                                        #
#########################################################################

# Libraries used by functions
require(INLA)
require(adaptMCMC)

## General sampling functions
	# Get samples of hyperparameters
	intPrior.addHyperSamples = function(modelObj,
										numSamples,
										inits,
										scale,
										acc.rate,
										adapt = TRUE){
		# Run
		hmm = MCMC(p = modelObj$target,
				   n = numSamples, 
				   init = inits, 
				   scale = scale, 
				   adapt = adapt,
				   acc.rate = acc.rate,
				   modelObj = modelObj)

		# Add result
		modelObj$hyper = hmm

		return(modelObj)
	}

## Exponential model
	# Log-target exponential model
	intPrior.Exp.logLikelihood = function(theta, modelObj){
		# Convert to correct scale
		nominalRange = exp(theta[1])
		processStdDev = exp(theta[2])

		# Set parameters
		nu = 0.5
		phi = sqrt(8*nu)/nominalRange
		varMod = processStdDev^2

		# Construct covariance matrix
		D = as.matrix(dist(modelObj$lik$data$loc))
		Sigma = varMod*exp(-phi*D)

		# Compute log-lihelihood
		obs = modelObj$lik$data$obs
		res = -0.5*length(obs)*log(2*pi)-0.5*determinant(Sigma, logarithm=TRUE)$modulus[1]-0.5*t(obs)%*%solve(Sigma, obs)

		return(res)
	}

	# Set likelihood for exponential model
	intPrior.Exp.makeModel = function(obs,
									  loc){
		# Likelihood
		data = list(obs = obs, loc = loc)
		fun  = intPrior.Exp.logLikelihood
		lik  = list(data = data, fun = fun)

		# Store in object
		modelObj = list(lik = lik)

		return(modelObj)
	}

	# Add prior for exponential model
	intPrior.Exp.addPrior = function(modelObj,
									 name,
									 pars = NULL){
		if(name == 'PC-prior'){
			# Data for prior
			lam1 = -pars$rho0*log(pars$a1)
			lam2 = -log(pars$a2)/pars$sig0
			data = list(lam1 = lam1, lam2 = lam2)

			# Function for prior
			fun = function(theta, modelObj){
				# Extract parameters
				nomR = exp(theta[1])
				stdDev = exp(theta[2])

				# Range prior
				lam1 = modelObj$prior$data$lam1
				logPrior = log(lam1)-2*log(nomR)-lam1*nomR^(-1)

				# Standard deviation
				lam2 = modelObj$prior$data$lam2
				logPrior = logPrior + log(lam2) - lam2*stdDev

				# Add jacobian of transformation
				logPrior = logPrior + theta[1] + theta[2]

				return(logPrior)
			}
		}
		if(name == 'Jeffrey'){
			# Data for prior
			data = list(NULL)

			# Function for prior
			fun = function(theta, modelObj){
				# Extract data
				loc = modelObj$lik$data$loc
				obs = modelObj$lik$data$obs

				# Extract parameters
				nomR = exp(theta[1])
				stdDev = exp(theta[2])

				# Correlation matrix
				nu = 0.5
				D = as.matrix(dist(loc))
				Sigma = exp(-sqrt(8*nu)/nomR*D)

				# Derivative
				SigmaDer = Sigma*(sqrt(8*nu)/nomR^2*D)

				# U
				U = t(solve(Sigma, SigmaDer))

				# Prior for std.dev.
				logPrior = -log(stdDev)

				# Prior for nominal range
				numObs = length(obs)
				logPrior = logPrior + 0.5*log(sum(sum(U*t(U)))-sum(diag(U))^2/numObs)

				# Add jacobian of transformation
				logPrior = logPrior + theta[1] + theta[2]

				return(logPrior)
			}
		}
		if(name == 'Uniform'){
			# Data
			data = list(L = pars$lower, U = pars$upper)

			# Function
			fun = function(theta, modelObj){
				# Extract parameters
				nomR = exp(theta[1])
				stdDev = exp(theta[2])

				# Range
				L = modelObj$prior$data$L
				U = modelObj$prior$data$U
				if((nomR < L) || (nomR > U)){
					logPrior = -Inf
				} else {
					logPrior = 0
				}

				# Standard deviation
				logPrior = logPrior - log(stdDev)

				# Correct for transformation
				logPrior = logPrior + theta[1] + theta[2]

				return(logPrior)
			}
		}
		if(name == 'Log-uniform'){
			# Data
			data = list(L = pars$lower, U = pars$upper)

			# Function
			fun = function(theta, modelObj){
				# Extract parameters
				nomR = exp(theta[1])
				stdDev = exp(theta[2])

				# Range
				L = modelObj$prior$data$L
				U = modelObj$prior$data$U
				if((nomR < L) || (nomR > U)){
					logPrior = -Inf
				} else {
					logPrior = -log(nomR)
				}

				# Standard deviation
				logPrior = logPrior - log(stdDev)

				# Correct for transformation
				logPrior = logPrior + theta[1] + theta[2]

				return(logPrior)
			}
		}

		# Store prior
		modelObj$prior = list(data = data, fun = fun)

		return(modelObj)
	}

	# Add target function
	intPrior.Exp.addTarget = function(modelObj){

		modelObj$target = function(theta, modelObj){
			logPost = modelObj$prior$fun(theta, modelObj) + modelObj$lik$fun(theta, modelObj)
			return(logPost)
		}
		
		return(modelObj)
	}	


## Logistic Exponential model
	# Log-target exponential model
	intPrior.Log.logLikelihood = function(theta, modelObj){
		# Extract spatial field
		uSpace = theta[-c(1,2)]
		loc = modelObj$lik$data$loc

		# Calculate pi(uSpace| nomR, stdDev)
		tmpObj = intPrior.Exp.makeModel(obs = uSpace, loc = loc)
		logLik = intPrior.Exp.logLikelihood(theta, tmpObj)

		# Calculate pi(y | uSpace, n)
		probs = pnorm(uSpace)
		obs = modelObj$lik$data$obs
		size = modelObj$lik$data$n
		logLik = logLik + sum(dbinom(obs, size, prob = probs, log = TRUE))

		return(logLik)
	}

	# Set likelihood for exponential model
	intPrior.Log.makeModel = function(obs,
									  n,
									  loc){
		# Likelihood
		data = list(obs = obs, n = n, loc = loc)
		fun  = intPrior.Log.logLikelihood
		lik  = list(data = data, fun = fun)

		# Store in object
		modelObj = list(lik = lik)

		return(modelObj)
	}

	# Add prior for exponential model
	intPrior.Log.addPrior = function(modelObj,
									 name,
									 pars = NULL){
		if(name == 'PC-prior'){
			# Data for prior
			lam1 = -pars$rho0*log(pars$a1)
			lam2 = -log(pars$a2)/pars$sig0
			data = list(lam1 = lam1, lam2 = lam2)

			# Function for prior
			fun = function(theta, modelObj){
				# Extract parameters
				nomR = exp(theta[1])
				stdDev = exp(theta[2])

				# Range prior
				lam1 = modelObj$prior$data$lam1
				logPrior = log(lam1)-2*log(nomR)-lam1*nomR^(-1)

				# Standard deviation
				lam2 = modelObj$prior$data$lam2
				logPrior = logPrior + log(lam2) - lam2*stdDev

				# Add jacobian of transformation
				logPrior = logPrior + theta[1] + theta[2]

				return(logPrior)
			}
		}
		if(name == 'Jeffrey'){
			# Data for prior
			data = list(NULL)

			# Function for prior
			fun = function(theta, modelObj){
				# Extract data
				loc = modelObj$lik$data$loc
				obs = modelObj$lik$data$obs

				# Extract parameters
				nomR = exp(theta[1])
				stdDev = exp(theta[2])

				# Correlation matrix
				nu = 0.5
				D = as.matrix(dist(loc))
				Sigma = exp(-sqrt(8*nu)/nomR*D)

				# Derivative
				SigmaDer = Sigma*(sqrt(8*nu)/nomR^2*D)

				# U
				U = t(solve(Sigma, SigmaDer))

				# Prior for std.dev.
				logPrior = -log(stdDev)

				# Prior for nominal range
				numObs = length(obs)
				logPrior = logPrior + 0.5*log(sum(sum(U*t(U)))-sum(diag(U))^2/numObs)

				# Add jacobian of transformation
				logPrior = logPrior + theta[1] + theta[2]

				return(logPrior)
			}
		}
		if(name == 'Uniform'){
			# Data
			data = list(L = pars$lower, U = pars$upper)

			# Function
			fun = function(theta, modelObj){
				# Extract parameters
				nomR = exp(theta[1])
				stdDev = exp(theta[2])

				# Range
				L = modelObj$prior$data$L
				U = modelObj$prior$data$U
				if((nomR < L) || (nomR > U)){
					logPrior = -Inf
				} else {
					logPrior = 0
				}

				# Standard deviation
				logPrior = logPrior - log(stdDev)

				# Correct for transformation
				logPrior = logPrior + theta[1] + theta[2]

				return(logPrior)
			}
		}
		if(name == 'Log-uniform'){
			# Data
			data = list(L = pars$lower, U = pars$upper)

			# Function
			fun = function(theta, modelObj){
				# Extract parameters
				nomR = exp(theta[1])
				stdDev = exp(theta[2])

				# Range
				L = modelObj$prior$data$L
				U = modelObj$prior$data$U
				if((nomR < L) || (nomR > U)){
					logPrior = -Inf
				} else {
					logPrior = -log(nomR)
				}

				# Standard deviation
				logPrior = logPrior - log(stdDev)

				# Correct for transformation
				logPrior = logPrior + theta[1] + theta[2]

				return(logPrior)
			}
		}

		# Store prior
		modelObj$prior = list(data = data, fun = fun)

		return(modelObj)
	}

	# Add target function
	intPrior.Log.addTarget = function(modelObj){

		modelObj$target = function(theta, modelObj){
			logPost = modelObj$prior$fun(theta, modelObj) + modelObj$lik$fun(theta, modelObj)
			return(logPost)
		}
		
		return(modelObj)
	}	

## Non-stationary model
	# Function to center covariates \int f_i(s) ds = 0, for i = 1,...,p
	intPrior.nonStat.centerCovariates = function(Cov,
												 mesh){
		# Extract integration matrix
		C0 = inla.mesh.fem(mesh)$c0

		# Make it integrate to zero
		h = h - sum(C0%*%h)/sum(diag(C0))
	}

	# Function to calculate KLD between to finite-dim approximations
	# KLD for kappa, tau --> kappa(.), tau
	intPrior.nonStat.kld.kappa = function(xKappa, xStat, xTau, modelObj){
		# Extract SPDE object
		spde = modelObj$lik$data.inla$spde

		# Precision matrix for stationary model
		Q0 = inla.spde.precision(spde = spde, theta = c(xStat, 0*xKappa, 0*xTau))

		# Precision matrix for non-stationary model
		Q1 = inla.spde.precision(spde = spde, theta = c(xStat, xKappa, 0*xTau))

		# Formula for KLD is
		#  D_{KL}(NonStat(1)||Stat(0)) = 0.5*(tr(Q0*Q1^{-1})-N-log|Q0|/|Q1|)

		# Trace-N
			Qinv1 = inla.qinv(Q1)
			part1 = sum(Qinv1*Q0)-dim(Q0)[1]
		
		# Logarithm
			require(spam)
			det0 = determinant(Q0, logarithm = TRUE)$modulus
			det1 = determinant(Q1, logarithm = TRUE)$modulus
			part2 = -det0+det1
		
		# Full KLD
			KLD = 0.5*(part1+part2)	

		return(KLD)
	}

	# KLD for tau
	intPrior.nonStat.kld.tau = function(xTau, xStat, xKappa, modelObj){
		# Formula for KLD is (in the limit)
		#    KLD(tau(.) || tau0) = 0.5*\int((tau0/tau(s))^2-1-log(tau0/tau(s))^2) dA
		
		# Extract
		spde = modelObj$lik$data.inla$spde
		C0   = inla.mesh.fem(modelObj$lik$data.inla$mesh)$c0

		# Find \tau(.)
		tau = exp(spde$param.inla$B0[,-1]%*%c(0*xStat, xKappa, xTau))

		# Integrate
		KLD = 0.5*sum(C0%*%(tau^-2-1-log(tau^-2)))	

		return(KLD)
	}

	# Function to make SPDE object
	intPrior.nonStat.makeSPDE = function(mesh, covKappa, covTau, extraconstr.int = NULL){
		# Fix input
		if(!is.null(covKappa)){
			covKappa = as.matrix(covKappa)
		}
		if(!is.null(covTau)){
			covTau   = as.matrix(covTau)
		}

		# Make B matrices
		alpha = 2
		d = 2
		nu = alpha - d/2
		kappa0 = log(sqrt(8*nu))
		tau0   = 0.5*(lgamma(nu)-lgamma(nu+d/2)-d/2*log(4*pi))-nu*kappa0
		B.tau   = cbind(tau0,   nu, -1, 0*covKappa,   covTau)
		B.kappa = cbind(kappa0, -1,  0,   covKappa, 0*covTau)

		# Assemble to spde object for INLA
		spde   = inla.spde2.matern(mesh = mesh,
								   alpha = alpha,
								   B.tau   = B.tau,
								   B.kappa = B.kappa,
								   extraconstr.int = extraconstr.int) 

		return(spde)
	}

	# Function to collect covariates and observations
	intPrior.nonStat.makeModel = function(obs,
										  loc,
										  mesh,
										  linCov,
										  covCov){
		# Check that library is loaded
		require(INLA)

		# Make the A matrix
		A = inla.spde.make.A(mesh = mesh, loc = loc)

		# Make SPDE object
		spde = intPrior.nonStat.makeSPDE(mesh, covCov, covCov)

		# Fix input
		covKappa = covCov
		if(!is.null(covKappa)){
			covKappa = as.matrix(covKappa)
		}
		covTau = covCov
		if(!is.null(covTau)){
			covTau   = as.matrix(covTau)
		}


		# Store data
		data = list(obs = obs,
					loc = loc,
					linCov = linCov,
					covKappa = covKappa,
					covTau = covTau)

		# Store SPDE stuff
		data.inla = list(mesh = mesh,
						 A = A,
						 spde = spde)

		# Put everything in target
		fun = NULL

		# Store everything
		modelObj = list(lik = list(data = data,
								   data.inla = data.inla,
								   fun = fun))

		return(modelObj)
	}

	# Prior for stationary part
	intPrior.nonStat.priorStat = function(x, par){
		# Prior for log(range)
		logRange = x[1]
		d = 2
		lam1 = par[1]
		logPrior1 = log(lam1*d/2) - (1+d/2)*logRange - lam1*exp(-d/2*logRange) + logRange

		# Prior for log(Std.Dev.)
		logStdDev = x[2]
		lam2 = par[2]
		logPrior2 = log(lam2)-lam2*exp(logStdDev)+logStdDev

		return(logPrior1+logPrior2)
	}

	# Prior for non-stationary part
	intPrior.nonStat.priorNonStat = function(xStat, xKappa, xTau, modelObj){
		# Extract stuff
		lookUp = modelObj$prior$data$nonStat$lookUp
		lambdaKappa = modelObj$prior$data$nonStat$lamKappa
		lambdaTau   = modelObj$prior$data$nonStat$lamTau

		if(length(xKappa) > 0){
		## Log(kappa(.))
			# Second-order approximation of KLD for
			#    kappa, tau --> kappa(.), tau
			#	require(numDeriv)
				#lookUp = NULL
				if(is.null(lookUp)){
					hess = hessian(intPrior.nonStat.kld.kappa,
								   x = 0*xKappa,
								   xStat = xStat,
								   xTau = xTau,
								   modelObj = modelObj)
				} else{
					idx = length(which(lookUp$x < xStat[1]))
					if(idx < 1){
						hess = lookUp$H[[1]]
					} else if(idx > 39){
						hess = lookUp$H[[40]]
					} else{
						d1 = xStat[1]-lookUp$x[idx]
						d2 = lookUp$x[idx+1]-xStat[1]
						hess = (lookUp$H[[idx]]*d2+lookUp$H[[idx+1]]*d1)/(d1+d2)
					}
				}

				# Use second-order approximation
				prior.kappa = inla.pc.multvar.sphere.d(x = xKappa,
													   lambda = lambdaKappa,
													   log = TRUE,
													   H = hess)
		} else{
			prior.kappa = 0
		}


		if(length(xTau) > 0){
		## Log(tau(.))
			# Second-order approximation of KLD for
			#    kappa(.), tau --> kappa(.), tau(.)
		#	hess = hessian(intPrior.nonStat.kld.tau,
		#				   x = 0*xTau,
		#				   xStat = xStat,
		#				   xKappa = xKappa,
		#				   modelObj = modelObj,
		#				   method.args = list(eps = 1e-1))
			hess = modelObj$prior$data$nonStat$hessTau

			# Use second-order approximation
			prior.tau = inla.pc.multvar.sphere.d(x = xTau,
												 lambda = lambdaTau,
												 log = TRUE,
												 H = hess)
		} else {
			prior.tau = 0
		}

		#print(paste('kappa: ', prior.kappa, '   tau: ', prior.tau, '   ', lambdaKappa, '    ', lambdaTau, sep = ""))
		return(prior.kappa+prior.tau)
	}

	# Function to calculate KLD prior
	intPrior.nonStat.prior = function(x, modelObj){
		# Load inla
		require(INLA)
	
		## Prior on nugget (log-precision)
		##    uN and aN loaded from file
			uN = modelObj$prior$data$noise$U
			aN = modelObj$prior$data$noise$alpha
			tau = exp(x[1])
			prior.nugget = inla.pc.dprec(tau, u = uN, alpha = aN, log = TRUE)+log(tau)
			x = x[-1]
			#prior.nugget = 0

		## Prior on stationary SPDE (log(range) and log(Std.Dev.)
			# Parameter for log(range)
			lR = modelObj$prior$data$range$rho0
			aR = modelObj$prior$data$range$alpha
			d = 2
			lambdaRange =-lR^(d/2)*log(aR)

			# Parameter for log(std.dev.)
			uS = modelObj$prior$data$stdDev$sigma0
			aS = modelObj$prior$data$stdDev$alpha
			lambdaStdDev = -log(aS)/uS

			# Penalisation rates for stationary model
			parStat = c(lambdaRange, lambdaStdDev)

			# Calculate stationary
			xStat = x[1:2]
			prior.stat = intPrior.nonStat.priorStat(xStat, parStat)
			x = x[-(1:2)]

		## Prior on non-stationary SPDE
			# Parameters for log(kappa(.)) and log(tau(.))

			# Calculate non-stationary
			covKappa = modelObj$lik$data$covKappa
			covTau = modelObj$lik$data$covTau
			if(!is.null(covKappa)){
				xKappa = x[1:dim(covKappa)[2]]
				x = x[-(1:dim(covKappa)[2])]
			} else {
				xKappa = NULL
			}
			if(!is.null(covTau)){
				xTau = x[1:dim(covTau)[2]]
				x = x[-(1:dim(covTau)[2])]
			} else{
				xTau = NULL
			}
			prior.nonStat = intPrior.nonStat.priorNonStat(xStat, xKappa, xTau, modelObj)

		## Assemble
			logDens = prior.nugget + prior.stat + prior.nonStat
			return(logDens)
	}

	# Wrapper for INLA
	intPrior.nonStat.prior.INLA = function(x){
		load('tmpModelObj.RData')
		return(intPrior.nonStat.prior(x, modelObj))
	}

	# Function to collect priors in a prior object
	intPrior.nonStat.addPrior = function(modelObj,
										 fixedPar,
										 noisePar,
										 rangePar,
										 stdDevPar,
										 nonStatPar){
		if(!is.null(modelObj$lik$data$covKappa)){
		## Precompute Hessian
			H = list()

			# Parameter in correlation
			## TODO: User specify this!
			xs = seq(-2, 10, length.out = 40)
			require(numDeriv)

			# Do all computations
			numKappa = dim(modelObj$lik$data$covKappa)[2]
			numTau   = dim(modelObj$lik$data$covTau)[2]
			for(idx in 1:length(xs)){
				hess = hessian(intPrior.nonStat.kld.kappa,
						   	   x = rep(0, numKappa),
						       xStat = c(xs[idx], 0),
						       xTau = rep(0, numTau),
						       modelObj = modelObj,
						       method.args = list(eps = 1e-1))
				print(idx)
				print(hess)
				H[[idx]] = hess
			}
			lookUp = list(x = xs, H = H)
		} else{
			lookUp = NULL
		}


		if(!is.null(modelObj$lik$data$covTau)){

			#    kappa(.), tau --> kappa(.), tau(.)
			hessTau = hessian(intPrior.nonStat.kld.tau,
							  x = rep(0, numTau),
						   	  xStat = c(0,0),
						   	  xKappa = rep(0, numKappa),
						   	  modelObj = modelObj,
						   	  method.args = list(eps = 1e-1))

		} else {
			hessTau = NULL
		}

		
		# Set prior function
		fun = intPrior.nonStat.prior

		# Collect everything
		data = list(fixed = fixedPar,
					noise = noisePar,
					range = rangePar,
					stdDev = stdDevPar,
					nonStat = c(nonStatPar, list(lookUp = lookUp, hessTau = hessTau)))

		# Add to model object
		modelObj$prior = list(data = data,
							  fun = fun)

		return(modelObj)
	}

	# Function to collect priors in a prior object
	intPrior.addTarget = function(modelObj,
								  constr = NULL){

		if(is.null(constr)){
			modelObj$constr = NULL
			modelObj$target = nonStat.logTarget
		} else {
			# Extract integration matrix
			C0 = inla.mesh.fem(modelObj$lik$data.inla$mesh)$c0

			# Set constraints
			modelObj$constr = C0%*%constr

			# Constraint is on spatial field
			pCol = dim(constr)[2]
			pRow = dim(modelObj$lik$data$linCov)[2]
			modelObj$constr = rBind(modelObj$constr, Matrix(0, ncol = pCol, nrow = pRow))

			# Correct rotation
			modelObj$constr = t(modelObj$constr)
					
			# Target
			modelObj$target = nonStat.logTarget
		}

		return(modelObj)
	}

	logTargetINLA = function(x, modelObj){
		# Prior on fixed effects
		tauB = modelObj$prior$data$fixed$tauB

		# Calculate model distribution
		theta = x[-1]
		Qxx = INLA:::inla.spde2.precision(spde = modelObj$spde, theta = theta)
		p = dim(modelObj$linCov)[2]
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

	# Target
	nonStat.logTarget = function(x, modelObj){
		# Prior on fixed effects
		tauB = modelObj$prior$data$fixed$tauB

		# Calculate model distribution
		theta = x[-1]
		Qxx = INLA:::inla.spde.precision(spde = modelObj$lik$data.inla$spde, theta = theta)
		p = dim(modelObj$lik$data$linCov)[2]
		m = dim(Qxx)[1]
		Qbb = diag(rep(tauB, p), ncol = p, nrow = p)
		Qxb = Matrix(0, ncol = p, nrow = m)
		Q = rBind(cBind(Qxx, Qxb), cBind(t(Qxb), Qbb))

		# Calculate A matrix
		A = modelObj$lik$data.inla$A
		A = cBind(A, modelObj$lik$data$linCov)

		# Calculate conditional distribution
		tau0 = exp(x[1])
		Qc = Q + tau0*t(A)%*%A

		# Precompute cholesky factors
		L = Cholesky(Q)
		Lc = Cholesky(Qc)

		# Find conditional mean
		obs = modelObj$lik$data$obs
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
		lPrior = modelObj$prior$fun(x, modelObj = modelObj)

		# Total density
		lDens = lLike+lPrior

		return(as.numeric(lDens))
	}


	# Get prediction samples
	intPrior.addPredSamples = function(modelObj, loc, covValues, thinning){
		# Convert locations to mesh combinations
		Apred = inla.spde.make.A(mesh = modelObj$lik$data.inla$mesh, loc = loc)
		Apred = cbind(Apred, covValues)

		# Extract thetas
		thetaAll = modelObj$hyper$samples

		# Thin
		thetaAll = thetaAll[seq(1, dim(thetaAll)[1], thinning), ]
		numSamples = dim(thetaAll)[1]

		# Store in
		numPred = dim(loc)[1]
		modelObj$pred = list(samples = matrix(0, nrow = numSamples, ncol = numPred))
		modelObj$pred$samplesNN = modelObj$pred$samples

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

				# Prior on fixed effects
				tauB = modelObj$prior$data$fixed$tauB

				# Calculate model distribution
				theta = x[-1]
				Qxx = inla.spde.precision(spde = modelObj$lik$data.inla$spde, theta = theta)
				p = dim(modelObj$lik$data$linCov)[2]
				m = dim(Qxx)[1]
				Qbb = diag(rep(tauB, p), ncol = p, nrow = p)
				Qxb = Matrix(0, ncol = p, nrow = m)
				Q = rBind(cBind(Qxx, Qxb), cBind(t(Qxb), Qbb))

				# Calculate A matrix
				A = cBind(modelObj$lik$data.inla$A, modelObj$lik$data$linCov)

				# Calculate conditional distribution
				tau0 = exp(x[1])
				Qc = Q + tau0*t(A)%*%A

				obs = modelObj$lik$data$obs
				numLatent = dim(Q)[1]

				# Precompute cholesky factors
				Lc = Cholesky(Qc, LDL = FALSE)

				# Find conditional mean
				muC = solve(Lc, tau0*t(A)%*%obs)

				# Permute
				muCperm = muC[Lc@perm+1,,drop=FALSE]

				# Sample
				meshValues = as.vector(solve(Lc, rnorm(numLatent), system = 'Lt') + muCperm)

				# Permute back
				meshValues = meshValues[invPerm(Lc@perm+1)]

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

				# Convert to prediction locations
				meshValues = as.vector(Apred%*%meshValues)
				
				# Calculate values at prediction locations without noise
				modelObj$pred$samplesNN[idx,] = meshValues

				# Predict with noise
				modelObj$pred$samples[idx,] = modelObj$pred$samplesNN[idx,] + rnorm(numPred)/sqrt(tau0)

			}
		}

		return(modelObj)
	}

	# Post-calculate leave-one-out predictions
	intPrior.LOOLcalc = function(modelObj, WAIC = FALSE){
		# Extract stuff
		obs = modelObj$lik$data$obs
		hSamples = modelObj$hyper$samples
		pSamples = modelObj$pred$samplesNN

		# Get dimensions
		modelObj$pred$LOOLpost = modelObj$pred$samples*0
		numSamples = dim(modelObj$pred$LOOLpost)[1]
		numPred = dim(modelObj$pred$LOOLpost)[2]

		# Calculate all densities \pi(y_i| x^s, \theta^s), where
		# (x^s, \theta^s) are drawn from \pi(x, \theta | y)
		for(oIdx in 1:numPred){
			for(idx in 1:numSamples){
				modelObj$pred$LOOLpost[idx, oIdx] = dnorm(obs[oIdx], mean = pSamples[idx, oIdx], sd = exp(-0.5*hSamples[idx,1]))
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

		## LOO score		
			modelObj$pred$LOOLpost = log(1/colMeans(1/modelObj$pred$LOOLpost))

		res = modelObj
	}

	## Check quality of samples
	checkQuality = function(res, burnin = 10000){
		cat(paste('Note: Using burn-in = ', burnin, '\n\n', sep = ''))
		for(idx in 1:length(res)){
			# Check if failed/success
			if(is.null(res[[idx]])){
				cat(paste('Model ', idx-1, ' failed!\n', sep = ''))
				next
			} else {
				cat(paste('Model ', idx-1, ':\n', sep = ''))
			}

			# Extract samples
			samples = res[[idx]]$hyper$samples[-c(1:burnin),]

			# Check effective sizes
			require(coda)
			eff = effectiveSize(samples)
			cat('Effective sizes:\n')
			for(par in 1:dim(samples)[2]){
				cat(paste('\tParameter ', par, ': ', sprintf('%.0f', eff[par]), '\n', sep = ''))
			}
			cat('\n')
		}
	}

	## Overlay different marginal posteriors of hyperparameters
	plotHyperPosterior = function(res, idx, models, new.dev = FALSE){
		if(new.dev){
			dev.new()
		}

		# Calculate densities
		dens = list()
		count = 1
		xlim = NULL
		ylim = NULL
		for(mod in models){
			dens[[count]] = density(res[[mod]]$hyper$samples[,idx])
			xlim = range(xlim, dens[[count]]$x)
			ylim = range(ylim, dens[[count]]$y)
			count = count + 1
		}
		
		plot(NULL, xlim = xlim, ylim = ylim)
		count = 1
		for(mod in models){
			lines(dens[[count]], col = count)
			count = count+1
		}
	}

	## Overlay different marginal posteriors of latent field
	plotLatentPosterior = function(res, idx, models, xlimOver = NULL, ylimOver = NULL, new.dev = FALSE){
		if(new.dev){
			dev.new()
		}

		# Calculate densities
		dens = list()
		count = 1
		xlim = NULL
		ylim = NULL
		for(mod in models){
			dens[[count]] = density(res[[mod]]$pred$samples[,idx])
			xlim = range(xlim, dens[[count]]$x)
			ylim = range(ylim, dens[[count]]$y)
			count = count + 1
		}
		
		if(!is.null(xlimOver)){
			xlim = xlimOver
		}
		if(!is.null(ylimOver)){
			ylim = ylimOVer
		}
		plot(NULL, xlim = xlim, ylim = ylim)
		count = 1
		for(mod in models){
			lines(dens[[count]], col = count)
			count = count+1
		}
	}
