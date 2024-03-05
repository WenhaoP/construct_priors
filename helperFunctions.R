
# Plot function
scatterplot = function(loc, val, ...){
	colors = rev(jet.colors(1001))
	zcolor = colors[(val-min(val))/diff(range(val))*1000+1]
	plot(loc[,1],loc[,2], col=zcolor, ...)
}

# Generate realization from spatial field
genRel = function(loc, nomR, varMod){
	numObs = dim(loc)[1]

	# Convert nominal range to phi parameter
	nu = 0.5
	phi = sqrt(8*nu)/nomR

	# Generate realization
	D = as.matrix(dist(loc))
	Sig = varMod*exp(-phi*D)
	return(mvrnorm(mu = rep(0, numObs), Sigma = Sig))
}

# Generate realization from logistic spatial field
genRelLogistic = function(loc, nomR, varMod, size){
	numObs = dim(loc)[1]

	# Generate realization of latent field
	obsX = genRel(loc, nomR, varMod)
	
	# Probit link
	probs = pnorm(obsX)
	
	# Binomial data
	obs = list(log.obs = rbinom(rep(1, numObs), size, probs), spatial.obs = obsX)
}

inla.qinv.mod = function(Q, constr, reordering = inla.reorderings(), numCov = 0) 
{
    Q = inla.sparse.check(Q)
    if (is(Q, "dgTMatrix")) {
    	# Storage types
    	nrow = dim(Q)[1]
        ncol = dim(Q)[2]
        elems = nrow * ncol
        datatype = 0
        valuetype = INLA:::inla.ifelse(is.integer(Q), integer(), double())
        matrixtype = 0
        storagetype = 1

        # Header
        h = integer(8)
    	valuetp = INLA:::inla.ifelse(identical(valuetype, integer()), 0, 1)
    	h = c(version, elems, nrow, ncol, datatype, valuetp, matrixtype, storagetype)

    	# Write file
    	fp = file(filename, "wb")
	    writeBin(as.integer(length(h)), fp)
    	writeBin(as.integer(h), fp)
    	if (datatype == 0) {
        	if (identical(valuetype, integer())) {
	            writeBin(as.integer(as.vector(A)), fp)
        	}
        	else {
	            writeBin(as.double(as.vector(A)), fp)
    	    }
    	}
   	    close(fp)

   	    # Continue as normal
        qinv.file = filename
        remove = TRUE
    }
    else if (is.character(Q)) {
        qinv.file = Q
        remove = FALSE
    }
    else {
        stop("This chould not happen.")
    }
    constr.file = inla.tempfile()
    if (!missing(constr) && !is.null(constr)) {
        stopifnot(is.list(constr))
        A = as.matrix(constr$A)
        e = as.numeric(constr$e)
        stopifnot(ncol(A) == ncol(Q))
        stopifnot(nrow(A) == length(e))
        xx = matrix(c(nrow(A), c(A), c(e)), ncol = 1)
        inla.write.fmesher.file(xx, filename = constr.file)
    }
    if (is.list(reordering)) {
        reordering = reordering$name
    }
    reordering = match.arg(reordering)
    out.file = inla.tempfile()
    if (inla.os("linux") || inla.os("mac")) {
        s = system(paste(shQuote(inla.getOption("inla.call")), 
            "-s -m qinv", "-r", reordering, qinv.file, constr.file, 
            out.file), intern = TRUE)
    }
    else if (inla.os("windows")) {
        s = system(paste(shQuote(inla.getOption("inla.call")), 
            "-s -m qinv", "-r", reordering, qinv.file, constr.file, 
            out.file), intern = TRUE)
    }
    else {
        stop("\n\tNot supported architecture.")
    }
    Qinv = inla.read.fmesher.file(out.file)
    if (remove) {
        unlink(qinv.file)
    }
    unlink(out.file)
    unlink(constr.file)
    return(Qinv)
}
