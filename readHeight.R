# Need INLA library
require(INLA)

# Get height at UTM33 coordinates
getHeight = function(loc){
    ## Read elevation from disk
        # Read lat/long gridded data from disk
        nRow = 4800
        nCol = 10800
        numEl = nRow*nCol
        
        # Check if file should be downloaded
        ex = file.exists("data/c10g")
        if(!ex)
            download.file("https://dl.dropboxusercontent.com/u/8401189/c10g", "data/c10g")
        
        # Read file
        inFile = file("data/c10g", "rb")
        hmm = readBin(inFile, integer(), n = numEl, size = 2)
        close(inFile)
        
        # Re-size
        height = matrix(hmm, ncol = nRow)
        height = t(height)
        
        # Properties of gridded data
        xStart = 0.004166666667
        yStart = 89.995833333333
        xStep = 0.008333333333
        yStep = 0.008333333333
        
        # Lat and lon values
        lat = seq(yStart, by = -yStep, length.out = nRow)
        lon = seq(xStart, by = xStep, length.out = nCol)
        
        # Grid of lat and long
        LAT = matrix(rep(lat, length(lon)), ncol = nCol)
        LON = matrix(rep(lon, each = length(lat)), ncol = nCol)
    
    ## Calculate gradients
        # Calculate gradients on spherical earth
        ds = xStep*(pi/180)
        nRow = dim(height)[1]
        nCol = dim(height)[2]
        Rearth = 6378
        
        # Longitude
        lonGrad = (height[2:(nRow-1), 3:nCol] - height[2:(nRow-1), 1:(nCol-2)])/(2*ds)
        lonGrad = lonGrad/(Rearth*sin(LAT[2:(nRow-1),2:(nCol-1)]*pi/180))
        
        # Latitude
        latGrad = (height[3:nRow, 2:(nCol-1)] - height[1:(nRow-2), 2:(nCol-1)])/(2*ds)
        latGrad = latGrad/Rearth
        
        # Combine
        gradMatrix = matrix(NA, nrow = nRow, ncol = nCol)
        gradMatrix[2:(nRow-1), 2:(nCol-1)] = sqrt(lonGrad^2+latGrad^2)
        
    ## Choose closest point
        # Convert input to lat/lon
        require(PBSmapping)
        xydata = data.frame(X = loc[,1], Y = loc[,2])
        attr(xydata, "projection") = "UTM"
        attr(xydata, "zone") = 33
        loc = convUL(xydata)
        
        hCov = loc[,1]
        gCov = hCov
        lonIdx = round((loc[, 1] - lon[1])/diff(range(lon))*nCol)+1
        latIdx = round((lat[1]-loc[, 2])/diff(range(lat))*nRow)+1
        
        for(idx in 1:dim(loc)[1]){
            if(lonIdx[idx] < 1){
#                cat('Lon outside range:', loc[idx,1], '\n')
                hCov[idx] = NA
                gCov[idx] = NA
            } else{
                hCov[idx] = height[latIdx[idx], lonIdx[idx]]
                gCov[idx] = gradMatrix[latIdx[idx], lonIdx[idx]]
            }
        }
        
        #levelplot(z~xs+ys, panel = panel.levelplot.points, col.regions = heat.colors(50))
    return(list(hCov = hCov, gCov = gCov))
}

getSmooth = function(loc, grid, mat){
    # Store output
    z = loc[,1]*NA
    
    # Get indicies into matrix
    xIdx = round((loc[, 1] - grid$x$start)/grid$x$ds)+1
    yIdx = round((loc[, 2] - grid$y$start)/grid$y$ds)+1
    
    for(idx in 1:dim(loc)[1]){
        z[idx] = mat[yIdx[idx], xIdx[idx]]
    }
    
    return(z)
}

getCovariates = function(mesh){
    # Grid size
    ds = 0.5
    easting = seq(-400, 800, by = ds)
    northing = seq(6100, 7500, by = ds)
    E = rep(easting, each = length(northing))
    N = rep(northing, length(easting))
    
    # Get data
    loc = cbind(E, N)
    covar = getHeight(loc)
    height = covar$hCov
    grad = covar$gCov
    
    # Replace sea with 0 and outside region with 0
    height[height<0] = 0
    grad[grad < 0] = 0
    height[is.na(height)] = 0
    grad[is.na(grad)] = 0
    
    # Create A matrix for gridded data
    A.proj = inla.spde.make.A(mesh = mesh, loc = loc)
    
    # Get approximate projection weights
    intHeight = t(A.proj)%*%height
    intUnit = t(A.proj)%*%rep(1, length(height))
    hCov = intHeight/intUnit 
    
	# Project to mesh
	intGrad = t(A.proj)%*%grad
	gCov = intGrad/intUnit
	
	
	return(list(covHeight = hCov, covGrad = gCov))
}
