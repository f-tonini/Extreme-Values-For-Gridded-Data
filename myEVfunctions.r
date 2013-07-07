#----------------------------------------------------------------------------------------------------------
# Name:         myEVfunctions.r
# Purpose:      External functions called by the EV_RL.r & EV_RP.r scripts 
# Author:       Francesco Tonini
# Email: 	    f_tonini@hotmail.com
# Created:      11/1/2009
# Copyright:    (c) 2009 by Francesco Tonini
# License:    	GNU General Public License (GPL)
# Software:     Tested successfully using R-2.14.0 64-bit version(http://www.r-project.org/)
#-------------------------------------------------------------------------------------------

##CALL PACKAGES MODULE
##This module is used to call all required packages
load.packages <- function()
{	
	#setRepositories(ind=1:2)
	pkg <- c("ismev","stats","tcltk2","RANN","maptools","evd")
	w <- which(pkg %in% row.names(installed.packages()) == FALSE)
	if (length(w) > 0) install.packages(pkg)[w]
	#install.packages(c("ismev","stats","tcltk2","RANN","maptools","evd")) alternatively
	
	#Load (call) specific packages from the existing library collection into the current R session
	require(ismev)       #An Introduction to Statistical Modeling of Extreme Values
	require(stats)       #R statistical functions
	require(tcltk2)      #Tkinter objects and classes
	require(RANN)        #spatial functions (e.g. NN interpolation)
	require(maptools)    #reading and writing shapefiles
	require(evd)		 #Functions for extreme value distributions
	#require(raster)     #functions for raster data type
	#require(extRemes)   #only for automatic return-values calculation
	
	cat('\nAll libraries have been loaded successfully!\n')
}

##GRADIENT FUNCTION MODULE
##This module is used within  the return level module to compute the confidence intervals
##with delta method
gevrlgradient<-function (z, p)
{
    scale <- z$mle[2]
    shape <- z$mle[3]
    if (shape < 0)
        zero.p <- p == 0
    else zero.p <- logical(length(p))
    out <- matrix(NA, nrow = 3, ncol = length(p))
    out[1, ] <- 1
    if (any(zero.p)) {
        out[2, zero.p & !is.na(zero.p)] <- rep(-shape^(-1), sum(zero.p,
            na.rm = TRUE))
        out[3, zero.p & !is.na(zero.p)] <- rep(scale * (shape^(-2)),
            sum(zero.p, na.rm = TRUE))
    }
    if (any(!zero.p)) {
        yp <- -log(1 - p[!zero.p])
        out[2, !zero.p] <- -shape^(-1) * (1 - yp^(-shape))
        out[3, !zero.p] <- scale * (shape^(-2)) * (1 - yp^(-shape)) -
            scale * shape^(-1) * yp^(-shape) * log(yp)
    }
    return(out)
}

##RETURN LEVEL CALCULATION MODULE
##This module computes return levels based on the distibution used (GEV, gumbel, GPD)
return.levels <- function (z, conf = 0.05, rperiods = c(10, 100, 210, 510, 810, 980), make.plot = TRUE)
{
    out <- list()
    out$conf.level <- conf
    eps <- 1e-06
    a <- z$mle

    #prova
    #a1<-z$mle[1] + z$mle[2]*rperiods
   
    std <- z$se
    mat <- z$cov
    dat <- z$data
	
    kappa <- qnorm(conf/2, lower.tail = FALSE)
    nx <- length(rperiods)
    cl <- 1 - conf

    if (class(z) == "gev.fit") {
        
		if (is.null(rperiods)) rperiods <- seq(1.1, 1000, , 200)
        if (any(rperiods <= 1))
            stop("return.level: this function presently only supports return periods >= 1")
        yp <- -log(1 - 1/rperiods)
        if (a[3] < 0)
            zero.p <- yp == 0
        else zero.p <- logical(length(rperiods))
        zero.p[is.na(zero.p)] <- FALSE
        q <- numeric(length(rperiods))
        if (any(zero.p))
            q[zero.p] <- a[1] - a[2]/a[3]
        if (any(!zero.p)) {
            if (a[3] != 0)
                q[!zero.p] <- a[1] - (a[2]/a[3]) * (1 - (yp[!zero.p])^(-a[3]))
            else if (a[3] == 0)
                q[!zero.p] <- a[1] - a[2] * log(yp[!zero.p])
        }
		
        d <- gevrlgradient(z = z, p = 1/rperiods)
        v <- apply(d, 2, q.form, m = mat)
        yl <- c(min(dat, q, na.rm = TRUE), max(dat, q, na.rm = TRUE))
        if (make.plot) {
            xp <- 1/yp
            plot(xp, q, log = "x", type = "n", xlim = c(0.1,
                1000), ylim = yl, xlab = "Return Period", ylab = "Return Level",
                xaxt = "n")
            axis(1, at = c(0.1, 1, 10, 100, 1000), labels = c("0.1",
                "1", "10", "100", "1000"))
            lines(xp, q)
            lines(xp, (q + kappa * sqrt(v)), col = "blue")
            lines(xp, (q - kappa * sqrt(v)), col = "blue")
            points(-1/log((1:length(dat))/(length(dat) + 1)),
                sort(dat))
        }
        out$return.level <- q
        out$return.period <- rperiods
        conf3 <- cbind(q - kappa * sqrt(v), q + kappa * sqrt(v))
        colnames(conf3) <- c("lower", "upper")
        out$confidence.delta <- conf3
        
	}


	if (class(z) == "gum.fit") {

        if (is.null(rperiods))
            rperiods <- seq(1.1, 1000,  length=200)
        if (any(rperiods <= 1))
            stop("return.level: this function presently only supports return periods >= 1")
        yp <- -log(1 - 1/rperiods)
        q <- a[1] - a[2] * log(yp)
        vq <- std[1]^2 + ((-log(yp))^2 * std[2]^2)
	      sq <- sqrt(vq)

		yl <- c(min(dat, q, na.rm = TRUE), max(dat, q, na.rm = TRUE))
        if (make.plot) {
            xp <- 1/yp
            plot(xp, q, log = "x", type = "n", xlim = c(0.1,
                1000), ylim = yl, xlab = "Return Period", ylab = "Return Level",
                xaxt = "n")
            axis(1, at = c(0.1, 1, 10, 100, 1000), labels = c("0.1",
                "1", "10", "100", "1000"))
            lines(xp, q)
            lines(xp, (q + kappa * sq), col = "blue")
            lines(xp, (q - kappa * sq), col = "blue")
            points(-1/log((1:length(dat))/(length(dat) + 1)),
                sort(dat))
        }

		out$return.level <- q
        out$return.period <- rperiods
        conf3 <- cbind(q - kappa * sq, q + kappa * sq)
        colnames(conf3) <- c("lower", "upper")
        out$confidence.delta <- conf3
    }

	if (class(z) == "gpd.fit") {
	
        u <- z$threshold
        la <- z$rate
        a <- c(la, a)
        n <- z$n
        npy <- z$npy
        xdat <- z$xdata
        if (is.null(rperiods)) {
            rperiods <- seq(0.1, 1000, , 200)
        }
        m <- rperiods * npy
        if (a[3] == 0)
            q <- u + a[2] * log(m * la)
        else q <- u + a[2]/a[3] * ((m * la)^(a[3]) - 1)
        d <- gpdrlgradient(z, m)
        mat <- matrix(c((la * (1 - la))/n, 0, 0, 0, mat[1, 1],
            mat[1, 2], 0, mat[2, 1], mat[2, 2]), nc = 3)
        v <- apply(d, 2, q.form, m = mat)
        yl <- c(u, max(xdat, q[q > u - 1] + kappa * sqrt(v)[q >
            u - 1], na.rm = TRUE))
        if (make.plot) {
            if (any(is.na(yl)))
                yl <- range(q, na.rm = TRUE)
            plot(m/npy, q, log = "x", type = "n", xlim = c(0.1,
                max(m)/npy), ylim = yl, xlab = "Return period (years)",
                ylab = "Return level", xaxt = "n")
            axis(1, at = c(0.1, 1, 10, 100, 1000), labels = c("0.1",
                "1", "10", "100", "1000"))
            lines(m[q > u - 1]/npy, q[q > u - 1])
            lines((m[q > u - 1]/npy), (q[q > u - 1] + kappa *
                sqrt(v)[q > u - 1]), col = "blue")
            lines((m[q > u - 1]/npy), (q[q > u - 1] - kappa *
                sqrt(v)[q > u - 1]), col = "blue")
            nl <- n - length(dat) + 1
            sdat <- sort(xdat)
            points((1/(1 - (1:n)/(n + 1))/npy)[sdat > u], sdat[sdat >
                u])
        }
        out$return.level <- q
        out$return.period <- m/npy
        conf3 <- cbind(q[q > u - 1] - kappa * sqrt(v)[q > u -
            1], q[q > u - 1] + kappa * sqrt(v)[q > u - 1])
        colnames(conf3) <- c("lower", "upper")
        out$confidence.delta <- conf3
        
    }

invisible(out)

}

##NN INTERPOLATION FOR WARNINGS MODULE (FOR RETURN LEVELS)
##This module serves the purpose of replacing pixels, flagged as warning in the main script
##because of convergence issues, with the median of the nearest 4 or 8 neighbors

##NOTE: this works better if the raster resolution is not too coarse and the terrain is 
##quite uniform to that of the central pixel. This is not recommended if that is not the case, or 
##if you are using a very coarse resolution raster (e.g. 10 Km or more).
median.interp_RL <- function()
{
	
	Matr_coord <- cbind(dataset[1:npixels,]$x, dataset[1:npixels,]$y)
	
	for (px in px.idWarning){
	
		#Rook's neighborhood (using the 'spdep' package)
		NN4 <- knearneigh(Matr_coord,4)
		Matr_NN4 <- NN4$nn
		
		#Queen's neighborhood (using the 'spdep' package)
		#NN8 <- knearneigh(Matr_coord,8)
		#Matr_NN8 <- NN8$nn

		#take a subset of the NN matrix only for the pixels flagged with a warning
		Sub4 <- Matr_NN4[px,]
		#Sub8 <- Matr_NN8[px,]
		
		Mat_param <- tab.parameters[Sub4,]
		Mat_rlevels <- tab.rlevels[Sub4,]
		
		tab.parameters[px,3:length(tab.parameters)] <- round(sapply(Mat_param[,-c(1:2)],median),2)
		tab.rlevels[px,2:length(tab.rlevels)] <- round(sapply(Mat_rlevels[,-1],median),2)
	}

}

##NN INTERPOLATION FOR WARNINGS MODULE (FOR PROBABILITY EXCEEDANCES)
##This module serves the purpose of replacing pixels, flagged as warning in the main script
##because of convergence issues, with the median of the nearest 4 or 8 neighbors

##NOTE: this works better if the raster resolution is not too coarse and the terrain is 
##quite uniform to that of the central pixel. This is not recommended if that is not the case, or 
##if you are using a very coarse resolution raster (e.g. 10 Km or more).
median.interp_RP <- function()
{
	
	Matr_coord <- cbind(dataset[1:npixels,]$x, dataset[1:npixels,]$y)
	
	for (px in px.idWarning){
	
		#Rook's neighborhood (using the 'spdep' package)
		NN4 <- knearneigh(Matr_coord,4)
		Matr_NN4 <- NN4$nn
		
		#Queen's neighborhood (using the 'spdep' package)
		#NN8 <- knearneigh(Matr_coord,8)
		#Matr_NN8 <- NN8$nn

		#take a subset of the NN matrix only for the pixels flagged with a warning
		Sub4 <- Matr_NN4[px,]
		#Sub8 <- Matr_NN8[px,]
		
		Mat_param <- tab.parameters[Sub4,]
		Mat_prob <- tab.prob[Sub4,]
		
		tab.parameters[px,3:length(tab.parameters)] <- round(sapply(Mat_param[,-c(1:2)],median),2)
		tab.prob[px,2:length(tab.prob)] <- round(sapply(Mat_prob[,-1],median),2)
	}

}

##SAVE SHAPEFILE MODULE (FOR RETURN LEVELS)
##This module is used to save a point dataset to a shapefile (.shp)
save_files_RL <- function()
{
	
	Matr_coord <- cbind(dataset[1:npixels,]$x, dataset[1:npixels,]$y)
	
	tab.parameters$X <- Matr_coord[,1]
	tab.parameters$Y <- Matr_coord[,2]
	tab.rlevels$X <- Matr_coord[,1]
	tab.rlevels$Y <- Matr_coord[,2]
	
	##Turn input point dataset into a SpatialPointsDataFrame
	coordinates(tab.parameters) <- ~X+Y
	coordinates(tab.rlevels) <- ~X+Y
	
	##Path to folder
	path <- file.path(mainDir, subDir)
	
	##Write .csv tables
	write.table(tab.parameters, paste(path,'/MLE_Parameters.csv',sep=''), row.names=F, sep=',')
	write.table(tab.rlevels, paste(path,'/ReturnLevels.csv',sep=''), row.names=F, sep=',')
	
	##Write/Save the shapefile (library 'maptools')
	writeSpatialShape( tab.parameters, paste(path,'/MLE_Parameters',sep='') )
	writeSpatialShape( tab.rlevels, paste(path,'/ReturnLevels',sep='') )
	
	##Write/Save rasters (library 'raster')
	
}

##SAVE SHAPEFILE MODULE (FOR PROBABILITY EXCEEDANCES)
##This module is used to save a point dataset to a shapefile (.shp)
save_files_RP <- function()
{
	
	Matr_coord <- cbind(dataset[1:npixels,]$x, dataset[1:npixels,]$y)
	
	tab.parameters$X <- Matr_coord[,1]
	tab.parameters$Y <- Matr_coord[,2]
	tab.prob$X <- Matr_coord[,1]
	tab.prob$Y <- Matr_coord[,2]
	
	##Turn input point dataset into a SpatialPointsDataFrame
	coordinates(tab.parameters) <- ~X+Y
	coordinates(tab.prob) <- ~X+Y
	
	##Path to folder
	path <- file.path(mainDir, subDir)
	
	##Write .csv tables
	write.table(tab.parameters, paste(path,'/MLE_Parameters.csv',sep=''), row.names=F, sep=',')
	write.table(tab.prob, paste(path,'/ProbExceedance.csv',sep=''), row.names=F, sep=',')
	
	##Write/Save the shapefile (library 'maptools')
	writeSpatialShape( tab.parameters, paste(path,'/MLE_Parameters',sep='') )
	writeSpatialShape( tab.prob, paste(path,'/ProbExceedance',sep='') )
	
	##Write/Save rasters (library 'raster')
	
}


##PROBABILITY OF EXCEEDANCE MODULE
##This module is used to calculate the probability that the event will be exceeded in any one month
##(or year, or season, depending on what is your maxima/minima time series unit). 
##The return period is the inverse of the probability of exceedance
prob_fun <- function(dataset,value,tab.extreme,nome,n=3){
	
	px.prob <- round(1 - pgev(-return_levels, loc=mu, scale=sig, shape=shp), 3)
	
	#Use the following line, and not the previous, if you are working with MAXIMA instead of MINIMA
	#px.prob <- round(1 - pgev(-return_levels, loc=mu, scale=sig, shape=shp), 3)
	
	#NOTE: the return period will then simply be 1/prob (in the time unit used, e.g. months)
	
	return(px.prob)

}




	
