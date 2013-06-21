
### DEFINE ALL FUNCTIONS TO BE USED IN THE MAIN ROUTINE ###

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


NN.interp <- function(input){
	
	out <- list()
	
	print('MLE did not converge in one or more cells: those values are being interpolated...')
	for (i in input){
		
		NN_fact <- 2
		#create a circular buffer around each NA cell to include 8 Near. Neighbors
		Buffer <- disc(res(layer_stack[[1]])[1] * NN_fact, centre=coords[i,])
		BufferPoly = data.frame(x = Buffer$bdry[[1]]$x, y = Buffer$bdry[[1]]$y)
		idx <- point.in.polygon(point.x = coords[,1], point.y = coords[,2], pol.x = BufferPoly[,1], pol.y = BufferPoly[,2])	
			
		while(!any(idx == 1)){
			NN_fact <- NN_fact + 1
			Buffer <- disc(res(layer_stack[[1]])[1] * NN_fact, centre=coords[i,])
			BufferPoly = data.frame(x = Buffer$bdry[[1]]$x, y = Buffer$bdry[[1]]$y)
			idx <- point.in.polygon(point.x = coords[,1], point.y = coords[,2], pol.x = BufferPoly[,1], pol.y = BufferPoly[,2])	
		}
				
		tab.rlevels.NN <- tab.rlevels[idx == 1,]
		tab.rlevels.NN <- tab.rlevels.NN[tab.rlevels.NN$id_pix!=i,]
		tab.rlevels.NN <- tab.rlevels.NN[,-1]
				
		#replace NA cell value with median of 8 NN
		if(useConfInt == 1 | length(rperiods) > 1) {
			newVal.rlevels <- apply(tab.rlevels.NN,2,FUN=median,na.rm=T) 
		}else{ 
			newVal.rlevels <- median(tab.rlevels.NN[!is.na(tab.rlevels.NN)])
		}
		tab.rlevels[i,2:ncol(tab.rlevels)] <- newVal.rlevels
		
		if(useCoeff == 1){

			tab.parameters.NN <- tab.parameters[idx == 1,]
			tab.parameters.NN <- tab.parameters.NN[tab.parameters.NN$id_pix!=i,]
			tab.parameters.NN <- tab.parameters.NN[,-c(1,2)]
				
			#replace NA cell value with median of 8 NN
			newVal.parameters <- apply(tab.parameters.NN,2,FUN=median,na.rm=T)
			tab.parameters[i,3:ncol(tab.parameters)] <- newVal.parameters

		}
	}
	
	out$tab.rlevels <- tab.rlevels
	if(useCoeff == 1) out$tab.parameters <- tab.parameters
	return(out)
}

save.raster <- function(input,prj){
	
	if(is.na(prj)) print('No projection defined for input rasters!')
	
	r <- layer_stack[[1]]
	
	if(useConfInt == 1 | length(rperiods) > 1){
	
		for (i in 1:ncol(input)){

			rlevel.raster <- rasterize(coords,r,field=input[,i])
			proj4string(rlevel.raster) <- prj
			writeRaster(rlevel.raster, filename=paste(outputWorkSpace,'\\',names(input)[i],sep=''), format='HFA', datatype='FLT4S', overwrite=TRUE)

		}
	
	}else{
	
		rlevel.raster <- rasterize(coords,r,field=input)
		proj4string(rlevel.raster) <- prj
		writeRaster(rlevel.raster, filename=paste(outputWorkSpace,'\\rl_',rperiods,sep=''), format='HFA', datatype='FLT4S', overwrite=TRUE)
	
	}
	
}

save.shp <- function(input,prj,fName){

	if(is.na(prj)) print('No projection defined for input rasters!')
	
	input$X <- coords[,1]
	input$Y <- coords[,2]

	coordinates(input) <- ~X+Y	
	proj4string(input) <- prj
	
	##Write .shp file
	writeOGR(input, outputWorkSpace, layer = fName, driver='ESRI Shapefile',overwrite_layer=TRUE)

}
