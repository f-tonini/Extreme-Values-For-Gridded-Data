#### Load Libraries ####
print("Loading Libraries....")

library(ismev)       #An Introduction to Statistical Modeling of Extreme Values
library(raster)
library(rgdal)
library(spatstat)
library(SDMTools) 

#### Get Arguments ####
Args <- commandArgs()
RSource <- Args[5]
inputWorkSpace <- Args[6]
outputWorkSpace <- Args[7]
rperiods <- Args[8]
useMinima <- as.integer(Args[9])
useConfInt <- as.integer(Args[10])
useCoeff <- as.integer(Args[11])

rperiods = strsplit(rperiods, ";")
rperiods = as.integer(c(unlist(rperiods)))

print(paste(commandArgs(), collapse=" "))

#Play with this command from the 'raster' package to avoid memory issues
rasterOptions(chunksize=1e+06, maxmemory=1e+07)

##Source to an external .R script to call a set of functions in the main routine (here) 
source(paste(dirname(RSource),'myEVfunctions.r',sep='/'))

### BEGINNING OF ROUTINE ###

files <- list.files(inputWorkSpace, full.names=F)
if(length(files) == 0) stop('Specified folder is empty!')

#let's assume for now it can only accept .img raster files in input
fileLst <- list.files(inputWorkSpace, pattern='\\.img$',full.names=T)

nbrs = sort(as.numeric(unique(sub("^([^.]*).*", "\\1", basename(fileLst))))) 

print('Reading and stacking rasters from input folder...')
layer_stack <- stack(paste(inputWorkSpace,'\\',nbrs,'.img',sep=''))

coords <- rasterToPoints(layer_stack[[1]])[,1:2]
X.data <- extract(layer_stack, coords)

##Number of pixels within the chosen dataset
npixels <- nrow(X.data)

if(useCoeff==1){
	tab.parameters <- data.frame(matrix(NA,npixels,11))
	colnames(tab.parameters) <- c('id_pix','gum_flag','mu','sig','shp','mu_inf','mu_sup','sig_inf','sig_sup','shp_inf','shp_sup')
}

print("Calculation of return levels begins...")

##Let's create some labels for the return levels estimates and their upper and lower confidence intervals
RL_labels <- paste('rl_',rperiods,sep='')

if(useConfInt == 1){
	RL.LOW_labels <- paste('rl_low_',rperiods,sep='')
	RL.UPPER_labels <- paste('rl_up_',rperiods,sep='')
	##Let's create a dataframe to store return levels estimates and their upper and lower confidence intervals
	tab.rlevels <- data.frame(matrix(NA,npixels,(length(rperiods) + (length(rperiods) * 2) + 1)))
	colnames(tab.rlevels) <- c("id_pix",RL_labels,RL.LOW_labels,RL.UPPER_labels)
}else{
	tab.rlevels <- data.frame(matrix(NA,npixels,length(rperiods) + 1))
	colnames(tab.rlevels) <- c("id_pix",RL_labels)
}

##Create an empty vector in which we would save the pixel IDs
##whose GEV/Gumbel parameter estimates did not converge	
px.idWarning <- c()  

##LOOP through each pixel, create its time series of maxima/minima, 
##estimate GEV parameters and return levels			

for (pix in 1:nrow(X.data)){

	if ( any(is.na(X.data[pix,])) ){
		
		w <- which(is.na(X.data[pix,]))
		cat('WARNING: one or more NA values in the time series. Values will be interpolated!\n')
		sapply(w, FUN=function(x){if(x==1){X.data[pix,x] <- X.data[pix,x+1]}else if(x==ncol(X.data)){X.data[pix,x] <- X.data[pix,x-1]}else{X.data[pix,x] <- (X.data[pix,x-1] + X.data[pix,x+1]) / 2}})
		
	}
	
	#estimate GEV model (if you are working with maxima, DO NOT use the '-' sign 
	if(useMinima==0){
		
		gev.est <- gev.fit(X.data[pix,],show=F)
		class(gev.est) <- "gev.fit"
			
		gum.est <- gum.fit(X.data[pix,],show=F)
		class(gum.est) <- "gum.fit"
		
	}else{
		
		gev.est <- gev.fit(-X.data[pix,],show=F)
		class(gev.est) <- "gev.fit"
			
		gum.est <- gum.fit(-X.data[pix,],show=F)
		class(gum.est) <- "gum.fit"
			
	}
		
	if (gev.est$conv > 0 | gum.est$conv > 0){
		
		#cat(paste('\nMLE for pixel #',pix,'did not converge!\n'))
		px.idWarning <- c(px.idWarning, pix)
		if(useCoeff==1) tab.parameters[pix,1] <- pix
		tab.rlevels[pix,1] <- pix
		next
	
	}
	
	# Likelihood Ratio method to compare two models
	lik.ratio <- 2 * (gum.est$nllh - gev.est$nllh)
	pvalue <- (1 - pchisq(lik.ratio,1))
		
	#if the p-value is not significant (>.05) we use the gumbel distr. (shape parameter --> 0)
	if (pvalue > 0.05){
			
		#use a variable to flag the use of a gumbel distr.
		gum.flag <- 1
		mu <- gum.est$mle[1]  #location
		sig <- gum.est$mle[2] #scale
		shp <- 0                 #shape
		#conf. intervals (based ~ N(0,1))
		mu.inf <- gum.est$mle[1] - 1.96 * gum.est$se[1] 
		mu.sup <- gum.est$mle[1] + 1.96 * gum.est$se[1]
		sig.inf <- gum.est$mle[2] - 1.96 * gum.est$se[2]
		sig.sup <- gum.est$mle[2] + 1.96 * gum.est$se[2]
		shp.inf <- 0
		shp.sup <- 0
			
		#call external function to calculate return levels and their conf. intervals
		rl <- return.levels(gum.est, conf = 0.05, rperiods, make.plot = F)
		#let's store the return levels. If you are working with maxima, remove the '-' sign.
		if(useMinima==1){
			rlevels <- -rl$return.level
			#conf. interval calculated using the delta method. If you are working with maxima, remove the '-' sign.
			if(useConfInt == 1){
				rlevels.int <- -rl$confidence.delta
				rlevels.int <- rlevels.int[,c(2:1)]
				colnames(rlevels.int) <- c('lower','upper')
			}
		}else{
			rlevels <- rl$return.level
			if(useConfInt == 1) rlevels.int <- rl$confidence.delta
		}
			
	}else{ #otherwise use the GEV distr.
			
		#use a variable to flag the NON use of a gumbel distr.			
		gum.flag <- 0
		mu <- gev.est$mle[1]  #location
		sig <- gev.est$mle[2] #scale
		shp <- gev.est$mle[3] #shape
		#conf. intervals (based ~ N(0,1))
		mu.inf <- gev.est$mle[1] - 1.96 * gev.est$se[1]
		mu.sup <- gev.est$mle[1] + 1.96 * gev.est$se[1]
		sig.inf <- gev.est$mle[2] - 1.96 * gev.est$se[2]
		sig.sup <- gev.est$mle[2] + 1.96 * gev.est$se[2]
		shp.inf <- gev.est$mle[3] - 1.96 * gev.est$se[3]
		shp.sup <- gev.est$mle[3] + 1.96 * gev.est$se[3]
			
		#call external function to calculate return levels and their conf. intervals			
		rl <- return.levels(gev.est, conf = 0.05, rperiods, make.plot = F)
		#let's store the return levels. If you are working with maxima, remove the '-' sign.
		if(useMinima==1){
			rlevels <- -rl$return.level
			#conf. interval calculated using the delta method. If you are working with maxima, remove the '-' sign.
			if(useConfInt == 1){
				rlevels.int <- -rl$confidence.delta
				rlevels.int <- rlevels.int[,c(2:1)]
				colnames(rlevels.int) <- c('lower','upper')
			}
		}else{
			rlevels <- rl$return.level
			if(useConfInt == 1) rlevels.int <- rl$confidence.delta
		}
				
	}

	#let's store the parameter and return level estimates in their respective dataset, for the row
	#corresponding to the current pixel in the LOOP
	if(useCoeff == 1) {
		tab.parameters[pix,] <- c(pix,gum.flag,round(mu,2),round(sig,2),round(shp,2),round(mu.inf,2),round(mu.sup,2),
							round(sig.inf,2),round(sig.sup,2),round(shp.inf,2),round(shp.sup,2))
	}
		
	
	if(useConfInt == 1){
		tab.rlevels[pix,] <- matrix(c(pix, round(rlevels,2), round(rlevels.int,2)),1,(length(rperiods) + (length(rperiods)*2) + 1))
	}else{
		tab.rlevels[pix,] <- matrix(c(pix, round(rlevels,2)),1,length(rperiods) + 1)
	}
		
}

if(length(px.idWarning)>0) out <- NN.interp(px.idWarning)

tab.rlevels <- out$tab.rlevels
tab.rlevels <- tab.rlevels[,-1]
save.raster(input=tab.rlevels,prj=proj4string(layer_stack[[1]]))
#save.shp(input=tab.rlevels,prj=proj4string(layer_stack[[1]]),fName="Return_Levels")

if(useCoeff == 1) {

	tab.parameters <- out$tab.parameters
	tab.parameters <- tab.parameters[,-1]
	save.shp(input=tab.parameters,prj=proj4string(layer_stack[[1]]),fName="MLE_coeff")
}

print("Calculations Complete...")