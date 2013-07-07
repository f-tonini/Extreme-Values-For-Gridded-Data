#----------------------------------------------------------------------------------------------------------
# Name:         EV_RL.r
# Purpose:      This main function estimates return levels using the block maxima/minima approach
#				on gridded data. All external modules are called within this script.
# Author:       Francesco Tonini
# Email: 	    f_tonini@hotmail.com
# Created:      11/1/2009
# Copyright:    (c) 2009 by Francesco Tonini
# License:    	GNU General Public License (GPL)
# Software:     Tested successfully using R-2.14.0 64-bit version(http://www.r-project.org/)
#-------------------------------------------------------------------------------------------

##Define the main working directory
##Always use either / or \\ to specify the path
mainDir <- 'C:/Temp'
#mainDir <- '~/Desktop/Temp'  #if you are on a MacOSX environment

##Let's set the working directory
setwd(file.path(mainDir))

##Define the subdirectory(ies) where your output files will be saved
subDir <- 'Output'

##Create a physical copy of the subdirectory folder(s) where you save your output
##If the directory already exists it gives a warning, BUT we can suppress it using showWarnings = FALSE
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)

##Use an external source file w/ all modules (functions) used within this script. 
##Use FULL PATH if source file is not in the same folder w/ this script
source('myEVfunctions.r')

##Path to folders in which you want to save all your output files
workdir_Output <- file.path(mainDir, subDir)

##Load all required libraries
print('Loading required libraries...')
load.packages()

##Let's set the simulation parameters:

##Alamata.tot is the dataset we built for the Alamata woreda, Tigray (Ethiopia).
##This can be replaced by an entire region or country (the bigger the extent, the more the code will take to finish)
dataset <- Alamata.tot

##Number of pixels within the chosen dataset
npixels <- length(unique(dataset$id))

##Return periods chosen for the analysis (in months, since our time series is monthy minima)	
rperiods <- c(10,100,1000)
		
##Let's create a dataframe to store the parameters of the GEV distributions and their confidence intervals	
tab.parameters <- data.frame(matrix(NA,npixels,11))
colnames(tab.parameters) <- c('id.pix','gum.flag','mu','sig','shp','mu.inf','mu.sup','sig.inf','sig.sup','shp.inf','shp.sup')
	
##Let's create some labels for the return levels estimates and their upper and lower confidence intervals
RL_labels <- paste('rl.',rperiods,sep='')
RL.LOW_labels <- paste('rl.low.',rperiods,sep='')
RL.UPPER_labels <- paste('rl.up.',rperiods,sep='')
	
##Let's create a dataframe to store return levels estimates and their upper and lower confidence intervals
tab.rlevels <- data.frame(matrix(NA,npixels,(length(rperiods) + (length(rperiods) * 2) + 1)))
colnames(tab.rlevels) <- c("id.pix",RL_labels,RL.LOW_labels,RL.UPPER_labels)

##Create an empty vector in which we would save the pixel IDs
##whose GEV/Gumbel parameter estimates did not converge	
px.idWarning <- c()  

######################################################################
##Flagged values (in DN) by VITO for SPOT VGT S10
#251=missing, 252=cloud, 253=snow, 254=sea, 255=background
#Converted to ADVI values using: ADVI = -1.25 + 0.01 * DN
#1.26=missing, 1.27=cloud, 1.28=snow, 1.29=sea, 1.30=background
######################################################################

##Create a Tk window element as a container
root <- tktoplevel()
tktitle(root) <- 'Extreme Value Analysis'
##Define the label of a progress bar
lab <- tk2label(root)
##Define the length and create a progress bar
pb <- tk2progress(root, length = 300)
#Configure 
tkconfigure(pb, value=0, maximum = npixels-1)
##Pack all labels and progress bars inside the Tk container
tkpack(lab)
tkpack(pb)
##Refresh the Tk window
tcl('update')

##LOOP through each pixel, create its time series of maxima/minima, 
##estimate GEV parameters and return levels			
for (pix in 1:npixels){
	
	#update progress bar according to the pixel currently being processed
	tkconfigure(lab, text = paste('Pixel',pix,'Out of',npixels, '    (',round(pix/npixels*100),'% )'))
    tkconfigure(pb, value = pix - 1)
	tcl('update')
		
	##Slice the main dataset by extracting only historical records for the current pixel
	pix.dataset <- dataset[dataset$id == pix,]
	
	##This section was developed specifically for the case of having dekadal (1 dekad = 10days) 
	##time series and extracting monthly (1 month = 3 dekads) MINIMA.
	##NOTE: If you are interested in taking monthly MAXIMA just use FUN = max 
	##NOTE2: If you wanna take YEARLY minima/maxima then you should use:
	##etiq <- factor(pix.dataset$year) ...BUT make sure the time series is longer than at least 30 years
	##otherwise the GEV estimates will not be reliable, as well as the return level estimates.
		
	##Identify the 'dekadal' variable in the dataset
	Dekads <- pix.dataset$decads
	
	##Identify the 'Year' variable in the dataset
	Years <- pix.dataset$year
	
	##Identify the 'Value' variable in the dataset
	Values <- pix.dataset$val

	##Let's store the labels (as factors) of all existing dekad-year pairs (regardless of the dekad number)
	etiq <- factor( paste(substring(Dekads,1,3), Years, sep='') )
	
	##Let's extract the minima (use FUN = max if you want maxima) for each  unique level of the factor
	monthExtract <- tapply(Values, etiq, FUN = min)
	
	##Since values in 'monthExtract' have an alphabetical order relative to the month (Apr1998,Apr1999,Apr2000,etc.)
	##we need to sort them so that the series of values follows the real time line
	
	##If your are taking YEARLY minima/maxima skip (1-2-3) and use only (4). This, because
	##years already have a natural ascending order:
	
	#(1) by taking unique factor levels, we have a natural time line order
	UniqueEtiq <- unique(etiq)
	
	#(2) let's now match the correct ordered position of unique factors with our monthly values
	PosIdx <- match(as.character(UniqueEtiq), names(monthExtract))
	
	#(3) now that we got all the correct positions, let's put all our values in a numerical vector. This
	#is the time series to be used with extreme value models
	px.min <- as.vector(monthExtract[PosIdx])
	
	#(4) Uncomment the next line and use it to build your YEARLY maxima/minima time series
	#px.min <- as.vector(monthExtract)  
	
	#If the pixel is coded as 1.29=sea or 1.30=background we simply put 999 
	if ( any(px.min == 1.29) | any(px.min == 1.30) ) {
		tab.parameters[pix,] <- 999 
		tab.rlevels[pix,] <- 999
		next
	}
	
	#Let's now look at values 1.26=missing, 1.27=cloud, 1.28=snow. They will need to be removed from the
	#numerical vector, otherwise the extreme models will be biased. Since in our study area there are very few 
	#or no cloud values, snow, or missing, we can either remove them from the numerical vector, or 
	#replace them with a simple average between N previous/following values. 
	Flags <- c(1.26, 1.27, 1.28)
	
	#replace flagged values with NAs (if there are any)
	px.min[ px.min %in% Flags] <- NA
	
	#if there are flagged values, count them and see if their number exceeds 15% of the total time series
	if (any(is.na(px.min))) {
		NA.sum <- sum(is.na(px.min))
		#print a warning on screen to let the user know that model estimates may be biased
		if (NA.sum / length(px.min) * 100 >= 15) cat('\nATTENTION: More than 15% of values are flagged.\nModel estimates may be biased')
	}
	
	#override the original time series by removing any NA (flagged) values
	px.min <- px.min[!is.na(px.min)]
	
	#estimate GEV model (if you are working with maxima, DO NOT use the '-' sign 
	minimi.gev <- gev.fit(-px.min,show=F)
	class(minimi.gev) <- "gev.fit"
	
	#estimate Gumbel model (if you are working with maxima, DO NOT use the '-' sign 
	minimi.gum <- gum.fit(-px.min,show=F)
	class(minimi.gum) <- "gum.fit"
	
	#If the MLE estimation procedure does not converge, save the pixel ID into
	#the warning vector
	if (minimi.gev$conv > 0 | minimi.gum$conv > 0){
	
		px.idWarning <- c(px.idWarning, pix)
	
	}else{

		# Likelihood Ratio method to compare two models
		lik.ratio <- 2 * (minimi.gum$nllh - minimi.gev$nllh)
		pvalue <- (1 - pchisq(lik.ratio,1))
		
		#if the p-value is not significant (>.05) we use the gumbel distr. (shape parameter --> 0)
		if (pvalue > 0.05){
			
			#use a variable to flag the use of a gumbel distr.
			gum.flag <- 1
			mu <- minimi.gum$mle[1]  #location
			sig <- minimi.gum$mle[2] #scale
			shp <- 0                 #shape
			#conf. intervals (based ~ N(0,1))
			mu.inf <- minimi.gum$mle[1] - 1.96 * minimi.gum$se[1] 
			mu.sup <- minimi.gum$mle[1] + 1.96 * minimi.gum$se[1]
			sig.inf <- minimi.gum$mle[2] - 1.96 * minimi.gum$se[2]
			sig.sup <- minimi.gum$mle[2] + 1.96 * minimi.gum$se[2]
			shp.inf <- 0
			shp.sup <- 0
			
			#call external function to calculate return levels and their conf. intervals
			rl <- return.levels(minimi.gum, conf = 0.05, rperiods, make.plot = F)
			#let's store the return levels. If you are working with maxima, remove the '-' sign.
			rlevels <- rl$return.level
			#conf. interval calculated using the delta method. If you are working with maxima, remove the '-' sign.
			rlevels.int <- rl$confidence.delta
			#comment the following two lines if you are working with MAXIMA
			rlevels.int <- rlevels.int[,c(2:1)]
			colnames(rlevels.int) <- c('lower','upper')
				
		}else{ #otherwise use the GEV distr.
			
			#use a variable to flag the NON use of a gumbel distr.			
			gum.flag <- 0
			mu <- minimi.gev$mle[1]  #location
			sig <- minimi.gev$mle[2] #scale
			shp <- minimi.gev$mle[3] #shape
			#conf. intervals (based ~ N(0,1))
			mu.inf <- minimi.gev$mle[1] - 1.96 * minimi.gev$se[1]
			mu.sup <- minimi.gev$mle[1] + 1.96 * minimi.gev$se[1]
			sig.inf <- minimi.gev$mle[2] - 1.96 * minimi.gev$se[2]
			sig.sup <- minimi.gev$mle[2] + 1.96 * minimi.gev$se[2]
		    shp.inf <- minimi.gev$mle[3] - 1.96 * minimi.gev$se[3]
			shp.sup <- minimi.gev$mle[3] + 1.96 * minimi.gev$se[3]
			
			#call external function to calculate return levels and their conf. intervals			
			rl <- return.levels(minimi.gev, conf = 0.05, rperiods, make.plot = F)
			#let's store the return levels. If you are working with maxima, remove the '-' sign.
			rlevels <- -rl$return.level
			#conf. interval calculated using the delta method. If you are working with maxima, remove the '-' sign.
			rlevels.int <- -rl$confidence.delta
			#comment the following two lines if you are working with MAXIMA
			rlevels.int <- rlevels.int[,c(2:1)]
			colnames(rlevels.int) <- c('lower','upper')
				
		}

		#let's store the parameter and return level estimates in their respective dataset, for the row
		#corresponding to the current pixel in the LOOP
		tab.parameters[pix,] <- c(pix,gum.flag,round(mu,2),round(sig,2),round(shp,2),round(mu.inf,2),round(mu.sup,2),
								 round(sig.inf,2),round(sig.sup,2),round(shp.inf,2),round(shp.sup,2))

		tab.rlevels[pix,] <- matrix(c(pix, round(rlevels,2), round(rlevels.int,2)),1,(length(rperiods) + (length(rperiods)*2) + 1))
		
	}      
					
}

#For those pixels whose return levels could not be computed, use the median of the 4 nearest neighbors
if(length(px.idWarning) > 0) median.interp_RL()

#Save files to output folder
save_files_RL()
	
tkmessageBox(title = "Message", message = 'Done! Output Files Saved To Folder', icon = "info", type = "ok")
##Suppress the progress bar
tkdestroy(root)





