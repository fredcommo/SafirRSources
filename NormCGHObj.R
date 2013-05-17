# Normalize CGH

#################
# Cy5/Cy3 adjustment and Log2Ratio calculation
CyAdjust <- function(Data, Fract = 1){
#	cat("Normalisation Cy5/Cy3\n")
	g <- log2(Data$gMedianSignal)				# Ref in Cy3
	r <- log2(Data$rMedianSignal)				# Test in Cy5

	if(!is.null(Fract)){							# Calculates weights to correct the dilution effect (due to tumor cell rate). No effect if Fract = 1. DO NOT USE until validation !
		Q <- quantile(r, Fract, na.rm = TRUE)
		w <- 1/(1+exp(1/sqrt(Fract)*(Q-r)))
		r <- r*(1 + (1-Fract)*w)
		}
	M <- r - g
	A <- (r + g)/2
	Loess <- loessFit(M, A)$fitted
	LR <- M - Loess
	Data$Log2Ratio <- LR - median(LR, na.rm = T)
	return (Data)
	}
#################


#################
# GC% adjustment
GCadjust <- function(Data){
	lr = Data$Log2Ratio
	GC <- Data$GCpercent
	adjLr <- lr - loessFit(lr, GC)$fitted
	Data$Log2Ratio = adjLr
	return(Data = Data)
	}
#################

#################
mergePeaks <- function(nG, m, s, p, MergeVal){
		Raw <- c(1, 1)
		while(length(Raw)!=0){
			Mdist <- matrix(0, nG, nG)
			for(i in 1:nG)
				for(j in 1:nG){
				Mdist[i,j] <- abs(m[i] - m[j])
				}
			diag(Mdist) <- NA
			Raw <- ceiling(which(Mdist<MergeVal)/nG)
			cat(Raw, "\n")
			if(length(Raw)!=0){
				C1 <- Raw[1]
				C2 <- Raw[2]
				pp = p[C1] * p[C2]
				m[C1] <- ((p[C1]-pp)*m[C1] + (p[C2]-pp)*m[C2])/(p[C1] + p[C2] - 2*pp)
				s[C1] <- ((p[C1]-pp)*s[C1] + (p[C1]-pp)*s[C2])/(p[C1] + p[C2] - 2*pp)
				p[C1] <- p[C1] + p[C2] - 2*pp
				m <- m[-C2]; s <- s[-C2]; p <- p[-C2]
				nG <- length(m)
				cat("means:", m, ", props:", p, "\n")
				}
			}
		retunr(list(nG = nG, m = m, s = s, p = p))
		}
#################

#################
EMcentr <- function(Data, useX11, Save){

	testLr <- Data$Log2Ratio
	if(getInfo(Data, 'platform') == 'Affymetrix'){
		testLr <- testLr[seq(1, length(testLr), by = 6)]
		}

	K = 51
	runLr <- runmed(testLr, k = K)	
	q1 <- cut[1]
	q2 <- cut[2]
	suppr <- which(runLr<=q1 | runLr>=q2)
	runLr2 <- runLr
	if(length(suppr)>0) runLr2 <- runLr[-suppr]
	len <- floor(length(runLr2)*0.25)
	
	model <- Mclust(runLr2[sample(1:length(runLr2), len)], G = G)
	nG <- model$G
	p <- model$parameters$pro
	m <- model$parameters$mean
	s <- model$parameters$variance$sigmasq
	if(length(s)<length(m)) s <- rep(s, length(m))
	p <- p[order(m)]
	s <- s[order(m)]
	m <- sort(m)

	# merge Classes: depending on MergePeaks and MergeVal
	#
	if(MergPeak){
		values <- mergePeaks(nG, m, s, p, MergeVal)
		nG <- values$nG
		m <- values$m
		s <- values$s
		p <- values$p
	}

	# compute densities
	denLr <- density(runLr)
	max.dlr <- max(denLr$y)
	maxy <- max(dlr$y)*1.25
	densList = list()
	dy <- cdy <- c()
	n <- length(lr)
	if(length(suppr)>0) n <- length(lr[-suppr])
	for(i in 1:length(m)){
		tmp <- rnorm(n*p[i], m[i], sqrt(s[i])) # rnorm(n*p[i], m[i], sqrt(s[i]))
		densList[[i]] <- tmpDens <- density(tmp, na.rm = T)
		dy <- c(dy, max(tmpDens$y*p[i]))
		cdy <- c(cdy, sum(tmpDens$y*p[i]))
		}

	if(Plot){
		if (useX11)
			X11(type="dbcairo")
			plot(dlr$x, dlr$y, type = "n", xlab = "Log2R", ylab = "Density", xlim = range(q1*Expand, q2*Expand), ylim = range(0, max.dlr*1.25))
			polygon(dlr$x, dlr$y, col = 'grey90')
			for (i in 1:length(densList)){
				L = densList[[i]]
				lines(L$x, I(L$y*p[i]), lwd = 1, col = rgb(i/nG, 0.2, (nG-i)/nG, 0.75))
				polygon(L$x, I(L$y*p[i]), col = rgb(i/nG, 0.2, (nG-i)/nG, 0.25))
				text(x = m[i], y = min(dy[i]+ max(dlr$y)/8, maxy), labels = round(m[i], 3),
						cex = ifelse(i == max(peakIndex), 1.5, 1.25),
						font =  ifelse(i == max(peakIndex), 2, 1))
				title(main = paste(getInfo(Data, 'sampleId'), "\nEM centralization: correction.factor =", round(correct, 5)))
				}
		if (Save){
			Path = getwd()
			fileName = paste(getInfo(Data, 'sampleId'), '_', getInfo(Data, 'barCode'), '.png', sep = '')
			savePlot(fileName, type = 'png')
			}
		}	

	peakIndex <- which.min(abs(m))
	if (method == "Left") peakIndex <- which(dy>=max(dy)*peakThresh & m<=0)
	if (method == "Zero") {peakIndex <- which(dy>=max(dy)*peakThresh); peakIndex <- peakIndex[which.min(abs(m)[peakIndex])]}
	if (method == "Right") peakIndex <- which(dy>=max(dy)*peakThresh & m>=m[which.max(dy)])
	peakIndex = peakIndex[which.max(m[peakIndex])]
	correct <- m[max(peakIndex)]

	cat("n.peaks = ", nG, "peakIndex = ", peakIndex, "\n") 
	cat("Valeur de centrage EM: correction factor = ", correct, "\n")

	# return the full cgh objet containing the Log2Ratios adjusted for Cy and CG, and centralized by EM
	Data$Log2Ratio = Data$Log2Ratio - correct
	return (Data)
	}	
#################





NormCGHObj <- function(cghObj, Cy = TRUE, GC = TRUE, Centr = TRUE, Fract = 1, span = 0.3,
										cut = c(-0.75, 0.75), G = 3:7, peakThresh = 0.75, MergePeaks = T, MergeVal = 0.075, method = "Left",
										Plot = TRUE, Expand = 1.75, useX11 = FALSE, Save = TRUE, Root = Root){

# cgh:				CGH data from loadAgilent or Load Affy
# Cy, CG, Centr:	Indicates if theses corrections have to be applied. If is.Agilent = FALSE, Cy = CG = FALSE
# Fract:			Indicates the approximative proportion of tumor cells (do not use before the method has been validated) 
# span:				Unused argument
# cut:				Quantiles thresholds in EM centralization
# G:				Number of groups to consider in EM centralization
# peakThresh: 		Proportion of the maximum density value to consider a peak as the central reference
# MergePeaks:		Allow to merge two peaks if there distance is lower than Mergeval
# MergeVal:			Minimum distance to consider two peaks as different.
# method:			Define the method to choose the reference peak. Left = left major peak, Zero = major peak closed to Zero, Right = right major peak.
# Plot:	 			Edit a density plot with EM peaks
# Expand:			A graphic parameter to expand the x axis
# Save:				To save automatically the density plot. The current default folder is Root/~/Safir CentralizeProfiles
# Root: 			The root to save the density plot.

	Data <- getCNset(cghObj)
	platform = getInfo(cghObj, 'platform')
	
	if (platform == 'Affymetrix') Cy = FALSE
	else{
		if(!Cy){
			g <- log2(Data$gMedianSignal)				# Ref in Cy3
			r <- log2(Data$rMedianSignal)				# Test in Cy5
			pos <- Data$genomicPos
			lr <- r - g
			Data$Log2Ratio <- lr - median(lr)
			}
		else
			Data <- CyAdjust(Data)
		}
		
	
	if(GC){
		adjust <- GCadjust(Data)
		Data <- adjust$Data
		# s2 <- adjust$s2
		# sumSpread <- adjust$sumSpread
		# medianSpread <- adjust$medianSpread
		}

	if(EM){
		Data <- EMcentr(Data)
		}
	cghObj@cgh = Data
	
	# To define in object class
	cghObj@param = c(CyAdjust = Cy, GCAdjust = GC, EMcentralized = EM, nPeak = nG, Means = m, Props = p, sDev = sqrt(S))
}
	