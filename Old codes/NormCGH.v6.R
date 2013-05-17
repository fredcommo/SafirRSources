NormCGH.v6 <- function(cgh, is.Agilent = TRUE, Cy = TRUE, GC = TRUE, Centr = TRUE, Fract = 1, span = 0.3, cut = c(-0.5, 0.5), G = 3:7, peakThresh = 0.75, MergePeaks = T, MergeVal = 0.075, method = "Left", Plot = TRUE, Expand = 1.75, Save = TRUE, Root = "E"){

# cgh:				CGH data from loadAgilent or Load Affy
# is.Agilent:		Type of platform. If FALSE (Affy) Cy correction and GC adjustement will not be used.
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

if(!Plot) Save = FALSE

Data <- cgh$cgh
PathPlot <- paste(Root, ":/Projet Safir/Data SAFIR/Safir Output/Safir CentralizeProfiles/", sep = "")

# For PreSafir only
# PathPlot <- "D:/Projet Safir/Data_PRESAFIR1/PreSafir_Output/PS_CentralizeProfiles/"

if(is.Agilent){
	g <- log2(Data$gMedianSignal)				# Ref in Cy3
	r <- log2(Data$rMedianSignal)				# Test in Cy5
	pos <- Data$GenomicPos
	lr <- r - g
	}
else {lr <- cgh$cgh$Log2Ratio
	Cy = FALSE
	GC = FALSE
	}

############################################
	# Cy5/Cy3 adjustment and Log2Ratio calculation
cgh$Cy.adjust <- "No"
if(Cy){
	cat("Normalisation Cy5/Cy3\n")
	if(!is.null(Fract)){							# Calculates weights to correct the dilution effect (due to tumor cell rate). No effect if Fract = 1. DO NOT USE until validation !
		Q <- quantile(r, Fract, na.rm = TRUE)
		w <- 1/(1+exp(1/sqrt(Fract)*(Q-r)))
		r <- r*(1 + (1-Fract)*w)
		}
	M <- r - g
	A <- (r + g)/2
	L <- loessFit(M, A)$fitted
	lr <- M - L
	lr <- lr - median(lr, na.rm = T)
	cgh$Cy.adjust <- "Yes"
	}

############################################
	# Ajustement GC
cgh$GC.Adjust <- "No"
if(GC){
	cat("Ajustement GC: ")
	K = 2*ceiling(length(lr)/800) + 1
	gc.lr <- lr
	Sspread <- Mspread <- S2 <- c()

	rmed <- runmed(lr, k = K)
	spread <- NA
	spread <- abs(diff(rmed))
	s2 <- sum((rmed - lr)^2)
	Sspread <- c(Sspread, sum(spread))
	Mspread <- c(Mspread, median(spread))
	S2 <- c(S2, s2)

	current <- Inf
	gc.crit <- NA
	is.GC <- which(substr(names(Data), 1, 2) == "GC")
	GC <- Data[,is.GC]

	for(i in 1:ncol(GC)){
		spread <- NA
		tmp <- lr - loessFit(lr, GC[,i])$fitted
		rmed <- runmed(tmp, k = K)
		spread <- abs(diff(rmed))
		s2 <- sum((rmed - tmp)^2)
		if(sum(spread) < current){
			current <- sum(spread)
			gc.lr <- tmp			
			gc.crit <- colnames(GC)[i]
			}
		Sspread <- c(Sspread, sum(spread))
		Mspread <- c(Mspread, median(spread))
		S2 <- c(S2, s2)
		}
	best.spread <- current
	GC.scores <- cbind.data.frame(c("none", colnames(GC)), Mspread, Sspread, S2)
	cat(gc.crit, "\n")
	lr <- gc.lr
	cgh$GC.Adjust <- gc.crit
	}


############################################
	# EM centralization
cgh$EM.centralization <- "No"
if(Centr){
	cat("Centrage EM\n")
	lr <- lr - median(lr, na.rm = T)
	test.lr <- lr
	if(!is.Agilent){
		test.lr <- lr[seq(1, length(lr), by = 6)]
		}
	# K = 2*ceiling(length(test.lr)/400) + 1
	# perform a runing median (k = 51) on Log2Ratio to estimate the number of peaks. More efficient than if applied on Log2ratios directly
	K = 51
	run.lr <- runmed(test.lr, k = K)	
	q1 <- cut[1]
	q2 <- cut[2]
	suppr <- which(run.lr<=q1 | run.lr>=q2)
	run.lr2 <- run.lr[-suppr]
	len <- floor(length(run.lr2)*0.25)
	model <- Mclust(run.lr2[sample(1:length(run.lr2), len)], G = G)

	p <- model$parameters$pro
	m <- model$parameters$mean
	s <- model$parameters$variance$sigmasq
	if(length(s)<length(m)) s <- rep(s, length(m))
	nG <- model$G
	p <- p[order(m)]
	cat("Props: ", p, "\n")
	s <- s[order(m)]
	m <- sort(m)
	pcor <- length(suppr)/length(test.lr)				# A vérifier ce truc !
	p <- p*(1-pcor)
	cat("Peak values: ", m, "\n")
	cat("Corrected Props: ", p, "\n")

	# merge Classes: depending on MergePeaks and MergeVal
	if(MergePeaks){
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
				m[C1] <- (p[C1]*m[C1] + p[C2]*m[C2])/(p[C1] + p[C2])
				s[C1] <- (p[C1]*s[C1] + p[C2]*s[C2])/(p[C1] + p[C2])
				p[C1] <- p[C1] + p[C2]
				m <- m[-C2]; s <- s[-C2]; p <- p[-C2]
				nG <- length(m)
				cat("means:", m, ", props:", p, "\n")
				}
			}
		}
	#

		dlr <- density(run.lr)
		max.dlr <- max(dlr$y)
		maxy <- max(dlr$y)*1.25
		plot(dlr$x, dlr$y, type = "h", col = "grey75", xlab = "Log2R", ylab = "Density", xlim = range(q1*Expand, q2*Expand), ylim = range(0, maxy))
		segments(x0 = c(q1, q2), y0 = c(0, 0), y1 = c(2, 2), lty = 2)
		dy <- cdy <- c()
		n <- length(lr[-suppr])
		for(i in 1:length(m)){
			tmp <- rnorm(n*p[i], m[i], sqrt(s[i])) # rnorm(n*p[i], m[i], sqrt(s[i]))
			tmp.d <- density(tmp, na.rm = T)
			dy <- c(dy, max(tmp.d$y*p[i]))
			cdy <- c(cdy, sum(tmp.d$y*p[i]))
			lines(I(tmp.d$y*p[i])~tmp.d$x, lwd = 2, col = i)
			# segments(x0 = m[i], x1 = m[i]+0.15, y0 = dy[i], lwd = 2)
			# text(x = m[i], y = dy[i]+max(dlr$y)/10, labels = round(m[i], 3), cex = 1.25)
			}

	cat("dy: ", dy, "\n")
	cat("Thresh: ", max(dy)*peakThresh, "\n")
	cat("Ratio dy/m: ", dy/m, "\n")
	peakIndex <- which.min(abs(m))
	# if (method == "Left") peakIndex <- which(dy>=max(dy)*peakThresh & m<=m[which.max(dy)])
	if (method == "Left") peakIndex <- which(dy>=max(dy)*peakThresh & m<=0)
	# if (method == "Zero") peakIndex <- which(dy>=max(dy)*peakThresh & abs(m)<abs(mean(test.lr))+2*sd(test.lr)/sqrt(length(test.lr)))
	if (method == "Zero") {peakIndex <- which(dy>=max(dy)*peakThresh); peakIndex <- peakIndex[which.min(abs(m)[peakIndex])]}
	if (method == "Right") peakIndex <- which(dy>=max(dy)*peakThresh & m>=m[which.max(dy)])

	# correct <- m[which.max(dy)]
	# if(length(peakIndex)>0) correct <- m[max(peakIndex)]
	if(length(peakIndex)==0) peakIndex <- which.max(dy)
	correct <- m[max(peakIndex)]

	symbols(correct, min(dy[max(peakIndex)] + max(dlr$y)/8, maxy), rectangle = matrix(c((q2-q1)/4, 0.2), 1, 2), inches = F, bg = "grey80", lty = 0, add = T)
	for(i in 1:length(m)){
		text(x = m[i], y = min(dy[i]+ max(dlr$y)/8, maxy), labels = round(m[i], 3), cex = 1.25, col = ifelse(i == max(peakIndex), "red3", "black"))
		}

	title(main = paste(cgh$Safir.Id, "\nEM centralization: correction.factor =", round(correct, 5)))
	if(Save) savePlot(filename = paste(PathPlot, cgh$Safir.Id, "_", cgh$BarCode, "_", cgh$analysisDate, "_CentralizeProfile", sep = ""), type="png")

	cat("n.peaks = ", nG, "peakIndex = ", peakIndex, "\n") 
	cat("Valeur de centrage EM: correction factor = ", correct, "\n")

	# return the full cgh objet containing the Log2Ratios adjusted for Cy and CG, and centralized by EM
	lr = lr - correct
	cgh$EM.grp <- nG
	cgh$EM.centralization <- correct
	}
	if(is.Agilent) cgh$cgh <- cbind.data.frame(cgh$cgh[,1:7], Log2Ratio = lr)	# for Agilent
	if(!is.Agilent) cgh$cgh$Log2Ratio <- lr										# for Affy
	rm(Data)
	return(cgh)
}
