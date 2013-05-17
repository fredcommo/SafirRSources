NormCGH.v3 <- function(cgh, is.Agilent = TRUE, Cy = TRUE, GC = TRUE, Centr = TRUE, span = 0.3, cut = c(0.01, 0.975), G = 3:5, MergePeaks = TRUE, Save = TRUE, Root = "E"){

# if(!Plot) Save = FALSE
data <- cgh$cgh
PathPlot <- paste(Root, ":/Projet Safir/Data SAFIR/Safir Output/Safir CentralizeProfiles/", sep = "")

	# PréSafir
# PathPlot <- "D:/Projet Safir/Data_PRESAFIR1/PreSafir_Output/PS_CentralizeProfiles/"

if(is.Agilent){
	g <- log2(data$gMedianSignal)				# Ref
	r <- log2(data$rMedianSignal)				# Test
	pos <- data$GenomicPos
	lr <- r - g
	}
else {lr <- cgh$cgh$Log2Ratio
	Cy = FALSE
	GC = FALSE
	}

############################################
	# Normalisation Cy5/Cy3 et calcul LogR
cgh$Cy.adjust <- "No"
if(Cy){
	cat("Normalisation Cy5/Cy3\n")
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
	is.GC <- which(substr(names(data), 1, 2) == "GC")
	GC <- data[,is.GC]

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
	# centrage à valider +++
cgh$EM.centralization <- "No"
if(Centr){
	cat("Centrage EM\n")
	Id <- cgh$Safir.Id
	lr <- lr - median(lr, na.rm = T)
	test.lr <- lr
	if(!is.Agilent){
		test.lr <- lr[seq(1, length(lr), by = 10)]
		}
	K = 21
	run.lr <- runmed(test.lr, k = K)
	q1 <- quantile(run.lr, cut[1])
	q2 <- quantile(run.lr, cut[2])
	suppr <- which(run.lr<q1 | run.lr>q2)
	run.lr2 <- run.lr[-suppr]
	len <- floor(length(run.lr2)*0.2)
	model <- Mclust(run.lr2[sample(1:length(run.lr2), len)], G = G)

	p <- model$parameters$pro
	m <- model$parameters$mean
	s <- model$parameters$variance$sigmasq
	if(length(s)<length(m)) s <- rep(s, length(m))
	nG <- model$G

	# merge Classes
	if(MergePeaks){
		Mdist <- matrix(0, nG, nG)
		for(i in 1:nG)
			for(j in 1:nG){
			Mdist[i,j] <- abs(m[i] - m[j])
			}
		diag(Mdist) <- NA
		Raw <- ceiling(which(Mdist<0.05)/nlevels(c))
		if(length(Raw)!=0){
			C1 <- Raw[1]
			C2 <- Raw[2]
			m[C1] <- (p[C1]*m[C1] + p[C2]*m[C2])/(p[C1] + p[C2])
			s[C1] <- (p[C1]*s[C1] + p[C2]*s[C2])/(p[C1] + p[C2])
			p[C1] <- p[C1] + p[C2]
			nG <- nG - length(Raw) + 1
			m <- m[-C2]; s <- s[-C2]; p <- p[-C2]
			}
		}
	#

	p <- p[order(m)]
	s <- s[order(m)]
	m <- sort(m)
	pcor <- length(suppr)/length(lr)
	p <- p*(1-pcor)

	par(mfrow = c(1, 2), cex.main = 1.5, cex.axis = 1.25, cex.lab = 1.25)
		dlr <- density(run.lr)
		max.dlr <- max(dlr$y)
		plot(dlr$x, dlr$y, type = "h", col = "grey75", xlab = "Log2R", ylab = "Density", xlim = range(-1.5, 1.5))
		segments(x0 = c(q1, q2), y0 = c(0, 0), y1 = c(2, 2), lty = 2)
		dy <- cdy <- c()
		n <- length(lr[-suppr])
		for(i in 1:length(m)){
			tmp <- rnorm(n*p[i], m[i], sqrt(s[i]))
			tmp.d <- density(tmp)
			dy <- c(dy, max(tmp.d$y*p[i]))
			cdy <- c(cdy, sum(tmp.d$y*p[i]))
			lines(I(tmp.d$y*p[i])~tmp.d$x, lwd = 2, col = i)
			}
	title(main = paste(Id, ": Peak assignment", sep = ""))

	peaks <- order(dy, decreasing = T)[1:min(nG, 3)]
	R <- 2^(m[peaks])
	R <- sort(R)
	CN <- seq(1, (min(nG, 3)))/2-1 
	plot(CN, R, xlim = range(-0.75, 0.75), ylim = range(0, max(R))*1.25, pch = 19, cex = 1.5, col = "red3", xlab = "Copy number", ylab = "Ratio Sample/ref")
	summary(lm.test <- lm(R ~ CN))
	estim <- coef(lm.test)[2]*2
	estim <- min(estim, 1)
	abline(lm.test, lwd = 3, col = "violetred")
	title(main = paste(Id, ": CN line fit", sep = ""))
	text(0, 0.25, labels = paste("Clonal fraction: ", round(estim, 2)), cex = 1.25, font = 2, bty = "n")
	par(op)
	if(Save) savePlot(filename = paste(PathPlot, cgh$Safir.Id, "_", cgh$BarCode, "_", cgh$analysisDate, "_CentralizeProfile", sep = ""), type="png")

	M <- log2(R[which(CN == 0)])

	lr <- lr - M
	cgh$EM.grp <- nG
	cgh$EM.centralization <- M
	cgh$ClonalFraction <- round(estim, 2)
	cat("EM factor = ", M, "\nClonale estimation = ", estim, "\n")
	}
	if(is.Agilent) cgh$cgh <- cbind.data.frame(cgh$cgh[,1:7], Log2Ratio = lr)
	if(!is.Agilent) cgh$cgh$Log2Ratio <- lr
	return(cgh)
}
