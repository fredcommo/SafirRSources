NormCGH.v2 <- function(cgh, is.Agilent = TRUE, Cy = TRUE, GC = TRUE, Centr = TRUE, span = 0.3, cut = c(0.05, 0.95), G = 3:7, peakThresh = 0.5, method = "Zero", Plot = TRUE, Save = TRUE){

if(!Plot) Save = FALSE
data <- cgh$cgh
PathPlot <- "D:/Projet Safir/Data SAFIR/Safir Output/Safir CentralizeProfiles/"

	# Pr�Safir
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
	# centrage � valider +++
cgh$EM.centralization <- "No"
if(Centr){
	cat("Centrage EM\n")
	lr <- lr - median(lr, na.rm = T)
	test.lr <- lr
	if(!is.Agilent){
		test.lr <- lr[seq(1, length(lr), by = 10)]
		}
	K = 2*ceiling(length(test.lr)/400) + 1
	# K = 51
	run.lr <- runmed(test.lr, k = K)	
	q1 <- quantile(run.lr, cut[1])
	q2 <- quantile(run.lr, cut[2])
	suppr <- which(run.lr<q1 | run.lr>q2)
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
	pcor <- length(suppr)/length(test.lr)
	p <- p*(1-pcor)
	cat("Peak values: ", m, "\n")
	cat("Corrected Props: ", p, "\n")

	dlr <- density(run.lr)
	subdlr <- density(run.lr2)
	max.dlr <- max(dlr$y)
	min.dx <- min(dlr$x, na.rm = T)
	if(Plot){
		plot(dlr$x, dlr$y, type = "h", lwd = 3, col = "grey95", xlab = "Log2R", ylab = "density", xlim = range(min.dx, -min.dx)*1.5, ylim = range(subdlr$y))
		points(subdlr$x, subdlr$y, type = "h", lwd = 3, col = "grey75")
		segments(x0 = c(q1, q2), y0 = c(0, 0), y1 = c(2, 2), lty = 2)
		}
	dy <- cdy <- c()
	n <- len

	for(i in 1:length(m)){
		tmp <- rnorm(n*p[i], m[i], sqrt(s[i]))
		tmp.d <- density(tmp)
		dy <- c(dy, max(tmp.d$y*p[i]))
		# cdy <- c(cdy, sum(tmp.d$y*p[i]))
		if(Plot) lines(I(tmp.d$y*p[i])~tmp.d$x, lwd = 2, col = i)
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

	correct <- m[which.max(dy)]
	if(length(peakIndex)>0) correct <- m[max(peakIndex)]
	# correct
	if(Plot){
		abline(v = correct, col = "darkred", lwd = 3)
		title(main = paste(Safir.Id, "\nEM centralization: correction.factor =", round(correct, 5)))
		if(Save) savePlot(filename = paste(PathPlot, cgh$Safir.Id, "_", cgh$BarCode, "_", cgh$analysisDate, "_CentralizeProfile", sep = ""), type="png")
		}
	cat(nG, "peaks, ", "Valeur de centrage EM: correction factor = ", correct, "\n")

	lr = lr - correct
	cgh$EM.grp <- nG
	cgh$EM.centralization <- correct
	}
	if(is.Agilent) cgh$cgh <- cbind.data.frame(cgh$cgh[,1:7], Log2Ratio = lr)
	if(!is.Agilent) cgh$cgh$Log2Ratio <- lr
	return(cgh)
}
