NormAgilent.v2 <- function(cgh, is.Agilent = T, Cy = T, GC = T, Centr = T, span = 0.3, cut = c(0.05, 0.95), G = 3:7, peakThresh = 0.8, Left = F, Plot = T){

data <- cgh$cgh

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
	g2 <- loessFit(r, g, span = span)$fitted
	lr <- r - g2
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
	GC <- data[,9:18]

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
	K = 2*ceiling(length(lr)/800) + 1
	if(!is.Agilent) K = 2*ceiling(length(lr)/32000) + 1
	run.lr <- runmed(lr, k = K)	
	q1 <- quantile(run.lr, cut[1])
	q2 <- quantile(run.lr, cut[2])
	suppr <- which(run.lr<q1 | run.lr>q2)
	model <- Mclust(run.lr[-suppr], G = G)

	p <- model$parameters$pro
	m <- model$parameters$mean
	s <- model$parameters$variance$sigmasq
	nG <- model$G
	p <- p[order(m)]
	s <- s[order(m)]
	m <- sort(m)
	pcor <- length(suppr)/length(lr)
	p <- p*(1-pcor)

	dlr <- density(run.lr)
	max.dlr <- max(dlr$y)
	max.dx <- max(dlr$x)
	if(Plot){
		plot(dlr$x, dlr$y, type = "h", col = "grey75", xlab = "Log2R", ylab = "density", xlim = range(-max.dx, max.dx)*1.1)
		segments(x0 = c(q1, q2), y0 = c(0, 0), y1 = c(2, 2), lty = 2)
		}
	dy <- cdy <- c()
	n <- length(lr[-suppr])

	for(i in 1:length(m)){
		tmp <- rnorm(n*p[i], m[i], sqrt(s[i]))
		tmp.d <- density(tmp)
		dy <- c(dy, max(tmp.d$y*p[i]))
		# cdy <- c(cdy, sum(tmp.d$y*p[i]))
		if(Plot) lines(I(tmp.d$y*p[i])~tmp.d$x, lwd = 2, col = i)
		}

	peakIndex <- which.min(abs(m))
	if (Left) peakIndex <- which(dy>=max(dy)*peakThresh & m<m[which.max(dy)])

	correct <- m[which.max(dy)]
	if(length(peakIndex)>0) correct <- m[max(peakIndex)]
	# correct
	if(Plot){
		abline(v = correct, col = "darkred", lwd = 3)
		title(main = paste(Safir.Id, "\nEM centralization: correction.factor =", round(correct, 5)))
		}
	cat(nG, "peaks, ", "Valeur de centrage EM: correction factor = ", correct, "\n")

	lr = lr - correct
	cgh$EM.centralization <- correct
	}
	cgh$cgh <- cbind.data.frame(cgh$cgh[,1:7], Log2Ratio = lr)
	return(cgh)
}
