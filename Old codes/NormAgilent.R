NormAgilent <- function(cgh, Cyan = T, GC = T, Centr = T, span = 0.3, cut = c(0.05, 0.95), G = 3){

data <- cgh$cgh

g <- log2(data$gMedianSignal)				# Ref
r <- log2(data$rMedianSignal)				# Test
pos <- data$GenomicPos
lr <- r - g

	# Normalisation Cy5/Cy3 et calcul LogR
if(Cyan){
	cat("Normalisation Cy5/Cy3\n")
	g2 <- loessFit(r, g, span = span)$fitted
	lr <- r - g2
	lr <- lr - median(lr, na.rm = T)
	}

	# Ajustement GC
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


	# centrage à valider +++
if(Centr){
	cat("Centrage\n")
	# cutoff = cut
	# drm <- gc.lr[which(abs(gc.lr) <= cutoff)]
	drm <- lr
	drm[lr<quantile(lr, cut[1])] <- NA
	drm[lr>quantile(lr, cut[2])] <- NA
	den <- density(drm, na.rm = T)

	model <- Mclust(as.numeric(drm[!is.na(drm)]), G = G)	# G = 3 ?
	means <- model$parameters$mean
	props <- model$parameters$pro
	v <- model$parameters$variance$sigmasq
	classif <- model$classification
	# centr1 <- median(drm[classif == which.min(abs(means))], na.rm = T)
	centr1 <- means[which.max(props)]
	centr2 = median(drm, na.rm = T)

	# plot(den)
	# n <- length(drm)
	# for(i in 1:length(means)){
	#	tmp <- rnorm(n*props[i], means[i], sqrt(v[i]))
	#	tmp.d <- density(tmp)
	#	lines(I(tmp.d$y*props[i])~tmp.d$x, lwd = 3, col = i)
	#	}

	# abline(v = c(centr1, centr2), col = c("blue", "red"))
	# legend("topleft", legend = paste(c("EM:", "Median:"), c(signif(centr1, 3), signif(centr2, 3))), lwd = 2, col = c("blue", "red"), bty = "n")
	cat("Means: ", round(means, 3), "\nProps: ", round(props, 3), "\n")
	# cat("Valeurs de centrage: EM = ", centr1, "\n", "Valeurs de centrage: mediane = ", centr2, "\n")

	lr = lr - centr1
	}
	cgh$cgh <- cbind.data.frame(cgh$cgh[,1:7], Log2Ratio = lr)
	return(cgh)
}
