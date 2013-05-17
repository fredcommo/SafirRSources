EM.Centr <- function(lr, Id = cgh$Safir.Id, cut = c(0.001, 0.999), G = 3:7, MergePeaks = TRUE, peakThresh = 0.85, MergeVal = 0.075, method = "Left", Expand = 1.75)
{
	cat("Centrage EM\n")
	newlr <- lr - median(lr, na.rm = T)
	test.lr <- newlr
	if(!is.Agilent){
		test.lr <- newlr[seq(1, length(newlr), by = 10)]
		}
	# K = 2*ceiling(length(test.lr)/400) + 1
	K = 51
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

	# merge Classes
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

	p <- p[order(m)]
	s <- s[order(m)]
	m <- sort(m)
	pcor <- length(suppr)/length(lr)
	p <- p*(1-pcor)

	par(mfrow = c(1, 2), cex.main = 1.5, cex.axis = 1.25, cex.lab = 1.25)
		dlr <- density(run.lr)
		max.dlr <- max(dlr$y)
		plot(dlr$x, dlr$y, type = "h", col = "grey75", xlab = "Log2R", ylab = "Density", xlim = range(q1*Expand, q2*Expand), ylim = range(0, max(dlr$y*1.25)))
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

	symbols(correct, dy[max(peakIndex)] + max(dlr$y)/6, rectangle = matrix(c((q2-q1)/4, 0.2), 1, 2), inches = F, bg = "grey80", lty = 0, add = T)
	for(i in 1:length(m)){
		text(x = m[i], y = dy[i]+ max(dlr$y)/6, labels = round(m[i], 3), cex = 1.25, col = ifelse(i == max(peakIndex), "red3", "black"))
		}
	# title(main = paste(Id, ": Peak assignment\n", "EM factor = ", round(correct, 4), sep = ""))
	return(list())

}