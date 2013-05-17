QCprofile.v2 <- function(CNA, Ncgh){

# En travaux

segtable <- CNA$output
lr <- Ncgh$cgh$Log2Ratio
pos <- Ncgh$cgh$GenomicPos
chr <- Ncgh$cgh$ChrNum
Rmed <- runmed(lr, k = floor(length(lr)/800)*2 + 1)

dL2R <- dRmed <- sSegm <- NSegm <- Pos <- Chr <- c()
for(i in 1:nrow(segtable)){
	x0 <- segtable$loc.start[i]
	x1 <- segtable$loc.end[i]
	y0 <- segtable$seg.med[i]
	tmp.lr <- lr[pos>=x0 & pos<=x1]
	tmp.rmed <- Rmed[pos>=x0 & pos<=x1]
	Pos <- c(Pos, pos[pos>=x0 & pos<=x1])
	Chr <- c(Chr, segtable$chrom[i])
	dlr <- tmp.lr - y0
	drmed <- tmp.rmed - y0
	dL2R <- c(dL2R, dlr)
	dRmed <- c(dRmed, drmed)
	sSegm <- c(sSegm, var(dlr))
	NSegm <- c(NSegm, rep(i, length(dlr)))
	}
	dL2Rs <- var(dL2R, na.rm = T)
	dRmeds <- nrow(segtable)*sum(dRmeds^2)		# en travaux
	Chr <- as.factor(Chr)

	#layout.show(nf)
	nf <- layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE), c(1,1,1), c(1, 1, 1), TRUE) 
	
	# mar = c(bottom, left, top, right)
	par(mar = c(5, 4, 1.5, 1), cex.main = 1.5, cex.axis = 1.25, cex.lab = 1.25) 
	plot(Pos, dL2R, type = "h", col = "grey75", ylim = range(-5, 5), xlab = "Genome positioning", ylab = "dist(Segment,Log2R)", main = paste(Ncgh$Safir.Id, "- QC profile"))
	abline(h = 0, lwd = 2, col = "red3")
	for(j in 1:nlevels(Chr)) abline(v = max(Pos[Ncgh$cgh$ChrNum == levels(Chr)[j]]), lty = 3, col = "grey45")
	legend("topleft", legend = paste("dL2Rs =", round(dL2Rs, 4)), bty = "n", cex = 1.5)

	par(mar=c(4.25, 4, 1.5, 1)) 
	# den <- density(dL2R, na.rm = T)
	plot(seq(1, nrow(segtable)), sSegm, type = "n", xlab = "Segment number", ylim = range(0, 2), ylab = "Var(dist)")
	abline(v = cumsum(table(Chr))+0.5, lty = 2, col = "darkslategray3")
	points(seq(1, nrow(segtable)), sSegm, type = "l", lwd = 3, col = "grey30")
	
	par(mar=c(4.25, 1.25, 1.5, 1))
	Col <- grey(sample(seq(0.1, 0.9, len = nlevels(Chr))))
	boxplot(dL2R~NSegm, border = "white", ylim = range(-5, 5))
	abline(v = cumsum(table(Chr))+0.5, lty = 2, col = "darkslategray3")
	boxplot(dL2R~NSegm, add = T, col = Col[Chr], outpch = 19, outcex = 25/nrow(segtable), outcol = Col[Chr], xlab = "Segment number", ylab = "dist(Segment,Log2R)")

	par(op) 	# default : par(mar=c(5,4,2,2))
	cat("dL2R = ", dL2Rs, "\n")
	Ncgh$nSegm <- nrow(segtable)
	Ncgh$dL2Rs <- dL2Rs
	return(Ncgh)
	}
