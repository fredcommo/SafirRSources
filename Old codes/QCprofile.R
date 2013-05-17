QCprofile <- function(CNA.output, Ncgh){

segtable <- CNA.output$output
lr <- Ncgh$cgh$Log2Ratio
pos <- Ncgh$cgh$GenomicPos
chr <- Ncgh$cgh$ChrNum

dSegm <- NSegm <- Pos <- Chr <- c()
for(i in 1:nrow(segtable)){
	x0 <- segtable$loc.start[i]
	x1 <- segtable$loc.end[i]
	y0 <- segtable$seg.mean[i]
	tmp <- lr[pos>=x0 & pos<=x1]
	Pos <- c(Pos, pos[pos>=x0 & pos<=x1])
	Chr <- c(Chr, segtable$chrom[i])
	dlr <- tmp - y0
	dSegm <- c(dSegm, dlr)
	NSegm <- c(NSegm, rep(i, length(dlr)))
	}
	dSegms <- var(dSegm, na.rm = T)
	Chr <- as.factor(Chr)

	#layout.show(nf)
	nf <- layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE), c(1,1,1), c(1, 1, 1), TRUE) 
	# layout.show(nf)	# à vérifier

	# mar = c(bottom, left, top, right)
	par(mar = c(5, 4, 1.5, 1), cex.main = 1.5, cex.axis = 1.25, cex.lab = 1.25) 
	plot(Pos, dSegm, type = "h", col = "grey75", ylim = range(-5, 5), xlab = "Genome positioning", ylab = "dist(Segment,Log2R)", main = paste(Ncgh$Safir.Id, "- QC profile"))
	abline(h = 0, lwd = 2, col = "red3")
	for(j in 1:nlevels(Chr)) abline(v = max(Pos[Ncgh$cgh$ChrNum == levels(Chr)[j]]), lty = 3, col = "grey45")
	legend("topleft", legend = paste("dSegms =", round(dSegms, 4)), bty = "n", cex = 1.5)

	par(mar=c(4.25, 4, 1.5, 1)) 
	den <- density(dSegm, na.rm = T)
	plot(den$x, den$y, type = "l", lwd = 3, col = "grey75", xlim = range(-2, 2), xlab = "dist(Segment,Log2R)", ylab = "density")
	
	par(mar=c(4.25, 1.25, 1.5, 1))
	Col <- grey(sample(seq(0.1, 0.9, len = nlevels(Chr))))
	boxplot(dSegm~NSegm, border = "white", ylim = range(-5, 5))
	abline(v = cumsum(table(Chr))+0.5, lty = 2, col = "darkslategray3")
	boxplot(dSegm~NSegm, add = T, col = Col[Chr], outpch = 19, outcex = 25/nrow(segtable), outcol = Col[Chr], xlab = "Segment number", ylab = "dist(Segment,Log2R)")

	par(op) 	# default : par(mar=c(5,4,2,2))
	cat("dSegm = ", dSegms, "\n")
	Ncgh$nSegm <- nrow(segtable)
	Ncgh$dSegms <- dSegms
	return(Ncgh)
	}
