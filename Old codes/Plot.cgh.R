Plot.cgh <- function(cgh, greq, hg19.info, cutGainLoss = 10, cutAmpDel = 85, probes.col = "grey90", gain.col = "dodgerblue2", loss.col = "red3", flat.col = "grey55"){

data <- cgh$cgh

# position genom., valeurs log2R, valeurs Segments
pos <- data$GenomicPos
lr <- data$Log2Ratio
seg.val <- data$Segm

# infos chrs, centres
cum.len <- cumsum(hg19.info$length)
centr <- hg19.info$centromere
cum.centr <- c(centr[1], cum.len[-24] + centr[2:24])*1000

# valeurs segments des gènes SAFIR
LR <- greq$seg.mean

# bornes axe y
max.val <- max(abs(lr))*0.75
expand <- max.val*5.5

		# Définition des cutoffs
	cutoff = log2(1 + cutGainLoss/100)
	cutoff2 = log2(1 + cutAmpDel/100)		# non utilisé (pour le moment)

par(cex.main = 1.5, cex.sub = 0.75)
	plot(lr~pos, pch = 19, cex = 0.25, col = probes.col, ylim = range(-max.val, max.val), xlab = "Genome positioning", ylab = "Log2(ratio)")	# Valeurs des sondes
	abline(v = cum.len[1:23]*1000, col = "lightblue", lty = 2)														# ligne centrale = 0
	lines(pos, seg.val, col = "royalblue2", lty = 2)															# visualisation de la segmentation
	pCol <- ifelse(seg.val < (-cutoff), loss.col, ifelse(seg.val > cutoff, gain.col, flat.col))								# définition des couleurs G/L/N
	points(pos, seg.val, pch = "+", col = pCol, cex = 0.5)														# visualisation des segments
	abline(h = 0, lty = 3)																				# repère central lr = 0
	Samp <- paste(	"Safir.Id: ", cgh$Safir.Id,
				"  /  Sample.Id: ", cgh$BarCode, "\n",
				"Scan.date: ", cgh$ScanDate,
				" / Analysis: ", cgh$analysisDate,
				" / Genome inst.(%): ", cgh$Prop, " (Cutoff = ", cutGainLoss, ")", sep = "")								# Définition du titre
	maInfos <- paste(	"Prod:", cgh$LabId,
				"/ MArray:", cgh$ArraySet,
				"/ GenomeDataBase: Human Feb. 2009(GRCh37/hg19)",
				"/ Script:", cgh$Script)												# définition du sous-titre
	title(main = Samp, sub = maInfos)

	# identification des Chr
	text(0, max.val, labels = "chr", cex = 0.75, col = "royalblue2")
	text(cum.len[1]/2*1000, max.val, labels = 1, cex = 0.75, col = "royalblue2")
	for(i in 2:length(cum.centr)){
		x <- (hg19.info$length[i]/2 + cum.len[i-1])*1000
		text(x, max.val, labels = i, cex = 0.75, col = "royalblue2")
		}

	# identification des gènes d'intérêt
	s.greq <- split(greq, greq$Chr)
	output <- c()
	prev.x = 1
	prev.y = 1
	for(i in 1:length(s.greq)){
		# i = 1
		tmp <- s.greq[[i]]
		chr <- unique(tmp$Chr)
		n <- nrow(tmp)
		xdecale = ydecale = 0
	
		if(n>1){
			xdecale <- (2e8/3 - 2e8*seq(0, 1, len = n))
			ydecale <- seq(0, max.val/4, len = n)
			}

		for(j in 1:n){
			tag <- as.character(tmp$Symb[j])
			g.start <- as.numeric(as.vector(tmp$Genom.start[j]))
			g.end <- as.numeric(as.vector(tmp$Genom.end[j]))
			# seg.values <- cbind(tmp$LR1, tmp$LR2)[j,]
			seg.value <- tmp$seg.mean[j]
			pCol <- ifelse(seg.value < (-cutoff), loss.col, ifelse(seg.value > cutoff, gain.col, "grey40"))
			xtext <- g.start - xdecale[j]

			ytext <- min(max.val*0.85, max(abs(seg.value)*expand, max.val*0.45))*seg.value/abs(seg.value)
			if(abs(seg.value)<cutoff) ytext <- ytext*(-1)^j
			if(abs(prev.x - xtext)<1e8 & abs(prev.y - ytext)<1) ytext <- ytext + 0.5*ytext/abs(ytext)
			ytext <- min(max.val*0.85, abs(ytext))*ytext/abs(ytext)
			text(xtext, ytext, labels = tag, cex = 0.75, col = "navy", font = 2)
			segments(x0 = g.start, x1 = g.start, y0 = ytext*0.7, y1 = seg.value)
			segments(x0 = xtext, x1 = g.start, y0 = ytext*0.7, y1 = ytext*0.7)
			segments(x0 = xtext, x1 = xtext, y0 = ytext*0.7, y1 = ifelse(ytext>0, ytext-0.2, ytext+0.2))		# y1 = ytext*0.9
			text(xtext, ifelse(ytext>0, ytext-0.15, ytext+0.15), labels = "--", cex = 1.5, col = pCol, font = 2)	# y1 = ytext*0.9
			prev.x <- xtext
			prev.y <- ytext
			}
		}
	par(op)
}
