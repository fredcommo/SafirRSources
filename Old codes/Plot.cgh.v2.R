
Plot.cgh.v2 <- function(cgh, greq, hg19.info, cutGainLoss = 10, cutAmpDel = 85, probes.col = "grey90", gain.col = "dodgerblue2", loss.col = "red3", flat.col = "grey55", expand = 0.5){

##### Marquage des gènes
tag.genes <- function(Greq, expand){
	lwd.seg = 2
	x1 <- quantile(Greq$loc.start, 0)*0.95
	x2 <- quantile(Greq$loc.end, 1)*1.05
	ntext <- nrow(Greq)
	g.start <- Greq$Genom.start

	gnames <- as.character(Greq$Symb)
	xtext <- seq(x1, x2, len = ntext)
	ytext <- seg.values <- Greq$seg.mean
	ytext <- ytext + max.val/2*abs(ytext)/ytext
	ytext <- ifelse(abs(ytext) < max.text, ytext, max.text*ytext/abs(ytext))
	textCol <- ifelse(seg.values < (-cutoff), loss.col, ifelse(seg.values > cutoff, gain.col, flat.col))

	# segments(x0 = g.start, x1 = g.start, y0 = seg.values, y1 = ifelse(abs(seg.values*expand)> abs(ytext), ytext, seg.values*expand), lwd = lwd.seg)
	# segments(x0 = xtext, x1 = g.start, y0 = ifelse(abs(seg.values*expand)> abs(ytext), ytext, seg.values*expand), lwd = lwd.seg)
	# segments(	x0 = xtext, x1 = xtext, 
	#		y0 = ifelse(abs(seg.values*expand)> abs(ytext), ytext, seg.values*expand), 
	#		y1 = ifelse(abs(ytext) == max.text, ytext, ifelse(ytext>0, ytext-0.2, ytext+0.2)),
	#		lwd = lwd.seg)

	# segments(x0 = g.start, x1 = g.start, y0 = seg.values, y1 = ifelse(abs(ytext)== max.text, ytext, ytext*expand), lwd = lwd.seg)
	segments(x0 = g.start, x1 = g.start, y0 = seg.values, y1 = ifelse(ytext== seg.values, ytext, ytext*expand), lwd = lwd.seg)
	# segments(x0 = xtext, x1 = g.start, y0 = ifelse(abs(ytext)== max.text, ytext, ytext*expand), lwd = lwd.seg)
	segments(x0 = xtext, x1 = g.start, y0 = ifelse(ytext== seg.values, ytext, ytext*expand), lwd = lwd.seg)
	segments(	x0 = xtext, x1 = xtext, 
			y0 = ifelse(ytext== seg.values, ytext, ytext*expand), # y0 = ifelse(abs(ytext)== max.text, ytext, ytext*expand),
			y1 = ifelse(ytext == seg.values, ytext, ifelse(ytext>0, ytext-max.val/10, ytext+max.val/10)), # y1 = ifelse(abs(ytext) == max.text, ytext, ifelse(ytext>0, ytext-max.val/10, ytext+max.val/10)),
			lwd = lwd.seg)

	symbols(xtext, ytext, rectangles = matrix(c(7e7, max.val/10), nrow = length(xtext), ncol = 2, byrow = T), inches = F, bg = "white", lty = 0, add = T)
	text(xtext, ytext, labels = gnames, cex = 0.75, font = 2)
	text(xtext, ifelse(ytext>0, ytext-max.val/20, ytext+max.val/20), labels = "---", cex = 1.5, col = textCol, font = 2)
	}
##### End function

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
	# LR <- greq$seg.mean

	# Définition des cutoff
		cutoff = log2(1 + cutGainLoss/100)
		cutoff2 = log2(1 + cutAmpDel/100)		# non utilisé (pour le moment)
		max.val <- max(abs(lr))*0.75
		max.text <- max.val*0.9 
		min.text <- max.val*0.25

	# Graphique principal
		plot(lr~pos, pch = 19, cex = 0.25, col = probes.col, ylim = range(-max.val, max.val), xlab = "Genome positioning", ylab = "Log2(ratio)")	# Valeurs des sondes
		abline(v = cum.len[1:23]*1000, col = "lightblue", lty = 2)														# ligne centrale = 0
		lines(pos, seg.val, col = "royalblue2", lty = 2)															# visualisation de la segmentation
		pCol <- ifelse(seg.val < (-cutoff), loss.col, ifelse(seg.val > cutoff, gain.col, flat.col))								# définition des couleurs G/L/N
		points(pos, seg.val, pch = "+", col = pCol, cex = 0.5)														# visualisation des segments
		abline(h = 0, lty = 3)

		# Définition du titre						
		Samp <- paste(	"Safir.Id: ", cgh$Safir.Id, "  /  Sample.Id: ", cgh$BarCode, "\n", "Scan.date: ", cgh$ScanDate,
				" / Analysis: ", cgh$analysisDate, " / Genome inst.(%): ", cgh$Prop, " (Cutoff = ", cutGainLoss, ")", sep = "")								# Définition du titre
		maInfos <- paste(	"Prod:", cgh$LabId, "/ MArray:", cgh$ArraySet,
				"/ GenomeDataBase: Human Feb. 2009(GRCh37/hg19)", "/ Script:", cgh$Script)												# définition du sous-titre
		title(main = Samp, sub = maInfos)

	# identification des Chr
		text(0, max.val, labels = "chr", cex = 0.75, col = "royalblue2")
		text(cum.len[1]/2*1000, max.val, labels = 1, cex = 0.75, col = "royalblue2")
		for(i in 2:length(cum.centr)){
			x <- (hg19.info$length[i]/2 + cum.len[i-1])*1000
			text(x, max.val, labels = i, cex = 0.75, col = "royalblue2")
			}

	# Identification des gènes
		tag.genes(greq[greq$seg.mean>=0,], expand)
		tag.genes(greq[greq$seg.mean<0,], expand)
}

