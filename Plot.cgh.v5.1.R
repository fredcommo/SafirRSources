
Plot.cgh.v5.1 <- function(Scgh, greq, hg19.info, Kmed = 3, cutGainLoss = 10, cutAmpDel = 85,
			probes.col = "grey90", gain.col = "dodgerblue3", loss.col = "red3", flat.col = "grey55",
			Y.expand = 0.5, Tag = TRUE, Dmin = 1.35e8, cex.text = 0.75, Tag.expand = 0.5, Equilibrate = TRUE,
			use.medians = T, useX11 = TRUE){

# Graphics device
if(useX11) X11(type = "dbcairo", width= 16, height= 8)

# cgh: the cgh object
# greq: the gene-requested table (see Safir list)
# hg19.info: hg19 data frame containing Chr lengths
# cutGainLoss: the gain/loss cutoff in %
# cutAmpDel: the Amplification/Deletion cutoff in %
# probes.col, gain.col, loss.col, flat.col: colors for Log2R probes and gain/loss/flat segments, respectively.
# Y.expand: a coefficient to define the Yaxis range, in proportion of Log2R range.
# Tag.expand: a coefficient to adjust gene names.
# Equilibrate: try to organize Norm-genes tags on both sides.


##### Marquage des genes
define.xpos <- function(x, xmin = 0, xmax = 3e9, dmin = Dmin, cex.text = cex.text){
	
	dmin = dmin*cex.text
	xmin = xmin
	xmax = xmax
	n <- length(x)
	if(dmin > (xmax - xmin)/(n-1)){
		dmin = (xmax - xmin)/(n-1)*0.9
		cat("Valeur minimale d'intervalles impossible a respecter\n")
		cat("Plus grande valeur possible: d = ", dmin, "\n")
		}
	d = x[-1]-x[-n]
	xprim <- x
	k = 0

	while(min(d)<dmin){
		for(i in 2:n){
			D <- abs(1/2*(xprim[i]-xprim[i-1]-dmin)) + dmin*0.25
			if(d[i-1]<dmin){
				xprim[i-1] <- max(xprim[i-1]-D, xmin)
				xprim[i] <- min(xprim[i]+D, xmax)
				}
			}
			k = k+1
			d = xprim[-1]-xprim[-n]
		}
	cat(k, "iterations\n")
	return(xprim)	
}

tag.genes <- function(Greq, Tag.expand, UseMed, dmin = Dmin, cex.text = cex.text){
	lwd.seg = 2
	x1 <- quantile(Greq$loc.start, 0)*0.75
	x2 <- quantile(Greq$loc.end, 1)*1.05
	ntext <- nrow(Greq)
	g.start <- Greq$Genom.start

	gnames <- as.character(Greq$Symb)
	xtext <- define.xpos(g.start, dmin = Dmin, cex.text = cex.text)
	seg.values <- Greq$seg.mean
	if(UseMed) seg.values <- Greq$seg.med

	w <- exp(-seq(-2, 2, len = length(seg.values))^2/(2*1))
	ytext <- seg.values + (w + max.val/3)*abs(seg.values)/seg.values

	ytext <- ifelse(abs(ytext) > min.text, ytext, min.text*ytext/abs(ytext))
	ytext <- ifelse(abs(ytext) < max.text, ytext, max.text*ytext/abs(ytext))
	textCol <- ifelse(seg.values < (-cutoff), loss.col, ifelse(seg.values > cutoff, gain.col, flat.col))

	segments(x0 = g.start, x1 = g.start, y0 = seg.values, y1 = ifelse(abs(ytext)<= abs(seg.values), seg.values, ytext*Tag.expand), lwd = lwd.seg)
	segments(x0 = xtext, x1 = g.start, y0 = ifelse(abs(ytext)<= abs(seg.values), seg.values, ytext*Tag.expand), lwd = lwd.seg)
	segments(	x0 = xtext, x1 = xtext, 
			y0 = ifelse(abs(ytext)<= abs(seg.values), seg.values, ytext*Tag.expand),
			y1 = ytext,
			lwd = lwd.seg)

	symbols(xtext, ytext, rectangles = cbind(nchar(gnames)*2.95e7*cex.text, max.val/(16*cex.text)), inches = F, bg = "white", lty = 0, add = T)
	text(xtext, ytext, labels = gnames, cex = cex.text, font = 2)
	text(xtext, ifelse(ytext>0, ytext-max.val/25, ytext+max.val/25), labels = "---", cex = 2, col = textCol, font = 2)
	}
##### End function

data <- Scgh$cgh

	# position genom., valeurs log2R, valeurs Segments
		pos <- data$GenomicPos
		lr <- data$Log2Ratio
		seg.val <- data$Segm
		chr = data$ChrNum

	# infos chrs, centres
		cum.len <- cumsum(hg19.info$length)
		centr <- hg19.info$centromere
		cum.centr <- c(centr[1], cum.len[-24] + centr[2:24])*1000

	# valeurs segments des genes SAFIR
	# LR <- greq$seg.mean
		if(Equilibrate){
			nSeg <- nrow(greq)
			nGain <- which(greq$Status == "Gain" | greq$Status == "Ampli")
			nNorm <- which(greq$Status == "Norm")
			nLoss <- which(greq$Status == "Loss")
			len <- length(nNorm)-length(nGain)
			if(len<=0){
				greq$seg.mean[nNorm] <- greq$seg.med[nNorm] <- -1e-3
				}
			if(len>0){
				index <- sample(nNorm, floor(nSeg/2)-length(nGain))
				greq$seg.mean[index] <- greq$seg.med[index] <- 1e-3
				greq$seg.mean[setdiff(nNorm, index)] <- greq$seg.med[setdiff(nNorm, index)] <- -1e-3
				}
			}

	# Definition des cutoff
		cutoff = log2(1 + cutGainLoss/100)
		cutoff2 = log2(1 + cutAmpDel/100)		# non utilise (pour le moment)
		max.val <- max(abs(lr))*Y.expand
		max.text <- max.val*0.9 
		min.text <- quantile(seg.val, probs = c(0.1, 0.9))
		min.text <- max(abs(min.text)) + 0.5

	# Graphique principal
		pSamp <- sample(1:length(lr), length(lr)/10)
		Rmed1 <- lr
		if(Scgh$Tech == "Affy") Rmed1 <- runmed(lr, k = 2*Kmed+1)
		Rmed2 <- runmed(lr, k = floor(length(lr)/800)*2 + 1)
		pCol <- ifelse(seg.val < (-cutoff), "thistle2", ifelse(seg.val > cutoff, "lightsteelblue1", "grey90"))																											# definition des couleurs G/L/N
		plot(Rmed1[pSamp] ~ pos[pSamp], pch = 19, cex = 0.25, col = pCol[pSamp], ylim = range(-max.val, max.val), xlab = "Genome positioning", ylab = "Log2(ratio)")	

		# Valeurs des sondes
		{
		if(any(chr == 24)) lines(pos[-which(chr == 24)], Rmed2[-which(chr == 24)], col = "black")
		else lines(pos, Rmed2, col = "black")
		}																								
				# visualisation de la segmentation
		pCol <- ifelse(seg.val < (-cutoff), loss.col, ifelse(seg.val > cutoff, gain.col, flat.col))																																
		# definition des couleurs G/L/N
		pSamp <- sample(1:length(seg.val), length(seg.val)/5)
		points(pos[pSamp], seg.val[pSamp], pch = "+", col = pCol[pSamp], cex = 0.75)																																	# visualisation des segments
		abline(h = 0, lty = 3)
		# cat("Rmed.sd =", sd(Rmed), "Dcoef =", sd((Rmed - seg.val)^2), "LR.sd = ", sd(lr), "\n")

		# Definition du titre						
		Samp <- paste(	"Safir.Id: ", Scgh$Safir.Id, "  /  Sample.Id: ", Scgh$BarCode, "\n", "Scan.date: ", Scgh$ScanDate,
				" / Analysis: ", Scgh$analysisDate, " / Genome inst.(%): ", Scgh$Prop, " (Cutoff = ", cutGainLoss, ")", sep = "")																			# Definition du titre
		maInfos <- paste(	"Prod:", Scgh$LabId, "/ MArray:", Scgh$ArraySet,
				"/ GenomeDataBase: Human Feb. 2009(GRCh37/hg19)", "/ Script:", Scgh$Script)																														# definition du sous-titre
		title(main = Samp, sub = maInfos)

	# identification des Chr
		abline(v = cum.len[1:23]*1000, col = "grey50", lty = 2)																																									# ligne centrale = 0
		chr.col = "grey30"
		text(0, max.val, labels = "chr", cex = 0.75, col = chr.col)
		text(cum.len[1]/2*1000, max.val, labels = 1, cex = 0.75, col = chr.col)
		for(i in 2:length(cum.centr)){
			x <- (hg19.info$length[i]/2 + cum.len[i-1])*1000
			text(x, max.val, labels = i, cex = 0.75, col = chr.col)
			}

	# Identification des genes
		if(Tag){
			if(nrow(greq[greq$seg.mean>0,])>0) tag.genes(greq[greq$seg.mean>0,], Tag.expand, UseMed = use.medians, dmin = Dmin, cex.text = cex.text)
			if(nrow(greq[greq$seg.mean<=0,])>0) tag.genes(greq[greq$seg.mean<=0,], Tag.expand, UseMed = use.medians, dmin = Dmin, cex.text = cex.text)
			}
}
