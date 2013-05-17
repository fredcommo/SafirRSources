
AddSegments <- function(seg.cna.obj, Ncgh, cutGainLoss = 10, use.medians = TRUE){

	# Ajoute les valeurs de segmentation pour chaque sonde/position

Data <- Ncgh$cgh

seg.start <- seg.cna.obj$output$loc.start
seg.end <- seg.cna.obj$output$loc.end			# !!!! verifier seg-end
seg.mean <- seg.cna.obj$output$seg.mean
if(use.medians) seg.mean <- seg.cna.obj$output$seg.med
seg.len <- seg.cna.obj$output$num.mark
cum.seg.len <- cumsum(seg.len)
s <- Data$ChrStart[cum.seg.len]
e <- Data$ChrEnd[cum.seg.len]
if(!is.null(e)) seg.end <- seg.end + (e-s)

# Proportion d'aberrations
cutoff <- log2(1 + cutGainLoss/100)
seg.delta <- seg.end - seg.start
Prop <- sum(seg.delta[which(abs(seg.mean)>=cutoff)])/sum(seg.delta)*100
Prop <- round(Prop, 2)

seg.val <- rep(0, nrow(Data))
for(i in 1:length(seg.mean)){
	index <- which(Data$GenomicPos>= seg.start[i] & Data$GenomicPos<=seg.end[i])
	seg.val[index] <- seg.mean[i]
	}
Data <- cbind.data.frame(Data, Segm = seg.val)
lr <- Data$Log2Ratio
Rmed <- runmed(lr, k = floor(length(lr)/800)*2 + 1)
segm <- Data$Segm
D <- Rmed - segm
Ncgh$ProfileQC <- sd(D^2)		# en travaux !
Ncgh$Prop <- Prop
Ncgh$cgh <- Data
rm(Data)
return(Ncgh)
}
