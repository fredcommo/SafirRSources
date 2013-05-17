
MSegments <- function(seg.cna.obj, Ncgh){

	# Ajoute les valeurs medianes de segmentation (tukey.biweight) pour chaque segment

	Data <- Ncgh$cgh
	
	seg.start <- seg.cna.obj$output$loc.start
	seg.end <- seg.cna.obj$output$loc.end			
	seg.len <- seg.cna.obj$output$num.mark
	cum.seg.len <- cumsum(seg.len)
	s <- Data$ChrStart[cum.seg.len]
	e <- Data$ChrEnd[cum.seg.len]
	if(!is.null(e)) seg.end <- seg.end + (e-s)


	seg.med <- c()
	for(i in 1:length(seg.start)){
		index <- which(Data$GenomicPos>= seg.start[i] & Data$GenomicPos<=seg.end[i])
		tmp <- tukey.biweight(Data$Log2Ratio[index])
		seg.med <- c(seg.med, tmp)
		}
	seg.cna.obj$output <- cbind.data.frame(seg.cna.obj$output, seg.med = seg.med)
	rm(Data)
	return(seg.cna.obj)
}