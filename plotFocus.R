
# Main function
plotFocus <- function(Scgh.obj, geneList, hg19.info, E= 2.5e-3, useX11 = TRUE,...){
	myTab <- buildTab(Scgh.obj, geneList, hg19.info)
	sub <- buildSub(Scgh.obj, myTab, E)
	for(i in 1:length(sub)){
		X11(type = 'dbcairo')
		plot(sub[[i]]$Log2Ratio~sub[[i]]$GenomicPos, pch = 19, cex = 0.5, col = rgb(0.2, 0.2, 0.2, 0.2),
			xlab = 'Position on genome', ylab = 'Log2Ratio',
			main = paste('Focus on ', myTab$Symb[i], '\nChr', myTab$Chr[i], ':', paste(myTab$Chr.start[i], myTab$Chr.end[i], sep = "-"), sep = ''),...)
		lines(sub[[i]]$Segm~sub[[i]]$GenomicPos, lwd = 2)
		locate(myTab[i,])
		}
	return (myTab)
}
# End main function

buildTab <- function(Scgh.obj, geneList, hg19.info){
	myTab <- try(Request(geneList, hg19.info))			
	for(i in c(5, 7:10))
		myTab[,i] <- as.numeric(as.vector(myTab[,i]))
	myTab <- myTab[order(myTab$Chr, myTab$Genom.start),]
	myTab <- RequestVal(Scgh.obj, myTab, cutGainLoss = 10, cutAmpDel = 85)
	return(myTab)
}

buildSub <- function(Scgh.obj, myTab, E){
	sub <- function(x, start, end){
		return(which(x>=min(start, end)*(1-E) & x<=max(start, end)*(1+E)))
		}
	subProbes = list()
	for(i in 1:nrow(myTab)){
		start = myTab$Genom.start[i]
		end = myTab$Genom.end[i]
		subProbes[[i]] = Scgh.obj$cgh[sub(Scgh.obj$cgh$GenomicPos, start, end), ]
		}
	return(subProbes)
}

locate <- function(geneTable){
	segments(x0 = geneTable$Genom.start, y0 = geneTable$seg.med, x1 = geneTable$Genom.end, lwd = 10, col = "red")
	text((geneTable$Genom.start+geneTable$Genom.end)/2-5e5, max(0.25, geneTable$seg.med*1.25), labels = geneTable$Symb, cex = 1.25, font = 2, col = 'red3')
}
