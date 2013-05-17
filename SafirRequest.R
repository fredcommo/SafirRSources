

SafirRequest <- function(Scgh, GeneList, HG = hg19.info, seg.cna.obj){

		# Requête web des coordonnées génomiques: GeneRequest(GeneList, hg19.ChrLen)
	greq <- GeneRequest.v2(GeneList, HG)					
	for(i in c(5, 7:10)) greq[,i] <- as.numeric(as.vector(greq[,i]))
	greq <- greq[order(greq$Chr, greq$Genom.start),]
	rownames(greq) <- seq(1, nrow(greq))

		# Affectation des valeurs de segmentation & statut G/L/N: greq.val(cgh, greq, cutmethod = "delta" or a value)
	greq <- greq.val.v2(Scgh, greq, seg.cna.obj)
	greq <- cbind.data.frame(Safir.Id = Scgh$Safir.Id, greq)
	return(greq)
}