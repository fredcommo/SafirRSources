SegCGH <- function(Ncgh, Undo.Splits = "sdundo", Alpha = 0.01, Undo.SD = 0.75){

	# Segmentation
	lr <- Ncgh$cgh$Log2Ratio
	chr <- Ncgh$cgh$ChrNum
	GenomPos <- Ncgh$cgh$GenomicPos

	cna.obj <- CNA(lr, chr, GenomPos, data.type = "logratio", sampleid = paste(Ncgh$Safir.Id, Ncgh$BarCode, sep = "_"))
	smooth.cna.obj <- smooth.CNA(cna.obj)
	seg.cna.obj <- segment(smooth.cna.obj, undo.splits = Undo.Splits, alpha = Alpha, undo.SD = Undo.SD)

	#Ajoute les valeurs de segmentation au fichier cgh: AddSegments(seg.cna.obj, cgh)
	Scgh <- AddSegments(seg.cna.obj, Ncgh)
	return(Scgh)
}