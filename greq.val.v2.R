greq.val.v2 <- function(Scgh, greq, cutGainLoss = 10, cutAmpDel = 85, use.medians = T){

data <- Scgh$cgh
seg.val <- data$Segm
lr <- data$Log2Ratio
seg.start <- seg.cna.obj$output$loc.start
seg.end <- seg.cna.obj$output$loc.end
seg.mean <- seg.cna.obj$output$seg.mean	


	# Definition des cutoffs
	GL = log2(1 + cutGainLoss/100)
	AmpDel = log2(1 + cutAmpDel/100)

	# identification des genes d'interet
	s.greq <- split(greq, greq$Chr)
	output <- new.greq <- c()

	for(i in 1:length(s.greq)){
		tmp <- s.greq[[i]]
		chr <- unique(tmp$Chr)
		n <- nrow(tmp)
	
		for(j in 1:n){
			g.start <- as.numeric(as.vector(tmp$Genom.start[j]))
			g.end <- as.numeric(as.vector(tmp$Genom.end[j]))

				# where is the gene start ?
			index.start <- which(seg.start <= g.start)						# !! Pb si gene a cheval sur 2 segments !!
			index.start <- index.start[length(index.start)]
			# start.value <- seg.cna.obj$output$seg.mean[index.start]

				# where is the gene end ?
			index.end <- which(seg.end >= g.end)						# !! Pb si gene a cheval sur 2 segments !!
			index.end <- index.end[1]
			# end.value <- seg.cna.obj$output$seg.mean[index.end]

			for(k in index.start:index.end){
				seg.num <- k
				seg.len <- seg.end[k] - seg.start[k]
				tmp.seg <- cbind(tmp[j,], SegmentNum = k, seg.cna.obj$output[k, 3:4], seg.len, seg.cna.obj$output[k, 5:7])
				new.greq <- rbind(new.greq, tmp.seg)
				}
			# output <- rbind(output, c(start.value, end.value))
			}
		}
	
	seg.values <- new.greq$seg.mean
	if(use.medians) seg.values <- new.greq$seg.med
	status <- ifelse(seg.values >= GL, "Gain", ifelse(seg.values <=(-GL), "Loss", "Norm"))
	status <- ifelse(seg.values >= AmpDel, "Ampli", ifelse(seg.values < (-AmpDel), "Del", status))
	FC <- ifelse(seg.values < 0, (-1)/2^seg.values, 2^seg.values)
	greq <- cbind.data.frame(new.greq, FC = FC, Status = status)
	return(greq)
}