greq.val <- function(cgh, greq, cutmethod = "delta"){

data <- cgh$cgh
seg.val <- data$Segm
lr <- data$Log2Ratio
seg.start <- seg.cna.obj$output$loc.start
seg.end <- seg.cna.obj$output$loc.end
seg.mean <- seg.cna.obj$output$seg.mean	


if(is.numeric(cutmethod)) cutoff = cutmethod
else {
	delta <- c()
	for(i in 2:length(seg.val)){
		delta <- c(delta, abs(lr[i] - lr[i-1]))
		}
	cutoff <- median(delta)*0.25
	}

	# identification des gènes d'intérêt
	s.greq <- split(greq, greq$Chr)
	output <- status <- c()
	for(i in 1:length(s.greq)){
		tmp <- s.greq[[i]]
		chr <- unique(tmp$Chr)
		n <- nrow(tmp)
	
		for(j in 1:n){
			g.start <- as.numeric(as.vector(tmp$Genom.start[j]))
			g.end <- as.numeric(as.vector(tmp$Genom.end[j]))
			index.start <- which(seg.start <= g.start)						# !! Pb si gène à cheval sur 2 segments !!
			index.end <- which(seg.end >= g.end)							# !! Pb si gène à cheval sur 2 segments !!
			start.value <- seg.mean[index.start[length(index.start)]]
			end.value <- seg.mean[index.end[1]]
			seg.value <- ifelse(abs(start.value)>abs(end.value), start.value, end.value)
			status <- c(status, ifelse(seg.value>cutoff, "Gain", ifelse(seg.value<(-cutoff), "Loss", "Norm")))
			output <- rbind(output, c(start.value, end.value))
			}
		}
	FC <- ifelse(output<0, (-1)/2^output, 2^output)
	greq <- cbind.data.frame(greq, LR1 = output[,1], LR2 = output[,2], FC1 = FC[,1], FC2 = FC[,2], GL.status = status)
	return(greq)
}