dLRsd <- function(LR){
	n <- length(LR)
	V1 <- LR[-1]
	V2 <- LR[-n]
	dLR <- V2-V1
	q1 <- quantile(dLR, 0.25, na.rm = TRUE)
	q3 <- quantile(dLR, 0.75, na.rm = TRUE)
	s <- sd(dLR[which(dLR > q1 & dLR < q3)], na.rm = TRUE)/sqrt(2)
	return(s)
	}
