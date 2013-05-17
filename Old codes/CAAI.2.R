# cumlen <- c(0, cumsum(hg19.info$length))
# cumlen <- cumlen[-length(cumlen)]
# cumcentr <- cumsum(hg19.info$centromere)
# cumcentr <- cumcentr[-length(cumcentr)]
# Centr <- (cumlen + cumcentr)*1000

CAAI.2 <- function(SegTable, alpha = 1e7, thetaP = 5e-3, thetaQ = 1.2){
	Score <- c()
	for(i in 1:24){
		pArm = paste(i, "p", sep = "")
		qArm = paste(i, "q", sep = "")
		p <- SegTable[which(tmp$chrom == i & tmp$loc.start <= Centr[i]), ]
		q <- SegTable[which(tmp$chrom == i & tmp$loc.start >= Centr[i]), ]
		npSeg <- nrow(p)
		nqSeg <- nrow(q)
		# if(npSeg == 0 | nqSeg == 0) pArm = qArm = i
		tmpScore <- c()
		caai = 0
		pScore = qScore = 0
		if(npSeg>1){
			for(s in 2:npSeg){
				size <- p$loc.end[s] - p$loc.start[s-1]
				P <- tanh(alpha/size)
				val <- p$seg.med[s] - p$seg.med[s-1]
				Q <- tanh(abs(val)/thetaQ)
				W <- 1/2*(1 + tanh(10*(P - 1/2))/tanh(5))
				tmpScore <- rbind(tmpScore, c(P, Q, W))
				}
			if(!is.null(tmpScore)){
				if(nrow(tmpScore)>1) pScore <- sum(apply(tmpScore[,1:2], 1, min)*tmpScore[,3])
				else pScore <- min(tmpScore[1:2])*tmpScore[3]
				}
			}
		Score <- rbind(Score, c(Chr = pArm, caai = pScore))

		if(nqSeg>1){
			for(s in 2:nqSeg){
				size <- q$loc.end[s] - q$loc.start[s-1]
				P <- tanh(alpha/size)
				val <- q$seg.med[s] - q$seg.med[s-1]
				Q <- tanh(abs(val)/thetaQ)
				W <- 1/2*(1 + tanh(10*(P - 1/2))/tanh(5))
				tmpScore <- rbind(tmpScore, c(P, Q, W))
				}
			if(!is.null(tmpScore)){
				if(nrow(tmpScore)>1) qScore <- sum(apply(tmpScore[,1:2], 1, min)*tmpScore[,3])
				else qScore <- min(tmpScore[1:2])*tmpScore[3]
				}
			}
	 	Score <- rbind(Score, c(Chr = qArm, caai = qScore))
		}
	Score <- as.data.frame(Score)
	# Score <- Score[-which(Score$Chr == "13p" | Score$Chr == "14p" | Score$Chr == "15p" | Score$Chr == "23p" | Score$Chr == "24p"),]	
	return(Score)
	}

# CAAI.2(tmp)
