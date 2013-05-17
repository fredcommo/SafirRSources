
# Major changes : calls of values with different names, e.g. genomicPos instead of GenomicPos

WAAI.2 <- function(SegTable, Centr){

	SD <- sd(SegTable$Log2Ratio, na.rm = T)
	Score <- c()
	for(i in 1:24){
		pArm = paste(i, "p", sep = "")
		qArm = paste(i, "q", sep = "")
		p <- SegTable[which(SegTable$chrNum == i & SegTable$genomicPos <= Centr[i]), ]
		q <- SegTable[which(SegTable$chrNum == i & SegTable$genomicPos >= Centr[i]), ]
		npProbe <- nrow(p)
		nqProbe <- nrow(q)
		# if(npSeg == 0 | nqSeg == 0) pArm = qArm = i
		# tmpScore <- c()
		# waai = 0
		pScore = qScore = 0
	
		rmed <- runmed(p$Log2Ratio, k = 501)
		mpArm <- mean(rmed, na.rm = T)
		spArm <- sd(rmed, na.rm = T)
		s = mpArm/abs(mpArm)
		Q <- ifelse(s<0, quantile(rmed, 0.25), quantile(rmed, 0.75))
		pScore <- Q/spArm 
		
		Score <- rbind(Score, c(Chr = pArm, waai = pScore))
		# cat("Chr", i, "pArm", pScore, "\n")

		rmed <- runmed(q$Log2Ratio, k = 501)
		mqArm <- mean(rmed, na.rm = T)
		sqArm <- sd(rmed, na.rm = T)
		s = mqArm/abs(mqArm)
		Q <- ifelse(s<0, quantile(rmed, 0.25), quantile(rmed, 0.75))
		qScore <- Q/sqArm 
		
		Score <- rbind(Score, c(Chr = qArm, waai = qScore))
		# cat("Chr", i, "qArm", qScore, "\n")

		}

	Score <- as.data.frame(Score)
	Score$waai <- as.numeric(as.vector(Score$waai))
	# Score <- Score[-which(Score$Chr == "13p" | Score$Chr == "14p" | Score$Chr == "15p" | Score$Chr == "23p" | Score$Chr == "24p"),]	
	return(Score)
	}

