
# Use org.Hs.eg.db package (June 9, 2011)

GeneRequest.Org <- function(genelist, hg19.info, DB = "gene", verbose = TRUE, PB = T){

require(org.Hs.eg.db)

	gene.annot <- as.data.frame(gene.annot)
	colnames(gene.annot) <- c("Query", "Symb", "Name", "Org", "Chr", "Cytoband", "Chr.start", "Chr.end", "Genom.start", "Genom.end", "RangeGB.Id", "GeneId")
	if(PB) close(pb)
	return(gene.annot)
}

myList <- c("ERCC1", "EGFR", "ALK")
Symb <- myList
EG <- unlist(mget(Symb, org.Hs.egALIAS2EG))		# EntrezGene Id
Name <- unlist(mget(EG, org.Hs.egGENENAME))		# Name
Org <- org.Hs.egORGANISM					# Organism
Chr <- unlist(mget(EG, org.Hs.egCHR))			# Chr
CytoBand <- unlist(mget(EG, org.Hs.egMAP))		# Cytoband
ChrLoc <- unlist(mget(EG, org.Hs.egCHRLOC))		# Start position

Len <- org.Hs.egCHRLENGTHS[1:24]
names(Len)[names(Len) == "X"] <- 23
names(Len)[names(Len) == "Y"] <- 24
chrNum <- as.numeric(as.vector(names(Len)))
Len <- Len[order(chrNum)]
cumsum(as.numeric(Len))

source("E:\\Projet Safir\\Data SAFIR\\Safir R Sources\\GeneRequest.v4.R")
GeneRequest.v4(myList)
