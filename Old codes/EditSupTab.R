EditSupTab <- function(stab = seg.table, CGH = cgh, Tech = "", Exclude.Y = T){

require(annotate)
require(tcltk)

setwd("D:/Projet Safir/Data SAFIR/Safir Output/Safir SupplGeneList")

census <- read.table("D:/Projet Safir/Data SAFIR/Arrays Infos/CensusTable_Annot_FC.txt", header = T, sep = "\t")
full <- read.table("D:/Projet Safir/Data SAFIR/Arrays Infos/022060_4x180K_hg19_20110628_GenomePos&GC_GNames_FC.txt", header = T, sep = "\t")

thresh = 1 	# Log2R cutoff
FileName <- paste(CGH$Safir.Id, CGH$BarCode, CGH$analysisDate, Tech, "SupplTable.html", sep = "_")

if(Exclude.Y) stab <- stab[-which(stab$chrom == 24),]

whichseg <- which(stab$seg.mean>=1)
nseg = length(whichseg)
PB <- tkProgressBar(title = "The R'atWork BaBar", min = 0, max = nseg, width = 500)
k = 0

Append = F
if(length(whichseg)>1) Append <- T
for(i in whichseg){

	Sys.sleep(0.1)
	# launch & increment the pBar
	k = k + 1
      setTkProgressBar(PB, k, label = paste("Segment #", i, ": ", k," of ", nseg, " done.", sep = ""))

	chr <- stab$chrom[i]
	start <- stab$loc.start[i]
	end <- stab$loc.end[i]
	seg.value <- stab$seg.mean[i]
	index.full <- which(full$GenomicPos>=start & full$GenomicPos<=end)
	Symbols <- as.character(unique(full$Symbol[index.full]))
	if(any(Symbols == "")) Symbols <- Symbols[-which(Symbols == "")]
	if(length(Symbols)>0){
		refseq <- c()
		for(j in 1:length(Symbols)) refseq <- c(refseq, as.character(full$RefSeq[which(full$Symbol == Symbols[j])][1]))

		is.census <- ifelse(as.character(Symbols) %in% census$Symb, "Yes", "")

		tmp.request <- GeneRequest.v2(as.character(Symbols), hg19.info)
		tmp.request <- cbind.data.frame(tmp.request, RefSeq = refseq, Cancer.Census.List = is.census)
		tmp.request <- tmp.request[order(tmp.request$Symb),]
		tmp.request <- rbind.data.frame(tmp.request, rep("", ncol(tmp.request)))
		GeneId <- as.data.frame(tmp.request$GeneId)
		htmlpage(GeneId, filename = FileName, title = paste("Segment on Chr", chr, ": ", start, "-", end, ", Log2R = ", seg.value, sep = ""),
				othernames = tmp.request[,c(2:3, 6:8, 11:14)],
				table.head = c("GeneId", colnames(tmp.request[,c(2:3, 6:8, 11:14)])), repository = list("en"),
				append = Append)
		}
	}
	close(PB)
	rm(census, full)
	cat("Supplementary table saved: ", paste(getwd(), FileName, sep = "/"), "\n")
}

