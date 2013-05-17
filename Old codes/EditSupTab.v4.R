EditSupTab.v4 <- function(stab = seg.table, CGH = Scgh, Tech = "", use.medians = TRUE, Select = c("Gain", "Loss", "Both"), Thresh.Gain = 1, Thresh.Loss = -1, CHR = 1:23, Restrict = FALSE, Root = "E"){

require(annotate)
require(tcltk)

# build a PGKB repository
	repofun <- function(ids, ...)
	paste("http://www.pharmgkb.org/gene/", ids, sep = "")
	setRepository("pgkb", repofun)

# build a Sanger Census Cancer repository
	repofun <- function(ids, ...)
	paste("http://www.sanger.ac.uk/perl/genetics/CGP/cosmic?action=gene&ln=", ids, sep = "")
	setRepository("census", repofun)

# build a Kegg repository
	repofun <- function(ids, ...)
	paste("http://www.genome.jp/dbget-bin/www_bget?hsa:", ids, sep = "")
	setRepository("kegg", repofun)


setwd(paste(Root, ":/Projet Safir/Data SAFIR/Arrays Infos", sep = ""))
census <- read.table("CensusTable_Annot_FC.txt", header = T, sep = "\t")
full <- read.table("022060_4x180K_hg19_20110628_GenomePos&GC_GNames_FC.txt", header = T, sep = "\t")
PGKB.db <- read.csv("relationships.tsv", header = T, sep = "\t")

setwd(paste(Root, ":/Projet Safir/Data SAFIR/Safir Output/Safir SupplGeneList", sep = ""))


Explore <- paste("Chr", CHR[1], "to", CHR[length(CHR)], sep = "")
FileName <- paste(CGH$Safir.Id, CGH$BarCode, CGH$analysisDate, Tech, Explore, Select, "SupplTable.html", sep = "_")

# Select what Chr to explore
stab <- stab[which(stab$chrom %in% CHR),]

# Select what segment values: means or medians
Values <- stab$seg.mean
if(use.medians) Values <- stab$seg.med

# Select what imbalances: Gain, Loss or Both
Select <- match.arg(Select)
switch(Select,  	Gain = (whichseg = which(Values >= Thresh.Gain)),
			Loss = (whichseg = which(Values <= Thresh.Loss)),
			Both = (whichseg = which(Values >= Thresh.Gain | Values <= Thresh.Loss)))


nseg = length(whichseg)
PB <- tkProgressBar(title = "The R'atWork BaBar", min = 0, max = nseg, width = 500)
k = 0

Append = F
if(length(whichseg)>1) Append <- T
for(i in whichseg){

	Sys.sleep(0.1)
	# launch & increment the pBar
	k = k + 1
	chr <- stab$chrom[i]
      setTkProgressBar(PB, k, label = paste("Segment #", i, " on Chr#", chr,": ", k," of ", nseg, " in progress...", sep = ""))

	cat(paste("Segment#", i, " on Chr", chr,": ", k," of ", nseg, sep = ""), "\n")
	start <- stab$loc.start[i]
	end <- stab$loc.end[i]
	seg.value <- round(Values[i], 3)
	index.full <- which(full$GenomicPos>=start & full$GenomicPos<=end)
	Symbols <- as.character(unique(full$Symbol[index.full]))
	if(any(Symbols == "")) Symbols <- Symbols[-which(Symbols == "")]
	if(length(Symbols)>0){
		cat("Symbols\n", Symbols, "\n")
			refseq <- c()
			for(j in 1:length(Symbols)) refseq <- c(refseq, as.character(full$RefSeq[which(full$Symbol == Symbols[j])][1]))

			# Search for PGKB annotations
			pgkb.genes <- ifelse(as.character(Symbols) %in% PGKB.db$Entity1_name, as.character(Symbols), "")
			cat("PgKB\n", ifelse(pgkb.genes!="", pgkb.genes, "---"), "\n")
			pgkb.Ids <- rep("", length(pgkb.genes))
			if(any(pgkb.genes != "")){
				for(pg in 1:length(pgkb.genes))
					if(pgkb.genes[pg] != ""){
						pgkb.Ids[pg] <- as.character(unique(PGKB.db$Entity1_id[which(PGKB.db$Entity1_name == pgkb.genes[pg])]))
						}
				}
			pgkb.Ids <- substr(pgkb.Ids, 6, 50)

			# Search for Census Cancer annotations
			is.census <- ifelse(as.character(Symbols) %in% census$Symb, "Yes", "")

			# Search for NCBI annotations
			if(Restrict){
				of.interest <- which(is.census == "Yes" | pgkb.Ids != "")
				Symbols <- as.character(Symbols[of.interest])
				refseq <- as.character(refseq[of.interest])
				is.census <- as.character(is.census[of.interest])
				pgkb.Ids <- as.character(pgkb.Ids[of.interest])
				cat(length(of.interest), "\n",
					ifelse(is.census == "Yes" & pgkb.Ids!="", paste(Symbols, "°", "+", sep = ""), ifelse(is.census == "Yes", paste(Symbols, "°", sep = ""), paste(Symbols, "+", sep = ""))),
					"\n")
				}

			if(length(Symbols)>0){
				tmp.request <- GeneRequest.v2(as.character(Symbols), hg19.info, verbose = FALSE)
				tmp.request <- cbind.data.frame(tmp.request, RefSeq = refseq, Cancer.Census.List = is.census, PgKB = pgkb.Ids)
				ord <- order(tmp.request$Chr.start)
				tmp.request <- tmp.request[ord,]
				tmp.request <- rbind.data.frame(tmp.request, rep("", ncol(tmp.request)))
				GeneId <- list(	tmp.request$GeneId,
							tmp.request$Symb, 
							tmp.request$PgKB, 
							ifelse(tmp.request$Cancer.Census.List == "Yes", as.character(tmp.request$Symb), ""))
				htmlpage(GeneId, filename = FileName, title = paste("Segment on Chr", chr, ": ", start, "-", end, ", Log2R = ", seg.value, ", nb genes = ", length(Symbols),sep = ""),
						othernames = tmp.request[,c(2:3, 6:8, 11:13)],
						table.head = c("EntrezGene", "Kegg", "PharmgKB", "CensusCancer", colnames(tmp.request[,c(2:3, 6:8, 11:13)])), 
						repository = list("en", "kegg", "pgkb", "census"), append = Append)
				}
		}
	}
	close(PB)
	rm(census, full, PKGB.db)
	cat("Supplementary table saved: ", paste(getwd(), FileName, sep = "/"), "\n")
}

