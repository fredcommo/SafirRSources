EditSupTab.v5 <- function(stab = seg.table, CGH = Scgh, Tech = "", use.medians = TRUE, Select = c("Gain", "Loss", "Both"), Thresh.Gain = 1, Thresh.Loss = -1, CHR = 1:23, Restrict = FALSE, Root = "E"){

require(annotate)
require(XML)
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

# build a KeggToDrug repository
	repofun <- function(ids, ...)
	paste("http://www.genome.jp/kegg-bin/get_htext?htext=br08303_target.keg&query=", ids, sep = "")
	setRepository("keggD", repofun)

gsearch <- function (geneSymb, database){

	# geneSymb : use official symbols. Multiple requests are accepted, e.g. "EGFR, Homo sapiens"
	# database : let's have a look at 'http://eutils.ncbi.nlm.nih.gov/entrez/query/static/eutils_help.html' for details
	# 		on available databases and other e-tools as well.
	# ! This function can return more than one Id !

	gsrch.stem <- "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
	gsrch.mode <- paste("db=", database, "&retmode=xml","&term=", sep = "")
	doc <- xmlTreeParse(paste(gsrch.stem, gsrch.mode, geneSymb, sep = ""), isURL = TRUE, useInternalNodes = TRUE)
	sapply(c("//Id"), xpathApply, doc = doc, fun = xmlValue)
	}

gsummary <- function (id, database){

	# id is provided by gsearch()
	# database : let's have a look at 'http://eutils.ncbi.nlm.nih.gov/entrez/query/static/eutils_help.html' for details
	# 		on available databases and other e-tools as well.

	sum.stem <- "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?"
	sum.mode <- paste("db=", database, "&id=", sep = "")
	doc <- xmlTreeParse(paste(sum.stem, sum.mode, id, sep = ""), isURL = TRUE, useInternalNodes = TRUE)
	sapply(c("//Item"), xpathApply, doc = doc, fun = xmlValue)
	}

# KeggToDrug function
KeggToDrug <- function(GeneSymb){
	srch.stem <- "http://www.genome.jp/kegg-bin/get_htext?"
	srch.mode <- "htext=br08303_target.keg&query="
	loc <- paste(srch.stem, srch.mode, GeneSymb, sep="")
	doc <- htmlParse(loc)
	tabNodes <- getNodeSet(doc, "//table")
	tabHtml <- readHTMLTable(tabNodes[[3]])
	drug.list <- as.character(tabHtml[-1,2])
	ndrug <- length(which(!is.na(drug.list) & drug.list != ""))
	return(ndrug)
	}


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

			# search for KeggToDrugs
			KG <- tkProgressBar(title = "Searching for Kegg's info", min = 0, max = length(Symbols), width = 500)
			nDrugs <- c()
			for(kg in 1:length(Symbols)){
				Sys.sleep(0.1)
				setTkProgressBar(KG, kg , label = paste(Symbols[kg], " in progress...", kg, " of", length(Symbols), " genes",sep = ""))
				tmp.id <- gsearch(paste(Symbols[kg], "homo sapiens"), "gene")
				tmp.id <- unlist(tmp.id)
				tmp.name = "NA"
				count = 1
				while(tmp.name != Symbols[kg] & count <= length(tmp.id)){
					id <- g.id[count]
					gsum <- unlist(gsummary(id, "gene"))
					tmp.name <- gsum[1]
					count = count + 1
					}
				ndrug = 0
				if(tmp.name == GeneSymbol) ndrug <- KeggToDrug(id)
				nDrugs <- c(nDrugs, ndrug)
				}
			close(KG)

			# Search for NCBI annotations
			if(Restrict){
				of.interest <- which(is.census == "Yes" | pgkb.Ids != "" | nDrugs != 0)
				Symbols <- as.character(Symbols[of.interest])
				refseq <- as.character(refseq[of.interest])
				is.census <- as.character(is.census[of.interest])
				pgkb.Ids <- as.character(pgkb.Ids[of.interest])
				nDrugs <- nDrugs[of.interest]
				cat(length(of.interest), "\n")
				}

			if(length(Symbols)>0){
				tmp.request <- GeneRequest.v2(as.character(Symbols), hg19.info, verbose = FALSE)
				tmp.request <- cbind.data.frame(	tmp.request,
										RefSeq = refseq,
										PgKB = pgkb.Ids,
										Cancer.Census.List = ifelse(is.census != "yes", as.character(Symbols), ""),
										KeggToDrugs = ifelse(nDrugs != 0, as.character(Symbols), "_"))

				ord <- order(tmp.request$Chr.start)
				tmp.request <- tmp.request[ord,]
				tmp.request <- rbind.data.frame(tmp.request, rep("", ncol(tmp.request)))
				GeneId <- list(	tmp.request$GeneId,
							tmp.request$Symb, 
							tmp.request$PgKB, 
							tmp.request$Cancer.Census.List,
							tmp.request$KeggToDrugs
							)

				htmlpage(GeneId, filename = FileName, title = paste("Segment on Chr", chr, ": ", start, "-", end, ", Log2R = ", seg.value, ", nb genes = ", length(Symbols),sep = ""),
						othernames = tmp.request[,c(2:3, 6:8, 11:13)],
						table.head = c("EntrezGene", "Kegg", "PharmgKB", "CensusCancer", "KeggToDrugs", colnames(tmp.request[,c(2:3, 6:8, 11:13)])), 
						repository = list("en", "kegg", "pgkb", "census", "keggD"), append = Append)
				}
		}
	}
	close(PB)
	rm(census, full, PKGB.db)
	cat("Supplementary table saved: \n\t", paste(getwd(), FileName, sep = "/"), "\n")
}

