EditSupTab.v6 <- function(stab = seg.table, CGH = Scgh, Tech = "", use.medians = TRUE, Select = c("Gain", "Loss", "Both"), Thresh.Gain = 1, Thresh.Loss = -1, CHR = 1:23, Restrict = FALSE, Root = "E"){

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

# build a ClinicalTrials repository
	repofun <- function(ids, ...)
	paste("http://clinicaltrials.gov/ct2/results?term=", ids, "+cancer&no_unk=Y", sep = "")
	setRepository("CT", repofun)


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


setwd(paste(Root, ":/Projet Safir/Data SAFIR/Arrays Infos", sep = ""))
census <- read.table("CensusTable_Annot_FC.txt", header = T, sep = "\t")
full <- read.table("022060_4x180K_hg19_20110628_GenomePos&GC_GNames_FC.txt", header = T, sep = "\t")
PGKB.db <- read.csv("relationships.tsv", header = T, sep = "\t")
KeggToDrug.db <- read.csv("Kegg_Drugs_Table.txt", header = T, sep = "\t")
CT.db <- read.csv("ClinicalTrials_Drugs_Table.txt", header = T, sep = "\t")

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

# htmlpage(GeneId, filename = FileName, title = paste("Segment on Chr", chr, ": ", start, "-", end, ", Log2R = ", seg.value, ", nb genes = ", length(Symbols),sep = ""),
#			othernames = tmp.request[,c(2:3, 6:8, 11:13)],
#			table.head = c("EntrezGene", "Kegg", "PharmgKB", "CensusCancer", "KeggToDrugs", "ClinicalTrials", colnames(tmp.request[,c(2:3, 6:8, 12:13)])), 
#			repository = list("en", "kegg", "pgkb", "census", "keggD", "CT"), append = Append)

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
			is.census <- ifelse(as.character(Symbols) %in% census$Symb, Symbols, NA)

			# search for KeggToDrugs
			is.Kegg <- ifelse(as.character(Symbols) %in% KeggToDrug.db$Symb, Symbols, NA)
			Kegg.GI <- apply(as.data.frame(is.Kegg), 1, function(x){
									if(as.character(x) %in% Symbols){
										kegg.index <- which(KeggToDrug.db$Symb == x)
										if(KeggToDrug.db$nDrug[kegg.index]!=0) return (KeggToDrug.db$GeneId[kegg.index])
										else return(NA)
										}
									else return (NA)
									}
									)


			# search for ClinicalTrials
			is.CT <- ifelse(as.character(Symbols) %in% CT.db$Symb, Symbols, NA)
			CT.Symb <- apply(as.data.frame(is.CT), 1, function(x){
									if(as.character(x) %in% Symbols){
										CT.index <- which(CT.db$Symb == x)
										if(CT.db$CTfound[CT.index]!=0) return (x)
										else return(NA)
										}
									else return (NA)
									}
									)
			

			# Search for NCBI annotations
			if(Restrict){
				of.interest <- which(!is.na(is.census) | pgkb.Ids != "" | !is.na(Kegg.GI) | !is.na(CT.Symb))
				Symbols <- as.character(Symbols[of.interest])
				refseq <- as.character(refseq[of.interest])
				is.census <- as.character(is.census[of.interest])
				pgkb.Ids <- as.character(pgkb.Ids[of.interest])
				Kegg.GI <- Kegg.GI[of.interest]
				CT.Symb <- CT.Symb[of.interest]
				cat(length(of.interest), "\n")
				}

			# Building output table
			if(length(Symbols)>0){
				tmp.request <- GeneRequest.v3(as.character(Symbols), hg19.info, verbose = FALSE)
				tmp.request <- cbind.data.frame(	tmp.request,
										RefSeq = refseq,
										PgKB = pgkb.Ids,
										Cancer.Census.List = ifelse(!is.na(is.census), as.character(is.census), "-"),
										KeggToDrugs = ifelse(!is.na(Kegg.GI), as.character(Kegg.GI), "-"),
										ClinicalTrials = ifelse(!is.na(CT.Symb), as.character(CT.Symb), "-")
										)

				ord <- order(tmp.request$Chr.start)
				tmp.request <- tmp.request[ord,]
				# tmp.request <- rbind.data.frame(tmp.request, rep("", ncol(tmp.request)))
				GeneId <- list(	tmp.request$GeneId,
							tmp.request$Symb, 
							tmp.request$PgKB, 
							tmp.request$Cancer.Census.List,
							tmp.request$KeggToDrugs,
							tmp.request$ClinicalTrials)

				# Edition html
				write.table("<meta http-equiv=\"content-type\" content=\"text/html bgcolor=#FFA500; charset=utf-8\" />", FileName, quote=F, row.names=F, col.names=F, append=T)
				htmlpage(GeneId, filename = FileName, title = paste("<b>Segment on Chr", chr, ": ", start, "-", end, ", Log2R = ", seg.value, ", nb genes = ", length(Symbols), "</b>", sep = ""),
						othernames = tmp.request[,c(2:3, 6:8, 12:13)],
						table.head = c("EntrezGene", "Kegg", "PharmgKB", "CensusCancer", "KeggToDrugs", "ClinicalTrials", colnames(tmp.request[,c(2:3, 6:8, 12:13)])), 
						repository = list("en", "kegg", "pgkb", "census", "keggD", "CT"), append = Append)
				write.table("<br>", FileName, quote=F, row.names=F, col.names=F, append=T)
				}

		## Design du html
		# html_tmp<-system(paste("sed 's/<TABLE/<TABLE width=90%/g'", FileName, "|sed 's/<TH>/<TH bgcolor=#FFA500>/g'", sep=" "), intern=T)
		
		#<meta http-equiv="content-type" content="text/html; charset=utf-8" />
		#html_tmp<-system(paste("sed 's/<TABLE/<TABLE width=1400/g'",filename, "|sed 's/<TD>/<TD width=175>/g' |sed 's/<TD align=/<TD width=175 align=/g'|sed 's/<TH>/<TH bgcolor=#C2A2A>/g'",   sep=" "), intern=T)
		# write.table(html_tmp, FileName, quote=F, row.names=F, col.names=F)

		}
	}
	close(PB)
	rm(census, full, PGKB.db, KeggToDrug.db, CT.db)
	cat("Supplementary table saved: \n\t", paste(getwd(), FileName, sep = "/"), "\n")
}

