#################################################################################################

# Exemple liste SAFIR

GeneRequest.v2 <- function(genelist, hg19.info, verbose = TRUE){

require(XML)
require(tcltk)

# Fonctions requêtes

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

# End fonctions requêtes

	Total = length(genelist)
	pb <- tkProgressBar(title = "The R'atWork BaBar", min = 0, max = Total, width = 500)

	cum.len <- cumsum(hg19.info$length)

	gene.annot <- c()
	for(i in 1:Total){

		Sys.sleep(0.1)

			# launch & increment the pBar
	    setTkProgressBar(pb, i, label = paste("Take it easy... I'm workin' for U... :p ", round(i/Total*100, 2), "% done"))

		Symb <- "Not found"
		Name <- Org <- Chr <- Cytoband <- Chr.start <- Chr.end <- RefSeq <- NA			# in case of error only !

		gname <- genelist[i]
		ids <- unlist(gsearch(paste(gname, "homo sapiens"), "gene"))					#homo sapiens
		if(is.null(ids)) cat("\n***", gname, "not found: ArétéConRy ! ***\n\n")

		if(!is.null(ids)){ 
			if(verbose) cat(gname, "found:", length(ids), "ids.\n")

			j = 1
			Id <- ids[j]
			gsum <- unlist(gsummary(paste(Id, "homo sapiens"), "gene"))					#homo sapiens
			Org <- gsum[3]
			status <- gsum[2]
			verifId <- as.numeric(gsum[5])
			Chr <- gsum[6]

			if(length(ids)>1){
				while ((Org != "Homo sapiens" | status=="reserved" | verifId != 0 | Chr =="") & j < length(ids)){
					j = j + 1
					Id <- ids[j]
					gsum <- unlist(gsummary(paste(Id, "homo sapiens"), "gene"))			#homo sapiens
					status <- gsum[2]
					Org <- gsum[3]
					verifId <- as.numeric(gsum[5])
					Chr <- gsum[6]
					}
					# cat(gname, "related to Homo sapiens found at occurence:", j, "\n")
				}
			Symb <- gsum[1]
			Name <- gsum[2]
			Chr <- gsum[6]
			Cytoband <- gsum[8]
			if(Chr == "") Chr <- NA
			if(!is.na(Chr) & Chr == "X") Chr <- 23
			if(!is.na(Chr) & Chr == "Y") Chr <- 24
			if(!is.na(Chr) & Chr == "X, Y") Chr <- 23
			Chr <- as.numeric(Chr)
			RefSeq <- gsum[20]
			Chr.start <- as.numeric(gsum[21])
			Chr.end <- as.numeric(gsum[22])
			if(length(gsum)!=27){
				RefSeq <- gsum[19]
				Chr.start <- as.numeric(gsum[20])
				Chr.end <- as.numeric(gsum[21])
				}
			Genom.start <- Chr.start
			Genom.end <- Chr.end
			if(Chr>=2){
				Genom.start <- as.numeric(Chr.start + cum.len[Chr - 1]*1000)
				Genom.end <- as.numeric(Chr.end + cum.len[Chr - 1]*1000)
				}
			}

		gene.annot <- rbind(gene.annot, c(gname, Symb, Name, Org, Chr, Cytoband, Chr.start, Chr.end, Genom.start, Genom.end, RefSeq, Id))
		}	

	gene.annot <- as.data.frame(gene.annot)
	colnames(gene.annot) <- c("Query", "Symb", "Name", "Org", "Chr", "Cytoband", "Chr.start", "Chr.end", "Genom.start", "Genom.end", "RangeGB.Id", "GeneId")
	close(pb)
	return(gene.annot)
}
