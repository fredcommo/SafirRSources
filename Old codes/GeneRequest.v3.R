#################################################################################################

# Exemple liste SAFIR

GeneRequest.v3 <- function(genelist, hg19.info, DB = "gene", verbose = TRUE, PB = T){

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
	if(PB) pb <- tkProgressBar(title = "The R'atWork BaBar", min = 0, max = Total, width = 500)
		

	cum.len <- cumsum(hg19.info$length)

	gene.annot <- c()
	for(i in 1:Total){

		if(PB) Sys.sleep(0.1)

			# launch & increment the pBar
		if(PB) setTkProgressBar(pb, i, label = paste("Take it easy... I'm workin' for U... :p ", round(i/Total*100, 2), "% done"))

		Symb <- "Not found"
		Name <- Org <- Chr <- Cytoband <- Chr.start <- Chr.end <- Genom.start <- Genom.end <- RefSeq <- Id <- NA				# in case of error only !

		gname <- genelist[i]
		try.gsearch <- try(gsearch(paste(gname, "homo sapiens"), database = DB), silent = T)							#homo sapiens
		ntries = 0
		while(class(try.gsearch) == "try-error" & ntries < 10){
			try.gsearch <- try(gsearch(paste(gname, "homo sapiens"), database = DB), silent = T)						#homo sapiens
			ntries = ntries + 1
			}
		ids <- unlist(try.gsearch)																	#homo sapiens

		if(is.null(ids)) cat("\n***", gname, "not found: ArétéConRy ! ***\n\n")

		if(!is.null(ids)){ 
			if(verbose) cat(gname, "found:", length(ids), "ids.\n")

			j = 1
			Id <- ids[j]
			try.gsum <- try(gsummary(paste(Id, "homo sapiens"), database = DB), silent = T)
			ntries = 0
			while(class(try.gsum) == "try-error" & ntries < 10){
				try.gsum <- try(gsummary(paste(Id, "homo sapiens"), database = DB), silent = T)							#homo sapiens
				ntries = ntries + 1
			}
			gsum <- unlist(try.gsum)																	#homo sapiens

			Org <- gsum[3]
			status <- gsum[2]
			verifId <- as.numeric(gsum[5])
			Chr <- gsum[6]

			if(length(ids)>1){
				while ((Org != "Homo sapiens" | status=="reserved" | (!is.na(verifId) & verifId != 0) | (!is.na(Chr) & Chr == "")) & j < length(ids)){
					j = j + 1
					Id <- ids[j]
					try.gsum <- try(gsummary(paste(Id, "homo sapiens"), database = DB), silent = T)
					ntries = 0
					while(class(try.gsum) == "try-error" & ntries < 10){
						try.gsum <- try(gsummary(paste(Id, "homo sapiens"), database = DB), silent = T)					#homo sapiens
						ntries = ntries + 1
						}
					gsum <- unlist(try.gsum)
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
			if(!is.na(Chr) & Chr == "") Chr <- NA
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
			if(!is.na(Chr) & Chr>=2){
				Genom.start <- as.numeric(Chr.start + cum.len[Chr - 1]*1000)
				Genom.end <- as.numeric(Chr.end + cum.len[Chr - 1]*1000)
				}
			}

		gene.annot <- rbind(gene.annot, c(gname, Symb, Name, Org, Chr, Cytoband, Chr.start, Chr.end, Genom.start, Genom.end, RefSeq, Id))
		# rm(ids, try.gsum)
		}	

	gene.annot <- as.data.frame(gene.annot)
	colnames(gene.annot) <- c("Query", "Symb", "Name", "Org", "Chr", "Cytoband", "Chr.start", "Chr.end", "Genom.start", "Genom.end", "RangeGB.Id", "GeneId")
	if(PB) close(pb)
	return(gene.annot)
}
