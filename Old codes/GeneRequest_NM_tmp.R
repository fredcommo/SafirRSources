#################################################################################################

# Exemple liste SAFIR

GeneRequest_NM_tmp <- function(genelist, hg19.info, DB = "unigene", verbose = TRUE){

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
		TITLE <- CHROMOSOME <- ORGANISM <- GENE <- CLUSTERID <- TAXID <- SEQ_COUNT <- EST_COUNT <- GENEID  <- RECORDTYPE  <- NA			# in case of error only !

		gname <- genelist[i]
		# gname <- "CD44"
		ids <- unlist(gsearch(paste(gname, "homo sapiens"), database = DB))					#homo sapiens
		if(is.null(ids)) cat("\n***", gname, "not found: ArétéConRy ! ***\n\n")

		if(!is.null(ids)){ 
			if(verbose) cat(gname, "found:", length(ids), "ids.\n")

			j = 1
			Id <- ids[j]
			gsum <- unlist(gsummary(paste(Id, "homo sapiens"), database = DB))					#homo sapiens
    			TITLE <- gsum[1]
			CHROMOSOME <- gsum[2]
			ORGANISM <- gsum[3]
    			GENE <- gsum[4]
			CLUSTERID  <- gsum[5]
			TAXID <- gsum[6]
			SEQ_COUNT <- gsum[7]
			EST_COUNT <- gsum[8]
			GENEID <- gsum[9]
			RECORDTYPE <- gsum[10]

			if(length(ids)>1){
				while (ORGANISM != "Homo sapiens" | GENE != gname & j < length(ids)){
					j = j + 1
					Id <- ids[j]
					gsum <- unlist(gsummary(paste(Id, "homo sapiens"), database = DB))					#homo sapiens
    					TITLE <- gsum[1]
					CHROMOSOME <- gsum[2]
					ORGANISM <- gsum[3]
    					GENE <- gsum[4]
					CLUSTERID  <- gsum[5]
					TAXID <- gsum[6]
					SEQ_COUNT <- gsum[7]
					EST_COUNT <- gsum[8]
					GENEID <- gsum[9]
					RECORDTYPE <- gsum[10]
					}
					# cat(gname, "related to Homo sapiens found at occurence:", j, "\n")
				}
			}

		gene.annot <- rbind(gene.annot, c(gname, GENE, TITLE, CHROMOSOME, ORGANISM, CLUSTERID , TAXID, SEQ_COUNT, EST_COUNT, GENEID, RECORDTYPE))
		}	

	gene.annot <- as.data.frame(gene.annot)
	colnames(gene.annot) <- c("Query", "Symb", "TITLE", "CHROMOSOME", "ORGANISM", "CLUSTERID", "TAXID", "SEQ_COUNT", "EST_COUNT", "GENEID", "RECORDTYPE")
	close(pb)
	return(gene.annot)
}
