### REMARQUES:
# FC 16/02/12:
#	- Ajout de CTDbase.
#	- Modif de la mise en page.
#	- Thresh.gain remplacé par ThreshGain, Thresh.loss remplacé par ThreshLoss
# ajout des test systeme: pour les setdir() ce n'est pas pareil si on est sur windows ou linux
#ajout d'un setwd() à la fin pour revenir au repertoir de travail si on doit re-executer le script.
#bug ligne 233



EditSupTab.v7.5 <- function(stab = seg.table, CGH = Scgh, ReqFunc = Request, use.medians = TRUE, Select = c("Gain", "Loss", "Both"), ThreshGain = 1, ThreshLoss = -1, CHR = 1:23, Restrict = FALSE, Root = "E"){

#require(annotate)
#require(XML)
require(tcltk)
##VERSION 7##
require(R2HTML)

	# Check for system
	system <- Sys.info()["sysname"]
	workingDir<-getwd()

	##############################
	## Liste des liens partiels ##
	##############################
	linkModel <- NULL
	linkModel <- list(Entrez = "http://www.ncbi.nlm.nih.gov/sites/entrez?db=gene&cmd=Retrieve&dopt=Graphics&list_uids=",	# Link to EntrezGene
				Kegg = "http://www.genome.jp/dbget-bin/www_bget?hsa:",														# Link to Kegg database
				SangerCensus = "http://www.sanger.ac.uk/perl/genetics/CGP/cosmic?action=gene&ln=",							# Link to Sanger Census Cancer
				PharmgKB = "http://www.pharmgkb.org/gene/",																	# Link to Pharmacogenomics Knowledge base
				KeggToDrugs = "http://www.genome.jp/kegg-bin/get_htext?htext=br08303_target.keg&query=",					# Link to Kegg database 'genes to drugs'
				CTDbase = "http://ctdbase.org/detail.go?type=gene&acc=",													# Link to CTdatabase (genes/drugs information, direct or indirect interactions)
				ClinicalTrials = "http://clinicaltrials.gov/ct2/results?term=")												# Link to ClinicalTrials.gov (official site for clinical trials registry)
	linkNames <- names(linkModel)

################################################################





################################################################


	####################
	## Bloc principal ##
	####################

	# Build path regarding system used (quite different for linux and windows)
	if(system=="Linux"){
		setwd(paste(Root, "/Projet Safir/Data SAFIR/Arrays Infos", sep = ""))
		FileDir <- paste(Root, "/Projet Safir/Data SAFIR/Safir Output/Safir SupplGeneList", sep = "")
	}
	if(system=="Windows"){
		setwd(paste(Root, ":/Projet Safir/Data SAFIR/Arrays Infos", sep = ""))
		FileDir <- paste(Root, ":/Projet Safir/Data SAFIR/Safir Output/Safir SupplGeneList", sep = "")
	}

	# load annotation tables
	census <- read.table("CensusTable_Annot_FC.txt", header = T, sep = "\t")
	full <- read.table("022060_4x180K_hg19_20110628_GenomePos&GC_GNames_FC.txt", header = T, sep = "\t")
	PGKB.db <- read.csv("relationships.tsv", header = T, sep = "\t")
	KeggToDrug.db <- read.csv("Kegg_Drugs_Table.txt", header = T, sep = "\t")
	CT.db <- read.csv("ClinicalTrials_Drugs_Table.txt", header = T, sep = "\t")
	CTDb <- read.csv("CTD_chem_gene_ixns_SansDuplic.txt", header = T, sep = "\t")
	CTDbList <- CTDb$GeneSymbol

	setwd(FileDir)

	# Select what Chr to explore
	stab <- stab[which(stab$chrom %in% CHR),]

	# Select what segment values: means or medians
	Values <- stab$seg.mean
	if(use.medians) Values <- stab$seg.med

	# What imbalances to consider: Gain, Loss or Both, and define both values and table title
	Select <- match.arg(Select)
	switch(Select,  Gain = (whichseg = which(Values >= ThreshGain)),
					Loss = (whichseg = which(Values <= ThreshLoss)),
					Both = (whichseg = which(Values >= ThreshGain | Values <= ThreshLoss))
					)
	gainSubTitl <- paste("gained (>", round(ThreshGain, 3), ")", sep = "")
	lossSubTitl <- paste("lost (<", round(ThreshLoss, 3), ")", sep = "")
	switch(Select,  Gain = (selecType = gainSubTitl),
					Loss = (selecType = lossSubTitl),
					Both = (selecType = paste(gainSubTitl, "or", lossSubTitl))
					)
					
	# Define File name and main table title
	Explore <- paste("Chr", CHR[1], "to", CHR[length(CHR)], sep = "")
	FileName <- paste(CGH$Safir.Id, CGH$BarCode, CGH$analysisDate, CGH$Tech, Explore, Select, "SupplTable", sep = "_")
	TitleName <- paste("Supplementary tables: ", CGH$Safir.Id, " / ", CGH$BarCode, "<br>",
						CGH$Tech, " / ", CGH$analysisDate, " / ", paste("Chr", CHR[1], "to", CHR[length(CHR)]), ": ", paste(selecType, "segments"))

	# Define progressBar
	nseg = length(whichseg)
	PB <- tkProgressBar(title = "The R'atWork BaBar", min = 0, max = nseg, width = 500)
	k = 0


	## Initialisation du html		
	html_init(init = T, filedir = FileDir, filename = FileName, titlename = TitleName)	# changer titlename = tabName
	
	# filename <- paste(FileName, "/", FileName,".html", sep="")
	# filename <- paste(FileName, ".html", sep = "")


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
			pgkb.genes <- ifelse(as.character(Symbols) %in% PGKB.db$Entity1_name, as.character(Symbols), "-")
			cat("PgKB\n", ifelse(pgkb.genes!="-", pgkb.genes, "---"), "\n")
			pgkb.Ids <- rep("-", length(pgkb.genes))
			if(any(pgkb.genes != "-")){
				for(pg in 1:length(pgkb.genes))
					if(pgkb.genes[pg] != "-"){
						tmp <- as.character(unique(PGKB.db$Entity1_id[which(PGKB.db$Entity1_name == pgkb.genes[pg])]))
						pgkb.Ids[pg] <- substr(tmp, 6, 50)
					}
			}

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
									else return (NA)})

			# Search for CTDbase
			is.ctd <- ifelse(as.character(Symbols) %in% CTDbList, Symbols, NA)
			
			# search for ClinicalTrials
			is.CT <- ifelse(as.character(Symbols) %in% CT.db$Symb, Symbols, NA)
			CT.Symb <- apply(as.data.frame(is.CT), 1, function(x){
									if(as.character(x) %in% Symbols){
										CT.index <- which(CT.db$Symb == x)
										if(CT.db$CTfound[CT.index]!=0) return (x)
										else return(NA)
									}
									else return (NA)})
			

			# Search for restrictions: if Restrict = TRUE
			if(Restrict){
				of.interest <- which(!is.na(is.census) | pgkb.Ids != "-" | !is.na(Kegg.GI) | !is.na(CT.Symb) | !is.na(is.ctd))
				Symbols <- as.character(Symbols[of.interest])
				refseq <- as.character(refseq[of.interest])
				is.census <- as.character(is.census[of.interest])
				pgkb.Ids <- as.character(pgkb.Ids[of.interest])
				Kegg.GI <- Kegg.GI[of.interest]
				is.ctd <- as.character(is.ctd[of.interest])
				CT.Symb <- CT.Symb[of.interest]
				cat(length(of.interest), "\n")
			}

			# Building output table
			if(length(Symbols)>0){
				# Search for NCBI annotations
				tmp.request <- ReqFunc(as.character(Symbols), hg19.info, verbose = FALSE)
				tmp.request <- cbind.data.frame(tmp.request,
												RefSeq = refseq,
												Cancer.Census.List = ifelse(!is.na(is.census), as.character(is.census), "-"),
												PgKB = pgkb.Ids,
												KeggToDrugs = ifelse(!is.na(Kegg.GI), as.character(Kegg.GI), "-"),
												CTDbase = ifelse(!is.na(is.ctd), as.character(tmp.request$GeneId), "-"),
												ClinicalTrials = ifelse(!is.na(CT.Symb), as.character(CT.Symb), "-")
												)

				ord <- order(tmp.request$Chr.start)
				tmp.request <- tmp.request[ord,]
				
				# Build a list of Ids to complete http links
				GeneId <- list(	tmp.request$GeneId,						# GeneId is used for EntrezGene
							tmp.request$Symb, 							# Symbol is used for Kegg
							tmp.request$Cancer.Census.List,				# Symbol is used for Sanger census list
							tmp.request$PgKB, 							# PharmgKB Id is used for PharmgKB
							tmp.request$KeggToDrugs,					# GeneId is used for Kegg to drug
							tmp.request$CTDbase,						# GeneId is used for CTDbase
							tmp.request$ClinicalTrials)					# Symbols is used for CliniclaTrials.gov

				# Edition html
				######################################
				##	Impression des data dans html	##
				######################################
				
				## Building complete links
				linkListe <- linkModel
				linkListe$Entrez <- http_link(linkListe$Entrez, tmp.request$GeneId)
				linkListe$Kegg <- http_link(linkListe$Kegg, tmp.request$GeneId)
				linkListe$SangerCensus <- http_link(linkListe$SangerCensus, tmp.request$Cancer.Census.List)
				linkListe$PharmgKB <- http_link(linkListe$PharmgKB, tmp.request$PgKB)
				linkListe$KeggToDrugs <- http_link(linkListe$KeggToDrugs, tmp.request$KeggToDrugs)
				linkListe$CTDbase <- http_link(linkListe$CTDbase, tmp.request$CTDbase)
				linkListe$ClinicalTrials <- http_link(linkListe$ClinicalTrials, tmp.request$ClinicalTrials)

				## Unlist linkList to build subtables of links
				# linkListe <- sapply(linkListe, function(x){paste("<a href=", x, " target=_blank >", tmp.request$Symb, "</a>", sep="")})
				linkList2 <- c()
				for(L in 1:length(linkListe)){
					tmplist <- sapply(as.data.frame(linkListe[[L]]), function(x){paste("<a href=", x, " target=_blank >", GeneId[[L]], "</a>", sep="")})
					linkList2 <- cbind(linkList2, tmplist)
				}
			
				linkList2 <- as.data.frame(linkList2)
				colnames(linkList2) <- linkNames

				# Liste des noms dans la table
				#  names(tmp.request)
				# [1] "Query"			"Symb"			"Name"				"Org"				"Chr"               
				# [6] "Cytoband"		"Chr.start"		"Chr.end"			"Genom.start"		"Genom.end"         
				# [11] "RangeGB.Id"		"GeneId"		"RefSeq"			"PgKB"				"Cancer.Census.List"
				# [16] "KeggToDrugs"	"CTDbase"		"ClinicalTrials"    

				# tab_name <- c(linkNames, colnames(tmp.request[,c(2:3, 6, 12:13)]))	##LA VARIABLE linkNames était remplacée par link_name
				
				# bug de formatage si il n'y a qu'un gÃšne sur le fragment...
				if(ncol(linkList2)>1){
					linkList2 <- cbind.data.frame(tmp.request[,c(2, 3, 6, 12, 13)], linkList2)			# Information first(see above for #col concordance), then links.
				}
				else{
						linkList2 <- cbind.data.frame(tmp.request[,c(2, 3, 6, 12, 13)], t(linkList2))		# Information first(see above for #col concordance), then links.
				}
					
					html_tab(filename = FileName, data = linkList2, caption = paste("<b>Segment on Chr", chr, ": ", round(start/1e6, 3), "-", round(end/1e6, 3), "(Mb), Log2R = ", seg.value, ", nb genes = ", length(Symbols),"</b>", sep = ""))
			}
		}
	}
	close(PB)
	rm(census, full, PGKB.db, KeggToDrug.db, CT.db, CTDb, CTDbList)
	cat("Supplementary table saved: \n\t", paste(getwd(), FileName, sep = "/"), "\n")
	
	## 
	html_css(filename = FileName, css_liste = T)
	html_init(end = T, filedir = FileDir, filename = FileName)
	setwd(workingDir)
}












	####################
	## Fonctions html ##
	####################


## Construction des liens complets. Recoit la requÃªte type (http) et les id pour complÃ©ter la requete.
http_link <- function(http, id)
{
	http <- sapply(http, function(x){paste(x, id, sep = "")})
	return(http)
}

## initialisation ou fermeture d'un fichier html.
# filedir : 	path
# filename: 	filename.
# titlename: 	main table title
html_init <- function(init = F, end = F, filedir, filename, titlename = NULL)
{
	path <- paste(filedir, "/", filename, ".html", sep="")

	if(init)
	{
		if(!file.exists(filedir)) system(paste("mkdir",filedir,sep=" "))
		HTMLInitFile(filedir, filename = filename, Title = titlename)
		HTML("<meta http-equiv=\"content-type\" content=\"text/html; charset=utf-8\" />",file=path)
		HTML(as.title(titlename),file=path)
	}
	if(end)	HTMLEndFile(file=path)
}

## insertion d'une image dans un fichier html filename déja initialisé.
html_img <- function(filename, img_path)
{
	HTMLInsertGraph(GraphFileName=img_path, Width="50%", file=filename)
}

## insertion de code css dans le html filename dÃ©ja initialisÃ©.
html_css <- function(filename, css_liste)
{
	filename <-  paste(filename, ".html", sep = "")
	write.table("<style type=\"text/css\">", file = filename, quote=F, append=T, col.names=F, row.names=F)
	write.table("table{width:100%;}", file = filename, quote=F, append=T,col.names=F, row.names=F)
	write.table("tr.firstline{background-color:#FFBD9D;}", file = filename, quote=F, append=T,col.names=F, row.names=F)				# fond entêtes de sous-tables
	write.table("a:link{text-decoration:none;color:blue;}", file = filename, quote=F, append=T,col.names=F, row.names=F)			# couleur du lien
	write.table("a:visited{text-decoration:none;color:#8A008A;}", file = filename, quote=F, append=T,col.names=F, row.names=F)		# couleur du lien après activation
	write.table("a:hover{text-decoration:underline;color:red;}", file = filename, quote=F, append=T,col.names=F, row.names=F)		# couleur du lien au passage de souris
	write.table("h2{background-color:#FFA366;text-align:center;}", file = filename, quote=F, append=T,col.names=F, row.names=F)		# fond entête principal
	write.table("</style> ",file = filename, quote=F, append=T,col.names=F, row.names=F)
}

## Ajout de donnÃ©es contenues dans data Ã  une page html existante.
##filename: nom du fichier.	Argument obligatoires!
#data: donnÃ©es Ã  inclure dans le fichier
#caption: titre du tableau de donnÃ©es (valable si data est un dataframe/matrix)

html_tab <- function(filename, data=NULL, caption=NULL)
{
	# path <- paste(filename, "/", filename, ".html", sep="")
	path <- paste(filename, ".html", sep = "")
	if(!is.null(data)){	HTML(data, row.names = F, innerBorder = 1, caption = caption, captionalign = "top", file = path)}
}