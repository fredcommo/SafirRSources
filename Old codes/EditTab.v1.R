### REMARQUES:
# FC 07/03/12:
# Ajout d'options pour le choix des databases supplémentaires

# FC 16/02/12:
#	- Ajout de CTDbase.
#	- Modif de la mise en page.
#	- Thresh.gain remplacé par ThreshGain, Thresh.loss remplacé par ThreshLoss
# ajout des test systeme: pour les setdir() ce n'est pas pareil si on est sur windows ou linux
#ajout d'un setwd() à la fin pour revenir au repertoir de travail si on doit re-executer le script.
#bug ligne 233



EditTab.v1 <- function(myTab = greq, Root = "E"){

require(annotate)
require(XML)
require(tcltk)
##VERSION 7.x##
require(R2HTML)

	# Check for system
	system <- Sys.info()["sysname"]
	workingDir <- getwd()

	##############################
	## Liste des liens partiels ##
	##############################
	linkModel <- NULL
	linkModel$GeneId = "http://www.ncbi.nlm.nih.gov/sites/entrez?db=gene&cmd=Retrieve&dopt=Graphics&list_uids="					# Link to EntrezGene
#	linkModel$Kegg = "http://www.genome.jp/dbget-bin/www_bget?hsa:"																# Link to Kegg database
#	if(sangercensus) linkModel$SangerCensus = "http://www.sanger.ac.uk/perl/genetics/CGP/cosmic?action=gene&ln="				# Link to Sanger Census Cancer
#	if(pharmgkb) linkModel$PharmgKB = "http://www.pharmgkb.org/gene/"															# Link to Pharmacogenomics Knowledge base
#	if(keggtodrug) linkModel$KeggToDrugs = "http://www.genome.jp/kegg-bin/get_htext?htext=br08303_target.keg&query="			# Link to Kegg database 'genes to drugs'
#	if(ctdbase) linkModel$CTDbase = "http://ctdbase.org/detail.go?type=gene&acc="												# Link to CTdatabase (genes/drugs information, direct or indirect interactions)
#	if(clintrials) linkModel$ClinicalTrials = "http://clinicaltrials.gov/ct2/results?term="										# Link to ClinicalTrials.gov (official site for clinical trials registry)
	
#	linkModel <- as.list(linkModel)								
#	linkNames <- names(linkModel)

################################################################





################################################################


	####################
	## Bloc principal ##
	####################

	# Build path regarding system used (quite different for linux and windows)
	if(system=="Linux"){
		setwd(paste(Root, "/Projet Safir/Data SAFIR/Arrays Infos", sep = ""))
		FileDir <- paste(Root, "/Projet Safir/Data SAFIR/Safir Output/Safir GenesOfInterest", sep = "")
	}
	if(system=="Windows"){
		setwd(paste(Root, ":/Projet Safir/Data SAFIR/Arrays Infos", sep = ""))
		FileDir <- paste(Root, ":/Projet Safir/Data SAFIR/Safir Output/Safir GenesOfInterest", sep = "")
	}

	# load annotation tables
	full <- read.table("022060_4x180K_hg19_20110628_GenomePos&GC_GNames_FC.txt", header = T, sep = "\t")
#	if(sangercensus) census <- read.table("CensusTable_Annot_FC.txt", header = T, sep = "\t")
#	if(pharmgkb) PGKB.db <- read.csv("relationships.tsv", header = T, sep = "\t")
#	if(keggtodrug) KeggToDrug.db <- read.csv("Kegg_Drugs_Table.txt", header = T, sep = "\t")
#	if(ctdbase){
#		CTDb <- read.csv("CTD_chem_gene_ixns_SansDuplic.txt", header = T, sep = "\t")
#		CTDbList <- CTDb$GeneSymbol
#	}	
#	if(clintrials) CT.db <- read.csv("ClinicalTrials_Drugs_Table.txt", header = T, sep = "\t")

	
	setwd(FileDir)

#	# Select what Chr to explore
#	stab <- stab[which(stab$chrom %in% CHR),]
#
#	# Select what segment values: means or medians
#	Values <- stab$seg.mean
#	if(use.medians) Values <- stab$seg.med
#
#	# What imbalances to consider: Gain, Loss or Both, and define both values and table title
#	Select <- match.arg(Select)
#	switch(Select,  Gain = (whichseg = which(Values >= ThreshGain)),
#					Loss = (whichseg = which(Values <= ThreshLoss)),
#					Both = (whichseg = which(Values >= ThreshGain | Values <= ThreshLoss))
#					)
#	gainSubTitl <- paste("gained (>", round(ThreshGain, 3), ")", sep = "")
#	lossSubTitl <- paste("lost (<", round(ThreshLoss, 3), ")", sep = "")
#	switch(Select,  Gain = (selecType = gainSubTitl),
#					Loss = (selecType = lossSubTitl),
#					Both = (selecType = paste(gainSubTitl, "or", lossSubTitl))
#					)
#					
#	# Define File name and main table title
#	Explore <- paste("Chr", CHR[1], "to", CHR[length(CHR)], sep = "")
	FileName <- paste(Scgh$Safir.Id, Scgh$BarCode, Scgh$analysisDate, Scgh$Tech, "GenesOfIntHTML", sep = "_")
	TitleName <- paste("Genes of interest: ", Scgh$Safir.Id, " / ", Scgh$BarCode, "<br>")
#						CGH$Tech, " / ", CGH$analysisDate)
#


	## Initialisation du html		
	html_init(init = T, filedir = FileDir, filename = FileName, titlename = TitleName)	# changer titlename = tabName
	Append <- T
	
	

	# Edition html
	######################################
	##	Impression des data dans html	##
	######################################
				
	## Building complete links
	linkListe <- linkModel
	linkListe$GeneId <- http_link(linkListe$GeneId, as.character(myTab$GeneId))
	linkList2 <- sapply(as.data.frame(linkListe$GeneId), function(x){paste("<a href=", x, " target=_blank >", as.character(myTab$GeneId), "</a>", sep="")})
	linkList2 <- as.data.frame(linkList2)
	colnames(linkList2) <- "GeneId"
 

		# tab_name <- c(linkNames, colnames(tmp.request[,c(2:3, 6, 12:13)]))	##LA VARIABLE linkNames était remplacée par link_name
	linkList2 <- cbind.data.frame(myTab[,3:4], linkList2, myTab[,c(6:9, 14:17, 20:22)])		# Information first(see above for #col concordance), then links.
	html_tab(filename = FileName, data = linkList2)
	cat("Genes of Interest table saved: \n\t", paste(getwd(), FileName, sep = "/"), "\n")
	
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
	write.table("table{width:90%;}", file = filename, quote=F, append=T, col.names=F, row.names=F)
	write.table("tr.firstline{background-color:#FFBD9D;}", file = filename, quote=F, append=T, col.names=F, row.names=F)				# fond entêtes de sous-tables
	write.table("a:link{text-decoration:none;color:blue;}", file = filename, quote=F, append=T, col.names=F, row.names=F)				# couleur du lien
	write.table("a:visited{text-decoration:none;color:#8A008A;}", file = filename, quote=F, append=T, col.names=F, row.names=F)			# couleur du lien après activation
	write.table("a:hover{text-decoration:underline;color:red;}", file = filename, quote=F, append=T, col.names=F, row.names=F)			# couleur du lien au passage de souris
	write.table("h2{background-color:#FFA366;text-align:center;}", file = filename, quote=F, append=T, col.names=F, row.names=F)		# fond entête principal
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