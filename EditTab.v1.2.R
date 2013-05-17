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


EditTab.v1.2 <- function(myTab = greq, Root = Root){

	# Check for system
	system <- Sys.info()["sysname"]
	workingDir <- getwd()

	####################
	## Fonctions html ##
	####################


## Construction des liens complets. Recoit la requÃªte type (http) et les id pour complÃ©ter la requete.
	http_link <- function(http, id){
		http <- sapply(http, function(x){paste(x, id, sep = "")})
		return(http)
	}

## initialisation ou fermeture d'un fichier html.
# filedir : 	path
# filename: 	filename.
# titlename: 	main table title
	html_init <- function(init = F, end = F, filedir, filename, titlename = NULL){
		path <- paste(filedir, "/", filename, ".html", sep="")
		if(init){
			if(!file.exists(filedir)) system(paste("mkdir",filedir,sep=" "))
			HTMLInitFile(filedir, filename = filename, Title = titlename)
			HTML("<meta http-equiv=\"content-type\" content=\"text/html; charset=utf-8\" />",file=path)
			HTML(as.title(titlename),file=path)
		}
		if(end)	HTMLEndFile(file=path)
	}

## insertion d'une image dans un fichier html filename déja initialisé.
	html_img <- function(filename, img_path){
		HTMLInsertGraph(GraphFileName=img_path, Width="50%", file=filename)
	}

## insertion de code css dans le html filename dÃ©ja initialisÃ©.
	html_css <- function(filename, css_liste){
		filename <-  paste(filename, ".html", sep = "")
		write.table("<style type=\"text/css\">", file = filename, quote=F, append=T, col.names=F, row.names=F)
		write.table("table{width:100%;}", file = filename, quote=F, append=T, col.names=F, row.names=F)
		write.table("tr.firstline{background-color:#FFBD9D;}", file = filename, quote=F, append=T, col.names=F, row.names=F)				# fond entêtes de sous-tables
		write.table("a:link{text-decoration:none;color:blue;}", file = filename, quote=F, append=T, col.names=F, row.names=F)				# couleur du lien
		write.table("a:visited{text-decoration:none;color:#8A008A;}", file = filename, quote=F, append=T, col.names=F, row.names=F)			# couleur du lien après activation
		write.table("a:hover{text-decoration:underline;color:red;}", file = filename, quote=F, append=T, col.names=F, row.names=F)			# couleur du lien au passage de souris
		write.table("h2{background-color:#FFA366;text-align:center;}", file = filename, quote=F, append=T, col.names=F, row.names=F)		# fond entête principal
		write.table("span{font-weight:bold;}",file=filename, quote=F, append=T,col.names=F, row.names=F)
		write.table("#Norm{color:#A1A0A0;}",file=filename, quote=F, append=T,col.names=F, row.names=F)
		write.table("#Gain{color:#3075ED;}",file=filename, quote=F, append=T,col.names=F, row.names=F)
		write.table("#Ampli{color:#1100A6;}",file=filename, quote=F, append=T,col.names=F, row.names=F)	
		write.table("#Loss{color:#E50810;}",file=filename, quote=F, append=T,col.names=F, row.names=F)												# 
		write.table("</style> ",file = filename, quote=F, append=T,col.names=F, row.names=F)
	}

## Ajout de donnÃ©es contenues dans data Ã  une page html existante.
##filename: nom du fichier.	Argument obligatoires!
#data: donnÃ©es Ã  inclure dans le fichier
#caption: titre du tableau de donnÃ©es (valable si data est un dataframe/matrix)

	html_tab <- function(filename, data=NULL, caption=NULL){
		# path <- paste(filename, "/", filename, ".html", sep="")
		path <- paste(filename, ".html", sep = "")
		if(!is.null(data)){	HTML(data, row.names = F, innerBorder = 1, caption = caption, captionalign = "top", file = path)}
	}
# END HTML FUNCTIONS

	##############################
	## Liste des liens partiels ##
	##############################
	linkModel <- NULL
	linkModel$GeneId = "http://www.ncbi.nlm.nih.gov/sites/entrez?db=gene&cmd=Retrieve&dopt=Graphics&list_uids="					# Link to EntrezGene

################################################################


	####################
	## Bloc principal ##
	####################
	#On Mac: /Users/fredcommo/Documents/Projet Safir/Examples

	# Build path regarding system used (quite different for linux and windows)

	if(system=="Linux"){
		setwd(paste(Root, "/Projet Safir/Data SAFIR/Arrays Infos", sep = ""))
		FileDir <- paste(Root, "/Projet Safir/Data SAFIR/Safir Output/Safir GenesOfInterestHTML", sep = "")
	}
	if(system=="Windows"){
		setwd(paste(Root, ":/Projet Safir/Data SAFIR/Arrays Infos", sep = ""))
		FileDir <- paste(Root, ":/Projet Safir/Data SAFIR/Safir Output/Safir GenesOfInterestHTML", sep = "")
	}
	if(system=="Darwin"){
		setwd(paste(Root, "Arrays Infos", sep = ""))
		FileDir <- paste(Root, "Examples", sep = "")
	}

	# load annotation tables
	full <- read.table("022060_4x180K_hg19_20110628_GenomePos&GC_GNames_FC.txt", header = T, sep = "\t")
	
	setwd(FileDir)
					
	# Define File name and main table title
	FileName <- paste(Scgh$Safir.Id, Scgh$BarCode, Scgh$analysisDate, Scgh$Tech, "GenesOfIntHTML", sep = "_")
	TitleName <- paste("Safir01 Genes of Interest<br>Sample: ", Scgh$Safir.Id, " / ", Scgh$BarCode, "<br>")

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
	myTab[,20] <- round(myTab[,20], 3)
	myTab[,21] <- signif(as.real(myTab[,21]), 6)
	myTab[, c(10:11, 15:17)] <- round(myTab[, c(10:11, 15:17)]/1e3, 3)
	linkList2 <- cbind.data.frame(myTab[,3:4], linkList2, myTab[,6:7], myTab[, c(10:11, 15:17)], myTab[,20:22])		# Information first(see above for #col concordance), then links.
	linkList2$Status <- paste("<span id=", linkList2$Status, "", ">", linkList2$Status, "</span>", sep = "")
	colnames(linkList2)[c(4, 6:12)] <- c("Chr#", "GenomicStart Kb", "GenomicEnd Kb", "Segm.Start Kb", "Segm.End Kb", "Segm.Length Kb", "Log2R", "FoldChg")
	
	# linkList2$Status <- paste("<p ", linkList2$Status, " color:", "red", ";</p>", sep = "")
	html_tab(filename = FileName, data = linkList2)
	cat("Genes of Interest table saved: \n\t", paste(getwd(), FileName, sep = "/"), "\n")
	
	## 
	html_css(filename = FileName, css_liste = T)
	html_init(end = T, filedir = FileDir, filename = FileName)
	# return(linkList2$FoldChg)
	setwd(workingDir)
}
