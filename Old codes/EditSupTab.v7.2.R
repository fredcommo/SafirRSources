### REMARQUES:
# ajout des test systeme: pour les setdir() ce n'est pas pareil si on est sur windows ou linux
#ajout d'un setwd() à la fin pour revenir au repertoir de travail si on doit re-executer le script.
#bug ligne 233



EditSupTab.v7.2 <- function(stab = seg.table, CGH = Scgh, use.medians = TRUE, Select = c("Gain", "Loss", "Both"), Thresh.Gain = 1, Thresh.Loss = -1, CHR = 1:23, Restrict = FALSE, Root = "E"){

#require(annotate)
#require(XML)
require(tcltk)
##VERSION 7##
require(R2HTML)


system<-Sys.info()["sysname"]
workingDir<-getwd()

	##############################
	## Liste des liens partiels ##
	##############################
	linkModel <- NULL
	linkModel <- list(Entrez = "http://www.ncbi.nlm.nih.gov/sites/entrez?db=gene&cmd=Retrieve&dopt=Graphics&list_uids=",
				Kegg = "http://www.genome.jp/dbget-bin/www_bget?hsa:",
				SangerCensus = "http://www.sanger.ac.uk/perl/genetics/CGP/cosmic?action=gene&ln=",
				PharmgKB = "http://www.pharmgkb.org/gene/",
				KeggToDrugs = "http://www.genome.jp/kegg-bin/get_htext?htext=br08303_target.keg&query=",
				ClinicalTrials = "http://clinicaltrials.gov/ct2/results?term=")
	linkNames <- names(linkModel)

################################################################





################################################################


	####################
	## Bloc principal ##
	####################

if(system=="Linux")
{

	setwd(paste(Root, "/Projet Safir/Data SAFIR/Arrays Infos", sep = ""))
	FileDir <- paste(Root, "/Projet Safir/Data SAFIR/Safir Output/Safir SupplGeneList", sep = "")
}
if(system=="Windows")
{
	setwd(paste(Root, ":/Projet Safir/Data SAFIR/Arrays Infos", sep = ""))
	FileDir <- paste(Root, ":/Projet Safir/Data SAFIR/Safir Output/Safir SupplGeneList", sep = "")
}

census <- read.table("CensusTable_Annot_FC.txt", header = T, sep = "\t")
full <- read.table("022060_4x180K_hg19_20110628_GenomePos&GC_GNames_FC.txt", header = T, sep = "\t")
PGKB.db <- read.csv("relationships.tsv", header = T, sep = "\t")
KeggToDrug.db <- read.csv("Kegg_Drugs_Table.txt", header = T, sep = "\t")
CT.db <- read.csv("ClinicalTrials_Drugs_Table.txt", header = T, sep = "\t")


setwd(FileDir)


Explore <- paste("Chr", CHR[1], "to", CHR[length(CHR)], sep = "")
# FileName <- paste(CGH$Safir.Id, CGH$BarCode, CGH$analysisDate, Tech, Explore, Select, "SupplTable.html", sep = "_")
FileName <- paste(CGH$Safir.Id, CGH$BarCode, CGH$analysisDate, CGH$Tech, Explore, Select, "SupplTable", sep = "_")





# Select what Chr to explore
stab <- stab[which(stab$chrom %in% CHR),]

# Select what segment values: means or medians
Values <- stab$seg.mean
if(use.medians) Values <- stab$seg.med

# Select what imbalances: Gain, Loss or Both
Select <- match.arg(Select)
switch(Select,  	Gain = (whichseg = which(Values >= Thresh.Gain)),
			Loss = (whichseg = which(Values <= Thresh.Loss)),
			Both = (whichseg = which(Values >= Thresh.Gain | Values <= Thresh.Loss))
			)


nseg = length(whichseg)
PB <- tkProgressBar(title = "The R'atWork BaBar", min = 0, max = nseg, width = 500)
k = 0



	############################
	## Initialisation du html ##
	############################

html_init(init = T, filedir = FileDir, filename = FileName, titlename = FileName)	# changer titlename = tabName
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
				of.interest <- which(!is.na(is.census) | pgkb.Ids != "-" | !is.na(Kegg.GI) | !is.na(CT.Symb))
				Symbols <- as.character(Symbols[of.interest])
				refseq <- as.character(refseq[of.interest])
				is.census <- as.character(is.census[of.interest])
				pgkb.Ids <- as.character(pgkb.Ids[of.interest])
				Kegg.GI <- Kegg.GI[of.interest]
				CT.Symb <- CT.Symb[of.interest]
				cat(length(of.interest), "\n")
				}

			# Building output table
			if(length(Symbols)>0)
			{
				tmp.request <- GeneRequest.v3(as.character(Symbols), hg19.info, verbose = FALSE)
				tmp.request <- cbind.data.frame(	tmp.request,
										RefSeq = refseq,
										Cancer.Census.List = ifelse(!is.na(is.census), as.character(is.census), "-"),
										PgKB = pgkb.Ids,
										KeggToDrugs = ifelse(!is.na(Kegg.GI), as.character(Kegg.GI), "-"),
										ClinicalTrials = ifelse(!is.na(CT.Symb), as.character(CT.Symb), "-")
										)

				ord <- order(tmp.request$Chr.start)
				tmp.request <- tmp.request[ord,]
				# tmp.request <- rbind.data.frame(tmp.request, rep("", ncol(tmp.request)))
				GeneId <- list(	tmp.request$GeneId,
							tmp.request$Symb, 
							tmp.request$Cancer.Census.List,
							tmp.request$PgKB, 
							tmp.request$KeggToDrugs,
							tmp.request$ClinicalTrials)

				# Edition html
				######################################
				##	Impression des data dans html	##
				######################################
				
				## CrÃ©ation des adresses completes
				linkListe <- linkModel
				linkListe$Entrez <- http_link(linkListe$Entrez, tmp.request$GeneId)
				linkListe$Kegg <- http_link(linkListe$Kegg, tmp.request$GeneId)
				linkListe$SangerCensus <- http_link(linkListe$SangerCensus, tmp.request$Cancer.Census.List)
				linkListe$PharmgKB <- http_link(linkListe$PharmgKB, tmp.request$PgKB)
				linkListe$KeggToDrugs <- http_link(linkListe$KeggToDrugs, tmp.request$KeggToDrugs)
				linkListe$ClinicalTrials <- http_link(linkListe$ClinicalTrials, tmp.request$ClinicalTrials)

				## CrÃ©ation des liens html
				# linkListe <- sapply(linkListe, function(x){paste("<a href=", x, " target=_blank >", tmp.request$Symb, "</a>", sep="")})
				linkList2 <- c()
				for(L in 1:length(linkListe)){
					tmplist <- sapply(as.data.frame(linkListe[[L]]), function(x){paste("<a href=", x, " target=_blank >", GeneId[[L]], "</a>", sep="")})
					linkList2 <- cbind(linkList2, tmplist)
					}
				# rm(linkListe)
				linkList2 <- as.data.frame(linkList2)
				colnames(linkList2) <- linkNames

				# Liste des noms dans la table
				#  names(tmp.request)
				# [1] "Query"              "Symb"               "Name"               "Org"                "Chr"               
				# [6] "Cytoband"           "Chr.start"          "Chr.end"            "Genom.start"        "Genom.end"         
				# [11] "RangeGB.Id"         "GeneId"             "RefSeq"             "PgKB"               "Cancer.Census.List"
				# [16] "KeggToDrugs"        "ClinicalTrials"    

				tab_name <- c(linkNames, colnames(tmp.request[,c(2:3, 6:8, 12:13)]))	##LA VARIABLE linkNames était remplacée par link_name
				
				# bug de formatage si il n'y a qu'un gÃšne sur le fragment...
				if(ncol(linkList2)>1)
				{
					linkList2 <- cbind.data.frame(linkList2, tmp.request[,c(2:3, 6:8, 12:13)])
				}else{
					linkList2 <- cbind.data.frame(t(linkList2), tmp.request[,c(2:3, 6:8, 12:13)])
				}
				
				html_tab(filename=FileName, data=linkList2, caption=paste("<b>Segment on Chr", chr, ": ", start, "-", end, ", Log2R = ", seg.value, ", nb genes = ", length(Symbols),"</b>", sep = ""))
			}
		}
	}
	close(PB)
	rm(census, full, PGKB.db, KeggToDrug.db, CT.db)
	cat("Supplementary table saved: \n\t", paste(getwd(), FileName, sep = "/"), "\n")
	
	##
	html_css(filename = FileName, css_liste = T)
	html_init(end=T, filedir = FileDir, filename = FileName)
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
#fichier dans le repertoire filename.
#nom du fichier: filename.
#titre du fichier: titlename
html_init <- function(init=F, end=F, filedir, filename, titlename=NULL)
{
	path <- paste(filedir, "/", filename, ".html", sep="")

	if(init)
	{
		if(!file.exists(filedir)) system(paste("mkdir",filedir,sep=" "))
		HTMLInitFile(filedir, filename=filename, Title=titlename)
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
	write.table("tr.firstline{background-color:#FFBD9D;}", file = filename, quote=F, append=T,col.names=F, row.names=F)			# fond entêtes de sous-tables
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
	path <- paste(filename, ".html", sep="")
	if(!is.null(data)){	HTML(data, row.names=F, innerBorder=1, caption=caption,captionalign="top", file=path)}
}