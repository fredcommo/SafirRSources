## htmlTab: script pour la génération de tableau html.
## USAGE: 	prend en entrée le fichier targets issus de build.target() (GC5fit.R) et un fichier d'annotation de gène.
##			Le fichier d'annotation est issu du script d'annotation de Frédéric COMMO GeneRequest.R
## AUTHOR: clement mazoyer. clement.mazoyer@igr.fr
##
## VERSION:
# v2.1: mise à jour pour lire le nouveau fichier target de GC5 v5.10
# v2.1: abandon de htmlpage() (library(annotate)) pour R2HTML, plus flexible

##TODO: html_css: faire un systeme de liste avec les différentes balises en sous liste.

#require(annotate)
#require(XML)
require(R2HTML)

mnt.proj<-"/mnt/proj/"

htmlRapport<- function(targets, genelist, chr=1:23){


	if(file.exists(targets))
	{
		target<-read.table(targets, h=T, stringsAsFactors=F, sep="\t");	## targetfile
	}else{
		stop("target file not found");
	}

	if(file.exists(genelist))
	{
		genes<-read.table(genelist, h=T, stringsAsFactors=F, sep="\t");
	}else{
		stop("genelist file not found")
	}

		cs <- chrload("hs", "19")
		
	
	###############
	## VERSION 2 ##
	###############
	
	##############################
	## Liste des liens partiels ##
	##############################
	linkModel<-NULL
	linkModel<-list(entrez="http://www.ncbi.nlm.nih.gov/sites/entrez?db=gene&cmd=Retrieve&dopt=Graphics&list_uids=",
	pgkb="http://www.pharmgkb.org/gene/",
	census="http://www.sanger.ac.uk/perl/genetics/CGP/cosmic?action=gene&ln=",
	kegg="http://www.genome.jp/dbget-bin/www_bget?hsa:",
	ctd="http://ctdbase.org/detail.go?type=gene&acc=")
	link_name<-names(linkModel)

	PGKB.db <- read.csv("~/script/MOSCATO/pgkb/relationships.tsv", header = T, sep = "\t", stringsAsFactors=F)
	pgkb_subset<-unique(PGKB.db[is.element(PGKB.db$Entity1_name, genes$Symb),c(1,2)])	#subset des annotations pgkb correspondant à la genelist
	pgkb_subset<-pgkb_subset[order(pgkb_subset[,2]),]
	
	## pour chaque échantillon
	for(f in 1:nrow(target))
	{
		samplename<- target$Test.barcode[f]
		patient<- target$Add.title[f]
		cat("Reading", samplename,"\n")
		nrf<-target$Normal.range.factor[f]

		##lecture des .gcx
		if(length(dir(pattern="bz2$"))!=0)	system("bunzip2 *bz2")
		x<-paste(samplename,"*","gcx",sep="")

		gcx_all<-read.table(system2(command="ls", args=x, stdout=T), h=T, stringsAsFactors=F)

		##lecture des .cbs NoCut
		x<-paste(samplename,"*NoCut*","cbs",sep="")
		cbs<-read.table(system2(command="ls", args=x, stdout=T), h=T, stringsAsFactors=F)
		cbs$genstart<- cbs$Start + cs$lchrsum[cbs$Chr]
		cbs$genend<- cbs$End + cs$lchrsum[cbs$Chr]

		## Ajout des logratio
		genelist<-genes[order(genes$Chr.start),]
		genelist$Log2Ratio<-sapply(genelist$Genom.start, FUN=function(x){cbs$Log2Ratio[max(which(x>cbs$genstart))]})

		## BUILDING THE L2R OBJECT
		l2rX <- gcx_all$LogRatio
		l2rX.mad <- median(abs(diff(l2rX)))
		current <- list(l2r=l2rX, l2r.mad=l2rX.mad)
		current_frame<-as.data.frame(current)
		scut <- nrf*current$l2r.mad

		## Récupération fes identifiants PGKB
		split<-paste("echo ",pgkb_subset[,1], " |cut -d: -f2")
		split<-sapply(split,function(x){system(x, intern=T)})
		pgkb_subset[,1]<-split
		
		## Regroupement des information de pgkb avec la liste de gènes d'intêret MOSCATO.
		## Garde que les gène qui sont en gain/perte. Élimination du "normal".
		genelist<-genelist[order(genelist[,2]),]
		genelist[which(genelist$Symb %in% intersect(genelist$Symb,pgkb_subset$Entity1_name)),"pgkb"]<-pgkb_subset$Entity1_id[which(pgkb_subset$Entity1_name %in% intersect(pgkb_subset$Entity1_name,genelist$Symb))]
		genelist[-which(genelist$Symb %in% intersect(genelist$Symb,pgkb_subset$Entity1_name)),"pgkb"]<-""
		genelist<-genelist[which(genelist$Log2Ratio> scut | genelist$Log2Ratio < -scut),]
		genelist<-genelist[order(genelist$Chr),]
		
		# Si le profil n'est pas plat
		if(nrow(genelist)!=0)
		{
			############################
			## Initialisation du html ##
			############################

			html_init(init=T,filename=patient,filedir=patient, titlename=paste("Patient",patient, sep=" "))
			filename<-paste(patient,"/",patient,".html",sep="")
			
			## inclusion du profil.png
			img<-system2(command="ls",args=paste(samplename,"*pang.png",sep=""),stdout=T)
			if(length(img)!=0)
			{
				system(paste("cp ", img," ", patient,"/",  sep=""))
				html_img(filename,img)
			}
			
			# Pour chaque liste (Moscato, Safir)
			for(l in unique(genelist$liste))
			{
				genelist_subset<-genelist[which(genelist$liste==l),]
				genelist_chr<-genelist_subset

				## Création des liens complets
				linkListe<-linkModel
				linkListe$entrez<-http_link(linkListe$entrez, genelist_chr$GeneId)
				linkListe$pgkb<-http_link(linkListe$pgkb, genelist_chr$pgkb)
				linkListe$census<-http_link(linkListe$census, genelist_chr$Symb)
				linkListe$kegg<-http_link(linkListe$kegg, genelist_chr$Symb)
				linkListe$ctd<-http_link(linkListe$ctd, genelist_chr$GeneId)
				tab_name<-c(link_name,"chromosome","description","essai","log2R","Status")
				
				
				if(nrow(genelist_chr)>1)
				{
					## Création des liens html
					linkListe<-sapply(linkListe,function(x){paste("<a href=",x," target=_blank >",genelist_chr$Symb,"</a>",sep="")})
					linkListe<-as.data.frame(linkListe)
					
				}else{	#corrige un bug si il n'y a qu'une seule ligne
					linkListe<-as.data.frame(t(paste("<a href=",linkListe," target=_blank >",genelist_chr$Symb,"</a>",sep="")))
				}
				
				
				if(!is.null(nrow(linkListe)))
				{
					linkListe<-cbind(linkListe, genelist_chr$Chr, genelist_chr$Name, genelist_chr$essai,round(genelist_chr$Log2Ratio,2),NA)
					colnames(linkListe)<-cbind(tab_name)
					linkListe$Status<-ifelse(linkListe$log2R>=scut,"<span id=\"gain\"> G </span>","<span id=\"loss\"> L </span>")
				}else{
					linkListe<-c(linkListe, genelist_chr$Chr,genelist_chr$Name, genelist_chr$essai,round(genelist_chr$Log2Ratio,2),NA)
					linkListe<-t(linkListe)
					colnames(linkListe)<-c(tab_name)
					linkListe$Status<-ifelse(linkListe$log2R>=scut,"<span id=\"loss\"> G </span>","<span id=\"gain\"> L </span>")
				}
				#coloration des logRatio
				linkListe$log2R<-sapply(linkListe$log2R,function(x)
				{
					if(x>=1) x<-paste("<span id=\"gain\">",x,"</span>",sep="")
					else if(x<=-1) x<-paste("<span id=\"loss\">",x,"</span>",sep="")
					else x
				})
				
				if(linkListe$essai=="")	linkListe<-linkListe[,-which(colnames(linkListe)=="essai")]
				linkListe<-linkListe[order(linkListe[,which(colnames(linkListe)=="chromosome"),]),]	#ordonne par chromosome
				html_tab(filename=patient,data=linkListe,caption=paste("Liste",l,sep=" "))
			}
		}
		html_css(filename=filename, css_liste=T)
		html_init(end=T, filename=patient, filedir=patient, titlename=patient)	#termine le .html
	}	# fin pour chaque sample
}


## Construction des liens complets. Recoit la requête type (http) et les id pour compléter la requete.
http_link<-function(http,id)
{
	http<-sapply(http,function(x){paste(x,id,sep="")})
	return(http)
}

## initialisation ou fermeture d'un fichier html.
#fichier dans le repertoire filename.
#nom du fichier: filename.
#titre du fichier: titlename
html_init<-function(init=F,end=F,filename,filedir,titlename=NULL)
{
	path<-paste(filedir,"/",filename,".html",sep="")

	if(init)
	{
		if(!file.exists(filedir)) system(paste("mkdir",filedir,sep=" "))
		HTMLInitFile(filedir, filename=filename, Title=titlename)
		HTML("<meta http-equiv=\"content-type\" content=\"text/html; charset=utf-8\" />",file=path)
		HTML(as.title(titlename),file=path)
	}
	if(end)	HTMLEndFile(file=path)
	#<hr size=1>
	#<font size=-1>
	#Generated on: <i>Thu Feb 16 11:26:30 2012</i> - <b>R2HTML</b> 
	#<hr size=1>
	#</body>
	#</html>
}

## insertion d'une image dans un fichier html filename déja initialisé.
html_img<-function(filename, img_path)
{
	HTMLInsertGraph(GraphFileName=img_path, Width="50%", file=filename)
}

## insertion de code css dans le html filename déja initialisé.
html_css<-function(filename, css_liste)
{
	write.table("<style type=\"text/css\">",file=filename, quote=F, append=T,col.names=F, row.names=F)
	write.table("table{width:100%;}",file=filename, quote=F, append=T,col.names=F, row.names=F)
	write.table("tr.firstline{background-color:#FFA500;}",file=filename, quote=F, append=T,col.names=F, row.names=F)
	write.table("a:link{text-decoration:none;color:blue;}",file=filename, quote=F, append=T,col.names=F, row.names=F)
	write.table("a:visited{text-decoration:none;color:#FFA500;}",file=filename, quote=F, append=T,col.names=F, row.names=F)
	write.table("a:hover{text-decoration:underline;color:red;}",file=filename, quote=F, append=T,col.names=F, row.names=F)
	write.table("h2{background-color:#FFA500;text-align:center;}",file=filename, quote=F, append=T,col.names=F, row.names=F)
	write.table("span{font-weight:bold;}",file=filename, quote=F, append=T,col.names=F, row.names=F)
	write.table("#gain{color:blue;}",file=filename, quote=F, append=T,col.names=F, row.names=F)
	write.table("#loss{color:red;}",file=filename, quote=F, append=T,col.names=F, row.names=F)
	write.table("</style> ",file=filename, quote=F, append=T,col.names=F, row.names=F)
}

## Ajout de données contenues dans data à une page html existante.
##filename: nom du fichier.	Argument obligatoires!
#data: données à inclure dans le fichier
#caption: titre du tableau de données (valable si data est un dataframe/matrix)

html_tab<-function(filename,data=NULL,caption=NULL)
{
	path<-paste(filename,"/",filename,".html",sep="")
	if(!is.null(data)){	HTML(data, row.names=F, innerBorder=1, caption=caption,captionalign="top", file=path)}
}



# IMPORT CHR DATA
chrload <- function(sp="hs", gb="19")
{
  
	sptxt <- sp
	if (sp == "hs") sptxt <- "hg"
	gv <- paste(sptxt, gb, sep="")
  
	## IMPORT BANDES CYTO ET PREPARATION POUR PLOT KARYO
	cat("Importing chromosomes data  ...\n")
	if (gv == "hg18") cytob <- read.table(paste(mnt.proj, "cgh/hg18/cytoBandIdeo.hg18", sep=""), sep="\t", header=T, stringsAsFactors=F, comment.char="_")
	else if (gv == "hg19") cytob <- read.table(paste(mnt.proj, "cgh/hg19/cytoBandIdeo.hg19", sep=""), sep="\t", header=T, stringsAsFactors=F, comment.char="_")
	else if (gv == "mm9") cytob <- read.table(paste(mnt.proj, "cgh/mm9/cytoBandIdeo.mm9", sep=""), sep="\t", header=T, stringsAsFactors=F, comment.char="_")
	
	cytob$chrA <- substr(cytob$X.chrom, 4, nchar(cytob$X.chrom))
	cytob$chr <- cytob$chrA
	if (sp == "hs")
	{
		cytob$chr[which(cytob$chr == "X")] <- 23
		cytob$chr[which(cytob$chr == "Y")] <- 24
	} else if (sp == "mm") {
		cytob$chr[which(cytob$chr == "X")] <- 20
		cytob$chr[which(cytob$chr == "Y")] <- 21
	}
  
	cytob$chr <- as.numeric(cytob$chr)
	cytob <- cytob[order(cytob$chr),]

	lchrxx <- sapply(unique(cytob$chr), function(x) { max((cytob[which(cytob$chr == x),])$chromEnd) })
	lchrtoadd <- c(0, lchrxx[1:length(lchrxx)-1])
	glen <- sum(as.numeric(lchrxx))
	lchrsum <- sapply(c(1:length(lchrtoadd)), function(x) { sum(lchrtoadd[1:x]) })

  
	return(list(cytob=cytob, lchrxx=lchrxx, lchrsum=lchrsum, lchrtoadd=lchrtoadd, glen=glen, cur.glen=glen))
}









