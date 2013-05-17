## genelist: fichier issu du script d'annotation de Frédéric COMMO GeneRequest.R

require(annotate)
require(XML)

mnt.proj<-"/mnt/proj/"

htmlTab<- function(targets, genelist, chr=1:23){


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

		cs <<- chrload("hs", "19")
		

# build a PGKB repository
	repofun <- function(ids, ...)
	paste("http://www.pharmgkb.org/gene/", ids, sep = "")
	setRepository("pgkb", repofun,verbose=F)

# build a Sanger Census Cancer repository
	repofun <- function(ids, ...)
	paste("http://www.sanger.ac.uk/perl/genetics/CGP/cosmic?action=gene&ln=", ids, sep = "")
	setRepository("census", repofun,verbose=F)

# build a Kegg repository
	repofun <- function(ids, ...)
	paste("http://www.genome.jp/dbget-bin/www_bget?hsa:", ids, sep = "")
	setRepository("kegg", repofun,verbose=F)


#build a CTD repository
	repofun<-function(ids,...)
	paste("http://ctdbase.org/detail.go?type=gene&acc=",ids,sep="")
	setRepository("ctd",repofun,verbose=F)



	PGKB.db <- read.csv("~/script/MOSCATO/pgkb/relationships.tsv", header = T, sep = "\t", stringsAsFactors=F)
	#pgkb.genes <- ifelse(as.character(genelist$Symb) %in% PGKB.db$Entity1_name, as.character(genelist$Symb), "")
	pgkb_subset<-unique(PGKB.db[is.element(PGKB.db$Entity1_name, genes$Symb),c(1,2)])
	pgkb_subset<-pgkb_subset[order(pgkb_subset[,2]),]
	
	# pour chaque échantillon
	for(f in 1:nrow(target))
	{
		samplename<- target$SampleName[f]
		patient<- target$Add.title[f]
		cat("Reading", samplename,"\n")
		nrf<-target$Normal.range.factor[f]
		
		##lecture des .gcx
		if(length(dir(pattern="bz2$"))!=0)	system("bunzip2 *bz2")
		x<-paste(samplename,"*","gcx",sep="")

		gcx_all<-read.table(system2(command="ls", args=x, stdout=T), h=T, stringsAsFactors=F)
		
		##lecture des .cbs NoCut
		x<-paste(samplename,"*","cbs",sep="")
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

		#pgkb_subset<-unique(PGKB.db[is.element(PGKB.db$Entity1_name, genelist$Symb),c(1,2)])
		#pgkb_subset<-pgkb_subset[order(pgkb_subset[,2]),]
		
		split<-paste("echo ",pgkb_subset[,1], " |cut -d: -f2")
		split<-sapply(split,function(x){system(x, intern=T)})
		pgkb_subset[,1]<-split
		
		genelist<-genelist[order(genelist[,2]),]
		
		genelist[which(genelist$Symb %in% intersect(genelist$Symb,pgkb_subset$Entity1_name)),"pgkb"]<-pgkb_subset$Entity1_id[which(pgkb_subset$Entity1_name %in% intersect(pgkb_subset$Entity1_name,genelist$Symb))]
		genelist[-which(genelist$Symb %in% intersect(genelist$Symb,pgkb_subset$Entity1_name)),"pgkb"]<-""
		
		genelist<-genelist[which(genelist$Log2Ratio> scut | genelist$Log2Ratio < -scut),]
		genelist<-genelist[order(genelist$Chr),]

		if(file.exists(paste(patient,"_",samplename,".html",sep="")))
		{
			app<-F
		}else{
			app<-T
		}
		## Pour chaque liste (Moscato, Safir)
		for(l in unique(genelist$liste))
		{
			genelist_subset<-genelist[which(genelist$liste==l),]
			
			## Pour chaque chromosome: si on veut une table par chromosome. Moins lisible.
#			for(chr in unique(genelist_subset$Chr))
#			{
				genelist_chr<-genelist_subset#[which(genelist_subset$Chr==chr),]

				GeneId <- list(genelist_chr$GeneId, genelist_chr$Symb, genelist_chr$pgkb, genelist_chr$GeneId)

				filename<-paste(patient,"_",samplename,".html",sep="")
				write.table("<meta http-equiv=\"content-type\" content=\"text/html; charset=utf-8\" />",filename,quote=F, row.names=F, col.names=F,append=T)
				htmlpage(GeneId, filename=filename,title=paste("<b>Patient ",patient,"- Liste de gènes", l, "</b>"), 
				table.head=c("EntrezGene","Kegg","pgkb","CTD","Symbole","Chromosome","Description","Essai","log2(ratio)"),
				repository = list("en","kegg","pgkb", "ctd"),
				othernames=list(genelist_chr$Symb, genelist_chr$Chr, genelist_chr$Name, genelist_chr$essai, genelist_chr$Log2Ratio), append = app, center.table=T)
				write.table("<br>",filename,quote=F, row.names=F, col.names=F,append=T)
				app<-T
#			}
			write.table("<hr>",filename,quote=F, row.names=F, col.names=F,append=T)
		}
		## Design du html
		html_tmp<-system(paste("sed 's/<TABLE/<TABLE width=90%/g'",filename,"|sed 's/<TH>/<TH bgcolor=#FFA500>/g'", sep=" "), intern=T)
		
		#<meta http-equiv="content-type" content="text/html; charset=utf-8" />
		#html_tmp<-system(paste("sed 's/<TABLE/<TABLE width=1400/g'",filename, "|sed 's/<TD>/<TD width=175>/g' |sed 's/<TD align=/<TD width=175 align=/g'|sed 's/<TH>/<TH bgcolor=#C2A2A>/g'",   sep=" "), intern=T)
		write.table(html_tmp,filename, quote=F, row.names=F, col.names=F)
	}
}


# IMPORT CHR DATA
chrload <- function(sp="hs", gb="19") {
  
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

	cytob$x1 <- 0.15
	cytob$x2 <- 0.85
	cytob$x1[which((cytob$gieStain == "acen") | (cytob$gieStain == "gvar") | (cytob$gieStain == "stalk"))] <- 0.25
	cytob$x2[which((cytob$gieStain == "acen") | (cytob$gieStain == "gvar") | (cytob$gieStain == "stalk"))] <- 0.75
	cytob$gieStain[which(cytob$gieStain == "gneg")] <- "white"
	cytob$gieStain[which(cytob$gieStain == "gpos25")] <- "grey75"
	cytob$gieStain[which(cytob$gieStain == "gpos33")] <- "grey66"
	cytob$gieStain[which(cytob$gieStain == "gpos50")] <- "grey50"
	cytob$gieStain[which(cytob$gieStain == "gpos66")] <- "grey33"
	cytob$gieStain[which(cytob$gieStain == "gpos75")] <- "grey25"
	cytob$gieStain[which(cytob$gieStain == "gpos100")] <- "black"
	cytob$gieStain[which(cytob$gieStain == "acen")] <- "yellow"
	cytob$gieStain[which(cytob$gieStain == "gvar")] <- "darkred"
	cytob$gieStain[which(cytob$gieStain == "stalk")] <- "darkolivegreen"
  
	return(list(cytob=cytob, lchrxx=lchrxx, lchrsum=lchrsum, lchrtoadd=lchrtoadd, glen=glen, cur.glen=glen))
}









