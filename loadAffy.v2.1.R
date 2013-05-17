
loadAffy.v2.1 <- function(hg19.info){	# Safir.Id, 

require(limma)
require(DNAcopy)
require(preprocessCore)
require(AMORE)
require(affy)
require(robustbase)
require(mclust)
require(tcltk)

op <- par(no.readonly = T)

# plot.new()

##########################################################################################

 fileName <- tclvalue(tkgetOpenFile()) # Very simple, isn't it?
 if (!nchar(fileName)) {
     tkmessageBox(message = "No file was selected!")
 } else {
     tkmessageBox(message = paste("The file selected was", fileName))
}

# head(fileName)	# return le filename as string

fnames <- unlist(strsplit(fileName, "/"))
nnames <- length(fnames)

path <- fnames[1]
for(i in 2:(nnames-1))
	path <- paste(path, fnames[i], sep = "/")
Folder <- fnames[nnames-1]
FileName <- fnames[nnames]

setwd(paste(path, "",sep = "/"))

##########################################################################################

	# Step1: Chargement des donnees
cat("Step1.a: Collecte des informations mArray:")
a.info <- try(read.csv(paste(FileName, sep = ""), header = F, fill = T, skip = 0, nrows = 1000, sep = "="), silent = T)

	{
	if(class(a.info)!="try-error"){
		ArrayRef = which(a.info[,1] == "#option-reference-input")	# "#reference-file"
		ArraySet = as.character(a.info[ArrayRef, 2])
		datadate <- which(a.info[,1] == "#state-time-start")		# "#state-time-start"
		ScanDate <- as.character(a.info[datadate, 2])
		ScanDate <- strsplit(ScanDate, " ")[[1]]
		m <- ScanDate[2]; d <- ScanDate[3]; y <- ScanDate[5]
		ScanDate = paste(m, d, y, sep = "-")

		# si presafir
		# ScanDate <- ArraySet <- NA
		
		a.info <- try(read.csv(paste(FileName, sep = ""), header = F, fill = T, skip = 0, nrows = 1000, sep = "\t"), silent = T)
		datastart <- which(a.info == "ProbeSet")
		cat("\nCGH info have been correctly downloaded:", "\n\tFolder: ", getwd(), "\n\tFileName:", FileName, "\n")
		}
	else stop("\t",FileName, ": No such file on directory (glup!). Please fix it!\n")
	}

#
#
#	ScanDate = as.character(a.info[188, 1])
#	ScanDate = strsplit(ScanDate, "=")[[1]][2]
#	ScanDate = strsplit(ScanDate, " ")[[1]]
#	d <- ScanDate[3]; m <- ScanDate[2]; y <- ScanDate[5]
#	ScanDate = paste(m, d, y, sep = "-")
#
	Safir.Id <- strsplit(FileName, "_")[[1]][1]
	BarCode <- strsplit(FileName, "_")[[1]][3]
#
	LabId = Folder
#
	analysisDate <- format(Sys.Date(), "%m-%d-%Y")
#
cat("\n\tSafir.Id:", Safir.Id, "\n\tBarCode:", BarCode, "\n\tArraySet:", ArraySet, "\n\tScanDate:", ScanDate, "\n\tLabId:", LabId, "\n\tAnalysis Date:", analysisDate, "\n\t")

cat("\nStep1.b: Chargement des donnees mArray:")
cgh <- read.csv(paste(FileName, sep = ""), header = T, skip = datastart - 1, sep = "\t")
cgh <- cgh[which(substr(cgh$ProbeSet, 1, 2) == "CN"),]
cat("\nCGH data have been correctly downloaded\n")

##########################################################################################

	# Step4: reannot selon hg19
cat("\nStep4: Calcul des positions genomiques selon hg19\n")

	# reannot genomic Pos selon hg19

levels(cgh$Chromosome)[23] <- 23
levels(cgh$Chromosome)[24] <- 24
cgh$Chromosome <- as.numeric(as.character(cgh$Chromosome))

cgh <- cgh[order(cgh$Chromosome, cgh$Chromosomal.Position),]

cum.len <- cumsum(hg19.info$length)

ch.len <- rep(0, length(which(cgh$Chromosome==1)))
for (i in 2:24) ch.len <- c(ch.len, rep(cum.len[i-1], length(which(cgh$Chromosome==i))))
cgh.gen.pos <- cgh$Chromosomal.Position + ch.len*1000

cgh <- cbind.data.frame(cgh[,1:3], GenomicPos = cgh.gen.pos, SmoothSignal = cgh[,6], cgh[,4:5])
colnames(cgh)[2:3] <- c("ChrNum", "ChrStart")
cgh <- cgh[order(cgh$ChrNum, cgh$ChrStart),]

	# Centrage du Y
lr <- cgh$Log2Ratio
chr <- cgh$ChrNum
lr[which(chr==24)] <- lr[which(chr==24)] - median(lr[which(chr==24)])
cgh$Log2Ratio <- lr

cat("\n\n\n\t and the winner is....:p\n\n\n")
return(list(cgh = cgh, Safir.Id = Safir.Id, BarCode = BarCode, ArraySet = ArraySet, ScanDate = ScanDate, LabId = LabId, analysisDate = analysisDate))
}
