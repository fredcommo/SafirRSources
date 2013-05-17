
loadAgilent.v4 <- function(Flags = TRUE, Root = "E"){

# The function activates an interactive selection of file to load
# Collect microarray informations and Cy5/Cy3 intensity values,
# and suppress flags and duplicated probes.
# Root : A letter indicating the hard disk partition to use. May run with other indications (~/home)

# Load required packages
require(limma)
require(DNAcopy)
require(preprocessCore)
require(AMORE)
require(affy)
require(robustbase)
require(mclust)
require(tcltk)

op <- par(no.readonly = T)

##########################################################################################

 fileName <- tclvalue(tkgetOpenFile()) 									# Open current the directory to select the file to load.
 if (!nchar(fileName)) {
     tkmessageBox(message = "No file was selected!")
 } else {
     tkmessageBox(message = paste("The file selected was", fileName))
}

# head(fileName)	# return le filename as string

fnames <- unlist(strsplit(fileName, "/"))
nnames <- length(fnames)

path <- fnames[1]														# create the path
for(i in 2:(nnames-1))
	path <- paste(path, fnames[i], sep = "/")
Folder <- fnames[nnames-1]
FileName <- fnames[nnames]

setwd(paste(path, "",sep = "/"))

##########################################################################################

	# Step1: Load the first block of data to collect array information
cat("Step1: Collecte des informations mArray:")
a.info <- try(read.csv(paste(FileName, sep = ""), header = F, fill = T, skip = 1, nrows = 8, sep = "\t"), silent = T)

ArraySet = "Agilent-022060 SurePrint G3 Human CGH Microarray 4x180K"

tmpNames <- as.vector(a.info[1,])
i <- which(tmpNames == "Protocol_Name"); Protocol = as.character(a.info[2, i])
i <- which(tmpNames == "Scan_ScannerName"); ScanName = as.character(a.info[2, i])
i <- which(tmpNames == "Scan_Date"); ScanDate = as.character(a.info[2, i])

ScanDate <- unlist(strsplit(ScanDate, " "))[1]
i <- which(tmpNames == "Grid_Name"); GridName = as.character(a.info[2, i])
i <- which(tmpNames == "FeatureExtractor_Barcode"); BarCode = as.character(a.info[2, i])
i <- which(tmpNames == "FeatureExtractor_UserName"); LabId = as.character(a.info[2, i])
i <- which(tmpNames == "GridPlacement_Version"); FEversion = as.character(a.info[2, i])

analysisDate <- format(Sys.Date(), "%m-%d-%Y")
Safir.Id <- unlist(strsplit(FileName, "_"))[1]

{
if(class(a.info)!="try-error") cat("\nCGH data have been correctly downloaded:", "\n\tFolder: ", getwd(), "\n\tFileName:", FileName, "\n\tSafir.Id:", Safir.Id, "\n")
else stop("\t",FileName, ": No such file on directory (glup!). Please fix it!\n")
}

	# Step2: Load the second block of data to collect Intensity values and QC columns
cat("\nStep2: Chargement des données mArray:")
cgh <- read.csv(paste(FileName, sep = ""), header = T, skip = 9, sep = "\t")

keepCol <- which(as.character(colnames(cgh)) %in% 
		c(	"ProbeName", "SystematicName",
			"gMedianSignal", "rMedianSignal",
			"gIsSaturated", "rIsSaturated",
			"gIsFeatNonUnifOL", "rIsFeatNonUnifOL",
			"gIsWellAboveBG", "rIsWellAboveBG"))


cgh <- cgh[which(substr(cgh$ProbeName, 1, 1)== "A"), keepCol]		# Select QC columns and rows containing A_xxx probes
cat("\nCGH data have been correctly downloaded\n")

		# Step3: Load Agilent annotation table(+ hg19 génomic position + CG% information)
cat("\nStep3: Chargement des tables d'annotation complémentaires:")
agil.hg19 <- try(read.table(paste(Root, ":/Projet Safir/Data SAFIR/Arrays Infos/022060_4x180K_hg19_20110628_GenomePos&GC_FC.txt", sep = ""), header = T, sep = "\t"), silent = T)
{
if(class(agil.hg19)!="try-error") cat("\nAgilent-hg19 information file have been correctly downloaded\n")
else stop("\tNo Agilent hg19 on this directory (glup!). Please fix it!\n")
}

agil.hg19 <- agil.hg19[order(agil.hg19$ProbeName),]


##########################################################################################

	# Step4: Flags suppression
pflags <- NA
if(Flags){
	cat("\nStep4: suppression des flags:\t")
	flags.info <- cgh[,5:10]
	flags <- which(	flags.info[,1]== 1 | 		# 1 = gTsSaturated invalide
			flags.info[,2]== 1 | 				# 1 = rTsSaturated invalide
			flags.info[,3]== 1 | 				# 1 = gIsFeatureNonUnifOL invalide
			flags.info[,4]== 1 | 				# 1 = rIsFeatureNonUnifOL invalide
			flags.info[,5]== 0 | 				# 0 = gIsWellAboveBG invalide
			flags.info[,6]== 0)					# 0 = rIsWellAboveBG invalide

	cgh <- cgh[-flags,]
	pflags <- round(length(flags)/nrow(cgh)*100, 2)
	cat(pflags, "%\n")
	}
	
##########################################################################################

	# Step5: suppression of duplicates probes
cat("\nStep5: Suppression des sondes dupliquées\n")
hg19.duplic <- which(duplicated(agil.hg19$ProbeName))
if(length(hg19.duplic>0)) agil.hg19 <- agil.hg19[-hg19.duplic,]
rownames(agil.hg19) <- seq(1, nrow(agil.hg19))

cgh <- cgh[order(cgh$ProbeName),]
cgh.duplic <- which(duplicated(cgh$ProbeName))
if(length(cgh.duplic>0)) cgh <- cgh[-cgh.duplic,]
rownames(cgh) <- seq(1, nrow(cgh))

##########################################################################################

	# Step6: add genomic positions (according to hg19)
cat("\nStep6: Calcul des positions génomiques selon hg19\n")

hg19 <- agil.hg19[which(is.element(agil.hg19$ProbeName, cgh$ProbeName)), ]
cgh.hg19 <- cgh[which(is.element(cgh$ProbeName, agil.hg19$ProbeName)), ]

verif <- ifelse(as.character(cgh.hg19$ProbeName) == as.character(hg19$ProbeName), "ok", "error")
n.error <- length(which(verif == "error"))

{
if(n.error == 0) {
	cgh <- cbind.data.frame(cgh.hg19[,1:2], hg19, cgh.hg19[,-c(1:2)])
	cat("\n\nHow lucky you are :-)\n\n\n")
	}
else stop("\tStep6: Errors on hg19 annotations")
}
colnames(cgh)[4:6] <- c("ChrNum", "ChrStart", "ChrEnd")
cgh <- cgh[order(cgh$ChrNum, cgh$ChrStart),]

cat("\nSafir.Id:", Safir.Id, "\nBarCode:", BarCode, "\nArraySet:", ArraySet, "\nProtocol:", Protocol, "\nScanName:", ScanName, "\nScanDate:", ScanDate, "\nLabId:", LabId, "\nFEversion:", FEversion, "\nAnalysis Date:", analysisDate, "\n")
return(list(cgh = cgh, Safir.Id = Safir.Id, BarCode = BarCode, ArraySet = ArraySet, Protocol = Protocol, ScanName = ScanName, ScanDate = ScanDate, LabId = LabId, FEversion = FEversion, analysisDate = analysisDate, flags = pflags))
}